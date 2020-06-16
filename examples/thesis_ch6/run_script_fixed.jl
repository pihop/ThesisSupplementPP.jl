# NOTE: 
# This is the same as run_script.jl except the 200 tested parameters are sampled
# once and used for all methods (lna, moment closure and ssa).
using JLD2
using DataFrames
using Random

function get_variance(idx, symbs, len, method)
    if method == "lna"
        var = Symbol("c_$(idx)$(idx)")
        return findall(x->x==var, symbs)[1]
    elseif method == "mc"
        array = zeros(Int, len)
        array[idx] = 2 
        str = string(array...)
        var = Symbol("M_$str")
        return findall(x->x==var, symbs)[1]
    else
        return 0
    end
end

scaling = parse(Int,ARGS[1])
sca = Int(scaling)

tspan = (0.0, 10.0)
sampling = 100
# 200 parametrisations
np = 200
iters = 5000
tplot = tspan[1]:tspan[2]/sampling:tspan[2]

include("rn.jl")

# SSA simulation is set up in the following file
include("ssa_simulation.jl")
u0_ssa = [0.0,scaling, 0.0, 0.0, 0.0]

## Fluid approximation related constructions are in the following file  .
## We are not going to run the fluid approximation here as no timing
# comparisions are given in chapter 6. Quantitative results are same are for lna
# below.
#include("hybrid_fluid.jl")
## Set up the initial conditions.
#f_u0_cpl = zeros(length(f_lhs_cpl))
#f_u0_iter = zeros(length(f_lhs_iter))
#f_u0_cpl[2] = 1.0
#f_u0_iter[2] = 1.0
#f_u0_cpl[end-1] = 1.0

# Moment-based constructions are in the following file  .
include("hybrid_mc.jl")
# Set up the initial conditions.
m_u0_iter = zeros(length(m_lhs_iter))
m_u0_cpl = zeros(length(m_lhs_cpl))
m_u0_iter[2] = scaling
m_u0_cpl[2] = scaling
m_u0_cpl[end-1] = 1.0

# Linear noise-based constructions are in the following file  .
include("hybrid_lna.jl")
# Set up the initial conditions.
l_u0_iter = zeros(length(l_lhs_iter))
l_u0_cpl = zeros(length(l_lhs_cpl))
l_u0_iter[2] = 1.0
l_u0_cpl[2] = 1.0
l_u0_cpl[end-1] = 1.0

function eval_ssa(p, u0, objective)
    mean, var = run_ssa_nosave(rn_full, u0, tspan, p, iters, sampling)
    mean = getindex.(mean.u, 4) ./ sum(u0[1:4])
    var = getindex.(var.u, 4) ./ sum(u0[1:4])^2
    return objective(mean, var)
end

function eval_lna(p, u0, sol_method, objective, symbs)
    var_idx = get_variance(4, symbs, 4, "lna")
    sol = sol_method(1, u0, p)
    mean = getindex.(sol.u, 4) 
    var = getindex.(sol.u, var_idx) ./ scaling 
    return objective(mean, var)
end

function eval_mc(p, u0, sol_method, objective, symbs)
    var_idx = get_variance(4, symbs, 4, "mc")
    mean, var = sol_method(1, u0, p)
    mean = getindex.(mean.u, 4) ./ scaling 
    var = getindex.(var.u, var_idx) ./ scaling^2 
    return objective(mean, var)
end

# Objective functions to evaluate. Provides goal_is_sat and prob_goal_is_sat.
include("opt_goals.jl")

Random.seed!(0)
# Fix movement to 1.0 and sensing to 0.1.
pm_fix = fill(1.0, np)
ps_fix = fill(0.1, np)
# randomly sample 200 parameter values for succp.
ps5 = rand(np)
ps7 = rand(np)
ps9 = rand(np)
# three parameters for congestion a.
pa_fix5 = fill(0.5, np)
pa_fix7 = fill(0.7, np)
pa_fix9 = fill(0.9, np)

scaling_raw = fill(scaling, np)
scaling_scaled = fill(1, np)

# movement rate r_m, sense r_s, fault prob p, parms of the congestion/interference
params_raw5 = zip(pm_fix, ps_fix, ps5, pa_fix5, scaling_raw) 
params_raw7 = zip(pm_fix, ps_fix, ps7, pa_fix7, scaling_raw) 
params_raw9 = zip(pm_fix, ps_fix, ps9, pa_fix9, scaling_raw) 
params_raw = vcat(collect(params_raw5), collect(params_raw7), collect(params_raw9))

params_scaled5 = zip(pm_fix, ps_fix, ps5, pa_fix5, scaling_scaled) 
params_scaled7 = zip(pm_fix, ps_fix, ps7, pa_fix7, scaling_scaled) 
params_scaled9 = zip(pm_fix, ps_fix, ps9, pa_fix9, scaling_scaled) 
params_scaled = vcat(collect(params_scaled5), collect(params_scaled7), collect(params_scaled9))

# Run the model then evaluate the objective.
println("First -- Coupling based LNA construction")
t_lcpl = @elapsed data_lcpl = map(x -> [x..., 
                                        eval_lna(x, l_u0_cpl, run_noise_cpl, prob_goal_is_sat, l_lhs_cpl),
                                        eval_lna(x, l_u0_cpl, run_noise_cpl, goal_is_sat, l_lhs_cpl)],
                                  params_scaled)
data_lcpl = convert(DataFrame, hcat(data_lcpl...)') 

println("Time taken $t_lcpl")

println("Second --- Iterative LNA construction")
t_liter = @elapsed data_liter = map(x -> [x..., 
                                          eval_lna(x, l_u0_iter, run_noise_iter, prob_goal_is_sat, l_lhs_iter),
                                          eval_lna(x, l_u0_iter, run_noise_iter, goal_is_sat, l_lhs_iter)],
                                   params_scaled)
data_liter = convert(DataFrame, hcat(data_liter...)') 
println("Time taken $t_liter")


println("Third --- Coupling based MC construction")
t_mcpl = @elapsed data_mcpl = map(x -> [x..., 
                                        eval_mc(x, m_u0_cpl, run_mc_cpl, prob_goal_is_sat, m_lhs_cpl),
                                        eval_mc(x, m_u0_cpl, run_mc_cpl, goal_is_sat, m_lhs_cpl)],
                                        params_raw)
data_mcpl = convert(DataFrame, hcat(data_mcpl...)') 
println("Time taken $t_mcpl")

println("Fourth --- Iterative LNA construction")
t_miter = @elapsed data_miter = map(x -> [x..., 
                                          eval_mc(x, m_u0_iter, run_mc_iter, prob_goal_is_sat, m_lhs_iter),
                                          eval_mc(x, m_u0_iter, run_mc_iter, goal_is_sat, m_lhs_iter)],
                                          params_raw)
data_miter = convert(DataFrame, hcat(data_miter...)') 
println("Time taken $t_miter")

println("Fifth --- Direct SSA")
t_ssa = @elapsed data_ssa = map(x -> [x..., 
                                      eval_ssa(x, u0_ssa, prob_goal_is_sat),
                                      eval_ssa(x, u0_ssa, goal_is_sat)],
                                      params_raw)
data_ssa = convert(DataFrame, hcat(data_ssa...)') 

println("Time taken $t_ssa")

# Write all the data to a file.
jldopen(joinpath("data/data_f$sca.jld2"), "w") do file
    write(file, "lcpl/init", l_u0_cpl)
    write(file, "lcpl/data", data_lcpl)
    write(file, "lcpl/time", t_lcpl)
    write(file, "liter/init", l_u0_iter)
    write(file, "liter/data", data_liter)
    write(file, "liter/time", t_liter)
    write(file, "mcpl/init", m_u0_cpl)
    write(file, "mcpl/data", data_mcpl)
    write(file, "mcpl/time", t_mcpl)
    write(file, "miter/init", m_u0_iter)
    write(file, "miter/data", data_miter)
    write(file, "miter/time", t_miter)
    write(file, "ssa/init", u0_ssa)
    write(file, "ssa/data", data_ssa)
    write(file, "ssa/time", t_ssa)
    write(file, "params/scaling", scaling)
    write(file, "params/tspan", tspan)
    write(file, "params/sampling", sampling)
    write(file, "params/iters", iters)
end
