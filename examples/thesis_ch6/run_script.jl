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
# Movement rate r_m, sense r_s, fault prob p, parms of the congestion/interference
params_raw = [1.0, 0.1, 1, 0.01, scaling]
params_scaled = [1.0, 0.1, 1, 0.01, 1]

tspan = (0.0, 10.0)
sampling = 100
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

function eval_ssa(par, p, u0, objective)
    par = par
    par[3] = p[1]
    mean, var = run_ssa_nosave(rn_full, 1, u0, tspan, par, iters, sampling)
    mean = getindex.(mean.u, 4) ./ sum(u0[1:4])
    var = getindex.(var.u, 4) ./ sum(u0[1:4])^2
    return objective(mean, var)
end

function eval_lna(par, p, u0, sol_method, objective, symbs)
    par = par
    par[3] = p[1]
    var_idx = get_variance(4, symbs, 4, "lna")
    sol = sol_method(1, u0, par)
    mean = getindex.(sol.u, 4) 
    var = getindex.(sol.u, var_idx) ./ scaling 
    return objective(mean, var)
end

function eval_mc(par, p, u0, sol_method, objective, symbs)
    par = par
    par[3] = p[1]
    var_idx = get_variance(4, symbs, 4, "mc")
    mean, var = sol_method(1, u0, par)
    mean = getindex.(mean.u, 4) ./ scaling 
    var = getindex.(var.u, var_idx) ./ scaling^2 
    return objective(mean, var)
end

# Objective functions to evaluate. Provides goal_is_sat and prob_goal_is_sat.
include("opt_goals.jl")

# Randomly sample 200 parameter values.
Random.seed!(0)
ps = rand(200)

# Run the model then evaluate the objective.
println("First")
t_lcpl = @elapsed data_lcpl = DataFrame(X = ps, 
                                        Y = map(x -> eval_lna(params_scaled, [x], l_u0_cpl, 
                                                              run_noise_cpl, 
                                                              prob_goal_is_sat, l_lhs_cpl), ps),
                                        Z = map(x -> eval_lna(params_scaled, [x], l_u0_cpl, 
                                                              run_noise_cpl, goal_is_sat, l_lhs_cpl), ps))
println("Time taken $t_lcpl")

ps = rand(200)
println("Second")
t_liter = @elapsed data_liter = DataFrame(X = ps, 
                                         Y = map(x -> eval_lna(params_scaled, [x], l_u0_iter, 
                                                               run_noise_iter,
                                                               prob_goal_is_sat, l_lhs_iter), ps),
                                         Z = map(x -> eval_lna(params_scaled, [x], l_u0_iter, run_noise_iter,
                                                               goal_is_sat, l_lhs_iter), ps))
println("Time taken $t_liter")

ps = rand(200)
println("Third")
t_mcpl = @elapsed data_mcpl = DataFrame(X = ps, 
                                        Y = map(x -> eval_mc(params_raw, [x], m_u0_cpl, 
                                                             run_mc_cpl, 
                                                             prob_goal_is_sat, m_lhs_cpl), ps),
                                        Z = map(x -> eval_mc(params_raw, [x], m_u0_cpl, run_mc_cpl, 
                                                             goal_is_sat, m_lhs_cpl), ps))
println("Time taken $t_mcpl")

ps = rand(200)
println("Fourth")
t_miter = @elapsed data_miter = DataFrame(X = ps, 
                                          Y = map(x -> eval_mc(params_raw, [x], m_u0_iter, 
                                                               run_mc_iter, 
                                                               prob_goal_is_sat, m_lhs_iter), ps),
                                          Z = map(x -> eval_mc(params_raw, [x], m_u0_iter, 
                                                               run_mc_iter, 
                                                               goal_is_sat, m_lhs_iter), ps))
println("Time taken $t_miter")

ps = rand(200)
println("Fifth")
t_ssa = @elapsed data_ssa = DataFrame(X = ps,
                                        Y = map(x -> eval_ssa(params_raw, [x], u0_ssa, prob_goal_is_sat), ps),
                                        Z = map(x -> eval_ssa(params_raw, [x], u0_ssa, goal_is_sat), ps))
println("Time taken $t_ssa")

# Write all the data to a file.
jldopen(joinpath("data/data_$sca.jld2"), "w") do file
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
