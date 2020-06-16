using Plots
pyplot()

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


scaling = 200.0
# Movement rate r_m, sense r_s, fault prob p, parms of the congestion/interference
p = 0.3
params_raw = [1.0, 0.1, p, 0.01, scaling]
params_scaled = [1.0, 0.1, p, 0.01, 1]
u0_ssa = [0.0,scaling, 0.0, 0.0, 0.0]
tspan = (0.0, 10.0)
sampling = 100
iters = 5000
tplot = tspan[1]:tspan[2]/sampling:tspan[2]

include("rn.jl")
include("ssa_simulation.jl")
include("hybrid_mc.jl")
include("hybrid_lna.jl")
include("opt_goals.jl")

# Initial conditions for the moment based approximation.
m_u0_iter = zeros(length(m_lhs_iter))
m_u0_cpl = zeros(length(m_lhs_cpl))
m_u0_iter[2] = scaling
m_u0_cpl[2] = scaling
m_u0_cpl[end-1] = 1.0

# Initial condition for the linear noise based approximation.
l_u0_iter = zeros(length(l_lhs_iter))
l_u0_cpl = zeros(length(l_lhs_cpl))
l_u0_iter[2] = 1.0
l_u0_cpl[2] = 1.0
l_u0_cpl[end-1] = 1.0

s_sol, s_sums = run_ssa(rn_full, 1, u0_ssa, tspan, params_raw, iters, sampling)
l_sol_cpl, l_sol_iter = run_noise(1, l_u0_cpl, l_u0_iter, params_scaled)
m_sol_cpl, m_sol_iter = run_mc(1, m_u0_cpl, m_u0_iter, params_raw)

var_idx_l = get_variance(4, l_lhs_cpl, 4, "lna")
var_idx_m = get_variance(4, m_lhs_cpl, 4, "mc")

plt_dist = plot()
var_l_cpl = sqrt.(getindex.(l_sol_cpl.u, var_idx_l) ./ scaling) 
mean_l_cpl = getindex.(l_sol_cpl.u, 4)
var_l_iter = sqrt.(getindex.(l_sol_iter[2][2].u, var_idx_l) ./ scaling) 
mean_l_iter = getindex.(l_sol_iter[2][2].u, 4)
var_m_cpl = sqrt.(getindex.(m_sol_cpl.u, var_idx_m)) ./ scaling
mean_m_cpl = getindex.(m_sol_cpl.u, 4) ./ scaling
var_m_iter = sqrt.(getindex.(m_sol_iter[2][2].u, var_idx_m))  ./ scaling
mean_m_iter = getindex.(m_sol_iter[2][2].u, 4) ./ scaling
var_ssa = sqrt.(getindex.(s_sums.v.u, 4)) ./ scaling
mean_ssa = getindex.(s_sums.u.u, 4) ./ scaling
plt = plot!(tplot, mean_l_cpl + var_l_cpl; color=1, line=(:dot))
plt = plot!(tplot, mean_l_cpl - var_l_cpl; color=1, line=(:dot))
plt = plot!(tplot, mean_l_iter + var_l_iter; color=5, line=(:dot))
plt = plot!(tplot, mean_l_iter - var_l_iter; color=5, line=(:dot))
plt = plot!(tplot, mean_ssa - var_ssa; color=2)
plt = plot!(tplot, mean_ssa + var_ssa; color=2)
plt = plot!(tplot, mean_m_cpl - var_m_cpl; color=1, line=(:dash,))
plt = plot!(tplot, mean_m_cpl + var_m_cpl; color=1, line=(:dash,))
plt = plot!(tplot, mean_m_iter - var_m_iter; color=3, line=(:dashdot,))
plt = plot!(tplot, mean_m_iter + var_m_iter; color=3, line=(:dashdot,))
