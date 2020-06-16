using JLD2
using DataFrames
using StatsPlots
using LaTeXStrings
using Plots

pyplot()
PyPlot.rc("text", usetex="true")

include("rn.jl")
include("ssa_simulation.jl")
include("opt_goals.jl")

sampling = 100
iters = 5000

function eval_ssa(par, p, u0, objective)
    par = par
    par[3] = p[1]
    mean, var = run_ssa_nosave(rn_full, 1, u0, tspan, par, iters, sampling)
    mean = getindex.(mean.u, 4) ./ sum(u0[1:4])
    var = getindex.(var.u, 4) ./ sum(u0[1:4])^2
    return objective(mean, var, p)
end

# Cost function for direct optimisation.
function discounted(mean, var, p)
    if goal_is_sat(mean, var)
        return -p[1]
    else 
        return p[1] 
    end
end

# Initial conditions.
tspan = (0.0, 10.0)
scaling = 100
u0_ssa = [0.0,scaling, 0.0, 0.0, 0.0]
params_raw = [1.0, 0.1, 1, 0.01, scaling]

ps = rand(200)
# This might take some time. ~16 min on my laptop.
t_ssa = @elapsed data_ssa = DataFrame(X = ps,
                                      Y = map(x -> eval_ssa(params_raw, [x], u0_ssa, discounted), ps))

plt = scatter(ylabel="Reward", xlabel="succp", legend=:right, fg_legend=:transparent,
              bg_legend=:transparent,
              size=(700,400),
              xtickfont=font(18), 
              ytickfont=font(18),
              guidefont=font(18), 
              legendfont=font(18))

plt = scatter!(data_ssa[:X], data_ssa[:Y]; labels="", marker = (stroke(0, 0.0, :black, :dot)))
savefig("plots/direct_reward.pdf")




