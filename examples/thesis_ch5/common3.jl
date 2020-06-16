using LinearAlgebra 
using DiffEqBase, OrdinaryDiffEq
using Interpolations
using Sundials

function run_single(func, u0, tspan, p, weights, saveat)
    # Runs a single solution given a weight.
    oprob = ODEProblem(func, u0, tspan, (p,weights))
    return solve(oprob, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, saveat=saveat)
end

function w_a(f, t)
    [1-f[1](t), f[1](t) - f[2](t), f[2](t)]
end

function run_iterative(func, u0, tspan, p, n_eqn, n_modes, weight_0, saveat)
    hs = []
    sols = []
    push!(sols, run_single(func, u0, tspan, p, weight_0, saveat))
    push!(hs, LinearInterpolation(tplot, [sols[1](t)[n_eqn-(n_modes-1)+1] for t in tplot]))
    push!(sols, run_single(func, u0, tspan, p, t -> w_a([t->hs[1](t), t -> 0.0], t), saveat))
    push!(hs, LinearInterpolation(tplot, [sols[2](t)[n_eqn-(n_modes-1)+2] for t in tplot]))
    push!(sols, run_single(func, u0, tspan, p, t -> w_a([t->hs[1](t), t->hs[2](t)], t), saveat))
    return (hs, sols)
end

function f(du,u,p,t)
    du[:] .= 0.
end

