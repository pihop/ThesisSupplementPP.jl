using LinearAlgebra 
using DiffEqBase, OrdinaryDiffEq
using Interpolations

function run_single(func, u0, tspan, p, weights, saveat)
    # Runs a single solution given a weight.
    oprob = ODEProblem(func, u0, tspan, (p,weights))
    return solve(oprob, Tsit5(), saveat=saveat)
end

function run_iterative(func, u0, tspan, p, n_eqn, n_modes, weight_0, saveat)
    hs = []
    sols = [] 
    push!(sols, run_single(func, u0, tspan, p, weight_0, saveat))
    push!(hs, LinearInterpolation(tplot, [sols[1](t)[n_eqn-(n_modes-1)+1] for t in tplot]))
    @inline function w_1(t) 
        return [1-hs[1](t),hs[1](t)]
    end
    push!(sols, run_single(func, u0, tspan, p, w_1, saveat))
    return hs, sols
end

function f(du,u,p,t)
    du[:] .= 0.
end

