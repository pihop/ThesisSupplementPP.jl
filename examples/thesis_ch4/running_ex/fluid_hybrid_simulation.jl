# Construct and simulate the hybrid fluid process. Running example of
# Chapter 4.
#
using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using JLD2
using LinearAlgebra
using LaTeXStrings
using Sundials
using SymEngine

include("common.jl")
include("rn.jl")

scaling = parse(Int, ARGS[1])
trajectories = parse(Int, ARGS[2])
p = (1.0, 0.2)
tspan = (0., 10.)
tplot = 0.0:0.01:10.0
u0 = [0., 1., 0., 0., 0.]
u0_ssa = [0., scaling, 0., 0., 0.]

rns = [rn_mode_1, rn_mode_2]

#fluid_odes = [rn.f_symfuncs for rn in rns] 
fluid_odes = [rn.f_symfuncs for rn in rns]
lhs = rn_mode_1.syms

fluid_odes_expr = [sym_to_expr(ode, lhs, rn_full.params, []) for ode in fluid_odes]

# Make the ODE function...
@generated function f_fluid_hyb(du, u, p, t)
    quote
        if p[2] == 1
            du .= $(fluid_odes_expr[1])
        else
            du .= $(fluid_odes_expr[2])
        end
    end
end

rate_switch(u,p,t) = p[1][2]*u[4]*scaling

function affect_switch!(integrator)
    integrator.p = (p, 0)
end

u0_fluid = zeros(length(lhs))
u0_fluid[2] = 1.0

jump_switch = VariableRateJump(rate_switch,affect_switch!)
oprob = ODEProblem(f_fluid_hyb, u0_fluid, tspan, (p, 1))
jprob = JumpProblem(oprob, Direct(), jump_switch)
eprob = EnsembleProblem(jprob)
println("The monte carlo simulation")
sim_time = @elapsed sol = solve(eprob, Tsit5(), EnsembleThreads();  trajectories=trajectories)
println("Time taken $sim_time")
sumr = EnsembleSummary(sol, tplot)

jldopen("data/hyridb_fluid_$scaling.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "initial", u0_fluid)
    write(file, "params", p)
    write(file, "sim/time", sim_time)
    # Only save the summary rather than all trajectories.
    # write(file, "sim/sol", sol)
    write(file, "sim/sum", sumr)
end
