using DiffEqBiological
using ThesisSupplementPP
using DiffEqBase, DiffEqJump, OrdinaryDiffEq
using Interpolations
using JLD2

include("parametrisation.jl")

function f(du,u,p,t)
    du[:] .= 0.
end

# As described in the thesis we use the empirical results from the stochastic
# simulation to get the weight function.
# This is just to demonstrate that the construction is "correct".
@load "data/ssa.jld2" sim_sum
@time spl = LinearInterpolation(sim_sum.t, getindex.(sim_sum.u.u, 5))
@reaction_func weight(t) = spl(t) 

rn_i = @reaction_network rnType begin
    (1-weight(t))*rm, x_1 --> x_2
    (1-weight(t))*0.5*rm, x_2 --> x_3
    (1-weight(t))*0.5*rm, x_3 --> x_4
    (1-weight(t))*rm, x_4 --> x_3
    (1-weight(t))*0.5*rm, x_3 --> x_2
    (1-weight(t))*0.5*rm, x_2 --> x_1

    weight(t)*rm, x_1 --> x_2
    weight(t)*rm, x_2 --> x_3
    weight(t)*rm, x_3 --> x_4
end rm rs

dprob = ODEProblem(f, u0_inhom, tspan, p)
jprob = JumpProblem(dprob, Direct(), jumps(rn_i)...; dep_graph=rxtospecies_depgraph(rn_i))
eprob = EnsembleProblem(jprob)

println("Monte carlo simulation of the model...")
sim_time = @elapsed sim_sol = solve(eprob, Tsit5(), EnsembleThreads(); 
                                        trajectories=1000, progress=true)
println("Time taken $sim_time s"...)
println("Constructing simulation summary...")
@time sim_sum = EnsembleSummary(sim_sol, tplot)

jldopen("data/ssa_inhom.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "initial", u0_inhom)
    write(file, "params", p)
    write(file, "sim_time", sim_time)
    write(file, "sim_sum", sim_sum)
end
