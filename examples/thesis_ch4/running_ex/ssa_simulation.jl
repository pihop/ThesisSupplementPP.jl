using DiffEqBiological
using DiffEqJump 
using OrdinaryDiffEq
using DiffEqBase.EnsembleAnalysis
using JLD2

include("rn.jl")

scaling = parse(Int, ARGS[1]) 
trajectories = parse(Int, ARGS[2])
tspan = (0., 10.)
tplot = 0.0:0.01:10.0
#scaling = 100
p = (1.0, 0.2)
u0_ssa = [0., scaling, 0., 0., 0.]

dprob = DiscreteProblem(rn_full, u0_ssa, tspan, p)
jprob = JumpProblem(dprob, DirectFW(), rn_full)
println("Monte carlo simulation of the model...")
sim_time = @elapsed sim_sol = solve(EnsembleProblem(jprob), FunctionMap(), EnsembleThreads(); trajectories = trajectories)
println("Time taken $sim_time s"...)
println("Constructing simulation summary...")
@time sim_sum = EnsembleSummary(sim_sol, tplot)

jldopen("data/ssa_$scaling.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "initial", u0_ssa)
    write(file, "params", p)
    write(file, "sim_time", sim_time)
    write(file, "sim_sum", sim_sum)
end
