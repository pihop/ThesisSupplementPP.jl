using DiffEqBiological
using DiffEqJump 
using OrdinaryDiffEq
using DiffEqBase.EnsembleAnalysis
using JLD2

include(joinpath(@__DIR__, "rn.jl"))

dprob = DiscreteProblem(rn_mode_1, u0_ssa, tspan, p)
jprob = JumpProblem(dprob, DirectFW(), rn_mode_1)
println("Monte carlo simulation of the model...")
sim_time = @elapsed sim_sol = solve(EnsembleProblem(jprob), FunctionMap(), EnsembleThreads(); trajectories=5000)
println("Time taken $sim_time s"...)
println("Constructing simulation summary...")
@time sim_sum = EnsembleSummary(sim_sol, tplot)

jldopen(joinpath(@__DIR__, "data/ssa_m1_res.jld2"), "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "initial", u0_ssa)
    write(file, "params", p)
    write(file, "sim_time", sim_time)
    write(file, "sim_sum", sim_sum)
end
