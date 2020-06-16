using DiffEqBiological
using DiffEqJump 
using OrdinaryDiffEq
using DiffEqBase.EnsembleAnalysis
using JLD2

function run_ssa(iter, traj, u0, param)
    dprob = DiscreteProblem(rn_full, u0, tspan, param)
    jprob = JumpProblem(dprob, DirectFW(), rn_full)
    println("Monte carlo simulation of the model...")
    sim_time = @elapsed sim_sol = solve(EnsembleProblem(jprob), FunctionMap(), EnsembleThreads(); trajectories=traj)
    println("Time taken $sim_time s"...)
    println("Constructing simulation summary...")
    sum_time = @elapsed sim_sum = EnsembleSummary(sim_sol, tplot)
    println("Time taken $sum_time s"...)

    str = string(iter, "_", Int(scaling))

    jldopen(joinpath("data/ssa_res_$str.jld2"), "w") do file
        write(file, "tspan", tspan)
        write(file, "sampling", sampling)
        write(file, "pop", scaling)
        write(file, "initial", u0)
        write(file, "params", param)
        write(file, "sim_time", sim_time)
        write(file, "total_time", sim_time+sum_time)
        write(file, "sim_sum", sim_sum)
    end
#    return sim_sol
end
