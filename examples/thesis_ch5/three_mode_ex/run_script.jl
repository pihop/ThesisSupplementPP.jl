import Random
scaling = parse(Int,ARGS[1])
end_time = parse(Float64,ARGS[2])
ssa_samples = parse(Int,ARGS[3])
iters = parse(Int,ARGS[4])
#ssa_samples = 5000
tspan = (0., end_time)
sampling = 200
tplot = range(0.0, stop=tspan[2], length=sampling) 

Random.seed!(1234)

println("Experiment set for population $scaling")

include("rn.jl")
include("../ssa.jl")
include("hybrid_fluid.jl")
include("hybrid_lna.jl")
include("hybrid_mc.jl")


function run_all()
    println("Population size $scaling")
    for i in 1:iters
        println("========= Iterations $i out of $iters ===============")
         
        u0_ssa = zeros(length(m_grid)+1)
        u0_ssa[CartesianIndex(2,1)] = scaling

        # Fluid inits 
        u0_coupl = zeros(length(m_grid)+3)
        u0_coupl[CartesianIndex(2,1)] = 1.0
        u0_coupl[length(m_grid)+1] = 1.0

        u0_iter = zeros(length(m_grid)+2)
        u0_iter[CartesianIndex(2,1)] = 1.0

        # LNA inits
        n_lna_covars = Int(length(m_grid)*(length(m_grid)+1)/2)
        u0_lna_coupl = zeros(length(m_grid) + n_lna_covars + 3)
        u0_lna_coupl[CartesianIndex(2,1)] = 1.0
        u0_lna_coupl[length(m_grid) + n_lna_covars + 1] = 1.0

        u0_lna_iter = zeros(length(m_grid) + n_lna_covars + 2)
        u0_lna_iter[CartesianIndex(2,1)] = 1.0

        # MC inits
        u0_mc_coupl = zeros(length(lhs_mc_coupl))
        u0_mc_coupl[CartesianIndex(2,1)] = scaling 
        u0_mc_coupl[end-2] = 1.0

        u0_mc_iter = zeros(length(lhs_mc_iter))
        u0_mc_iter[CartesianIndex(2,1)] = scaling
        u0_mc_iter[end-1] = 1.0

#        param = rand(2) 
#        run_ssa(i,5000, u0_ssa, param)

        str = string(i, "_", Int(scaling))
        path = joinpath(@__DIR__, "data/ssa_res_$str.jld2")
        file = jldopen(path, "r")
        param = file["params"] 
#
#        run_fluid(i, u0_coupl, u0_iter, param)
#        run_lna(i, u0_lna_coupl, u0_lna_iter, param)
        run_mc(i, u0_mc_coupl, u0_mc_iter, param)
        println("=====================================================")
    end
end
 
run_all()
