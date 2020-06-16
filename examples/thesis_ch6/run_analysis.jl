using Plots
using GLM
using JLD2
using DataFrames
using StatsPlots
pyplot()

function rel_err(x, y)
    return abs((x-y))/x
end

filename = ARGS[1]
#filename = "data/data_f100.jld2"

file = jldopen(filename, "r")

for param in [0.5, 0.7, 0.9]
    df_lna_cpl = file["lcpl/data"]
    df_lna_iter = file["liter/data"] 
    df_mc_cpl = file["mcpl/data"] 
    df_mc_iter = file["miter/data"] 
    df_ssa = file["ssa/data"] 

    df_lna_cpl = df_lna_cpl[df_lna_cpl.x4 .== param, :]
    df_lna_iter = df_lna_iter[df_lna_iter.x4 .== param, :] 
    df_mc_cpl = df_mc_cpl[df_mc_cpl.x4 .== param, :]  
    df_mc_iter = df_mc_iter[df_mc_iter.x4 .== param, :]   
    df_ssa = df_ssa[df_ssa.x4 .== param, :]    

    # First objective
    logist_l_cpl = glm(@formula(x6 ~ x3), df_lna_cpl, Binomial(), LogitLink())
    logist_l_iter = glm(@formula(x6 ~ x3), df_lna_iter, Binomial(), LogitLink())
    logist_m_cpl = glm(@formula(x6 ~ x3), df_mc_cpl, Binomial(), LogitLink())
    logist_m_iter = glm(@formula(x6 ~ x3), df_mc_iter, Binomial(), LogitLink())
    logist_ssa = glm(@formula(x6 ~ x3), df_ssa, Binomial(), LogitLink())

    # Second objective
    mlogist_l_cpl = glm(@formula(x7 ~ x3), df_lna_cpl , Binomial(), LogitLink())
    mlogist_l_iter = glm(@formula(x7 ~ x3), df_lna_iter, Binomial(), LogitLink())
    mlogist_m_cpl = glm(@formula(x7 ~ x3), df_mc_cpl, Binomial(), LogitLink())
    mlogist_m_iter = glm(@formula(x7 ~ x3), df_mc_iter, Binomial(), LogitLink())
    mlogist_ssa = glm(@formula(x7 ~ x3), df_ssa, Binomial(), LogitLink())

    # Decision bounaries. First objecttive.
    mb_ssa = -mlogist_ssa.model.pp.beta0[1] / mlogist_ssa.model.pp.beta0[2]
    mb_lcpl = -mlogist_l_cpl.model.pp.beta0[1] / mlogist_l_cpl.model.pp.beta0[2]
    mb_liter = -mlogist_l_iter.model.pp.beta0[1] / mlogist_l_iter.model.pp.beta0[2]
    mb_mcpl = -mlogist_m_cpl.model.pp.beta0[1] / mlogist_m_cpl.model.pp.beta0[2]
    mb_miter = -mlogist_m_iter.model.pp.beta0[1] / mlogist_m_iter.model.pp.beta0[2]

    # Decision bounaries. Second objecttive.
    b_ssa = -logist_ssa.model.pp.beta0[1] / logist_ssa.model.pp.beta0[2]
    b_lcpl = -logist_l_cpl.model.pp.beta0[1] / logist_l_cpl.model.pp.beta0[2]
    b_liter = -logist_l_iter.model.pp.beta0[1] / logist_l_iter.model.pp.beta0[2]
    b_mcpl = -logist_m_cpl.model.pp.beta0[1] / logist_m_cpl.model.pp.beta0[2]
    b_miter = -logist_m_iter.model.pp.beta0[1] / logist_m_iter.model.pp.beta0[2]


    println("Results for congestion parameter $param")
    println("Decision boundaries. Obj_1")
    println("Noise coupling \t\t Noise iterative \t MC coupling \t\t MC iterative \t\t SSA")
    println(mb_lcpl, "\t", mb_liter, "\t", mb_mcpl, "\t", mb_miter, "\t", mb_ssa)

    println("Decision boundaries. Obj_2")
    println("Noise coupling \t\t Noise iterative \t MC coupling \t\t MC iterative \t\t SSA")
    println(b_lcpl, "\t", b_liter, "\t", b_mcpl, "\t", b_miter, "\t", b_ssa)

    println("Relative errors. Obj_1")
    println("Noise coupling \t\t Noise iterative \t MC coupling \t\t MC iterative \t\t SSA")
    println(rel_err(mb_ssa, mb_lcpl), "\t", rel_err(mb_ssa, mb_liter), "\t", 
            rel_err(mb_ssa, mb_mcpl), "\t", rel_err(mb_ssa, mb_miter), "\t", 
            rel_err(mb_ssa, mb_ssa))

    println("Relative errors. Obj_2")
    println("Noise coupling \t\t Noise iterative \t MC coupling \t\t MC iterative \t\t SSA")
    println(rel_err(b_ssa, b_lcpl), "\t", rel_err(b_ssa, b_liter), "\t\t\t",
            rel_err(b_ssa, b_mcpl), "\t", rel_err(b_ssa, b_miter), "\t\t\t", 
            rel_err(b_ssa, b_ssa))
    println()
end
