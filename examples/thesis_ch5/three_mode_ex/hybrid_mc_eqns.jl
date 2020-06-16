using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using LinearAlgebra
using SymEngine

include("rn.jl")
include("../common_run.jl")
include("../common3.jl")

rs = symbols(:rs)
d_1 = symbols(:d_1)
d_2 = symbols(:d_2)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
w_3 = symbols(:w_3)
ut1 = symbols(t1) 
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)
pi_3 = symbols(:pi_3)

rns = [rn_mode_1, rn_mode_2, rn_mode_1]

mc_modes = [make_symbolic_mc(rn,2) for rn in rns]
syms = mc_modes[1][1]
to_save = [convert.(Expr, eqn) for eqn in getindex.(mc_modes, 2)]

jldopen("eqns/mc_eqns.jld2", "w") do file
    write(file, "symb_expr", to_save)
    write(file, "symbs", convert.(Expr, syms))
end


