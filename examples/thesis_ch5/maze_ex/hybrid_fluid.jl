using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using LinearAlgebra
using SymEngine

include("../common_run.jl")
include("../common3.jl")

rs = symbols(:rs)
d_1 = symbols(:d_1)
d_2 = symbols(:d_2)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
w_3 = symbols(:w_3)
ut1 = symbols(t1) 
ut2 = symbols(t2) 
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)
pi_3 = symbols(:pi_3)

rns = [rn_mode_1, rn_mode_2, rn_mode_3]
# Full symbolic equation

symbolic_iter = sum(map((rn,w) -> w*rn.f_symfuncs, rns, [w_1, w_2, w_3]))
# Mode switching dist.
jump_rates = [(1-d_1)*rs*ut1*scaling, (1-d_2)*d_1*rs*ut2*scaling]
symbolic_iter = vcat(symbolic_iter, jump_rates)
odes_iter = sym_to_expr(symbolic_iter, vcat(rns[1].syms, [:d_1, :d_2]), rns[1].params, [:w_1, :w_2, :w_3])
lhs_fluid_iter = vcat(rns[1].syms, [:d_1, :d_2])

symbolic_coupl = (rn_mode_1.f_symfuncs*(pi_1) + rn_mode_2.f_symfuncs*(pi_2) + rn_mode_3.f_symfuncs*(pi_3))
filtering = [-pi_1*rs*ut1*scaling, pi_1*rs*ut1*scaling - pi_2*rs*ut2*scaling, pi_2*rs*ut2*scaling]
symbolic_coupl= vcat(symbolic_coupl, filtering)
odes_coupl = sym_to_expr(symbolic_coupl, vcat(rns[1].syms, [:pi_1, :pi_2, :pi_3]), rns[1].params, [])
lhs_fluid_cpl = vcat(rns[1].syms, [:pi_1, :pi_2, :pi_3])

# Make the ODE function...
@generated function f_iter(du, u, p, t)
    quote
        du .= $(odes_iter)
    end
end

@generated function f_coupl(du, u, p, t)
    quote
        du .= $(odes_coupl)
    end
end

function w0(t)
    return [1.0, 0.0, 0.0]
end
