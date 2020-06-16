using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using LinearAlgebra
using SymEngine

include("../common_run.jl")
include("../common3.jl")

# The moment equations are generated separately by mc_eqns.jl. 
file = jldopen("eqns/mc_eqns.jld2", "r") 

mc_modes = file["symb_expr"]
mc_modes = [convert.(Basic, eqn) for eqn in mc_modes]
syms = convert.(Basic, file["symbs"])

jump_rates = [(1-d_1)*rs*ut1*scaling, (1-d_2)*d_1*ut1*0.1*rs*scaling]
symbolic_iter_mc = sum(map((mc,w) -> w * mc, mc_modes, [w_1, w_2, w_3]))

symbolic_iter_mc = vcat(symbolic_iter_mc, jump_rates)
lhs_mc_iter = vcat(Symbol.(vec(syms)), [:d_1, :d_2])
mc_odes_iter = sym_to_expr(symbolic_iter_mc, lhs_mc_iter, rn_full.params, [:w_1, :w_2, :w_3])

symbolic_coupl_mc = (mc_modes[1]*pi_1 + mc_modes[2]*pi_2 + mc_modes[3]*pi_3)

filtering = [-pi_1*rs*ut1*scaling, 
             pi_1*rs*ut1*scaling - 0.1*pi_2*rs*ut1*scaling, 
             0.1*pi_2*rs*ut1*scaling]

symbolic_coupl_mc = vcat(symbolic_coupl_mc, filtering)
lhs_mc_coupl = vcat(Symbol.(vec(syms)), [:pi_1, :pi_2, :pi_3])
mc_odes_coupl = sym_to_expr(symbolic_coupl_mc, lhs_mc_coupl, rn_full.params, [])

# Make the ODE function...
@generated function f_mc_iter(du, u, p, t)
    quote
        du .= $(mc_odes_iter)
    end
end

@generated function f_mc_coupl(du, u, p, t)
    quote
        du .= $(mc_odes_coupl)
    end
end

function w0(t)
    return [1.0, 0.0, 0.0]
end
