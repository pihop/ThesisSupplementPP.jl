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
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)
pi_3 = symbols(:pi_3)

rns = [rn_mode_1, rn_mode_2, rn_mode_1]
covar_syms = gen_covar_syms(rn_mode_1)

# Full symbolic equation

jump_rates = [(1-d_1)*rs*ut1*scaling, (1-d_2)*d_1*ut1*0.1*rs*scaling]
symbolic_iter = sum(map((rn,w) -> w*rn.f_symfuncs, rns, [w_1, w_2, w_3]))
symb_covar_expr = sum(map((rn,w) -> w*gen_noise_covar(rn, covar_syms), rns, [w_1, w_2, w_3]))
symbolic_iter = vcat(symbolic_iter, collect(Iterators.flatten(symb_covar_expr)), jump_rates)

lhs_lna_iter = unique(vcat(rns[1].syms, vec(covar_syms), [:d_1, :d_2]))
odes_iter_lna = sym_to_expr(symbolic_iter, lhs_lna_iter, rn_full.params, [:w_1, :w_2, :w_3])

symbolic_coupl = (rn_mode_1.f_symfuncs*(pi_1) + rn_mode_2.f_symfuncs*(pi_2) + rn_mode_1.f_symfuncs*(pi_3))

symbolic_coupl_cov = (gen_noise_covar(rn_mode_1, covar_syms)*pi_1 + 
                      gen_noise_covar(rn_mode_2, covar_syms)*pi_2 +
                      gen_noise_covar(rn_mode_1, covar_syms)*pi_3)

filtering = [-pi_1*rs*ut1*scaling, pi_1*rs*ut1*scaling - 0.1*pi_2*rs*ut1*scaling, 0.1*pi_2*rs*ut1*scaling]

lhs_lna_coupl = unique(vcat(rns[1].syms, vec(covar_syms), [:pi_1, :pi_2, :pi_3]))
symbolic_coupl = vcat(symbolic_coupl, collect(Iterators.flatten(symbolic_coupl_cov)), filtering)
odes_coupl_lna = sym_to_expr(symbolic_coupl, lhs_lna_coupl, rns[1].params, [])

# Make the ODE function...
@generated function f_lna_iter(du, u, p, t)
    quote
        du .= $(odes_iter_lna)
    end
end

@generated function f_lna_coupl(du, u, p, t)
    quote
        du .= $(odes_coupl_lna)
    end
end

function w0(t)
    return [1.0, 0.0, 0.0]
end
