using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using LinearAlgebra
using SymEngine

include("common.jl")

# Setting up the symbols.
rs = symbols(:rs)
d_1 = symbols(:d_1)
d_2 = symbols(:d_2)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)
pi_3 = symbols(:pi_3)
x_4 = symbols(:x_4)

l_rns = [rn_mode_1, rn_mode_2]
l_cov_syms = gen_covar_syms(rn_mode_1)

#collect(Iterators.flatten(symb_covar))
# Mode switching jump rates.
l_jump_rates = [(1-d_1)*rs*x_4*scaling]
# Full symbolic equation for iterative.
l_sym_iter = vcat(sum(map((rn,w) -> w*rn.f_symfuncs, l_rns, [w_1, w_2])),
                  sum(map((rn,w) -> w*gen_noise_covar(rn, l_cov_syms), l_rns, [w_1, w_2])),
                  l_jump_rates)
l_lhs_iter = unique(vcat(vec(rn_mode_1.syms),vec(l_cov_syms), [:d_1]))
l_odes_iter = sym_to_expr(l_sym_iter, l_lhs_iter, l_rns[1].params, [:w_1, :w_2])
# Full symbolic equations for coupling
l_filter_eqns = [-pi_1*rs*x_4*scaling, pi_1*rs*x_4*scaling]
l_sym_cpl = vcat((rn_mode_1.f_symfuncs*pi_1 + rn_mode_2.f_symfuncs*pi_2),
                 (gen_noise_covar(rn_mode_1, l_cov_syms)*pi_1 + 
                  gen_noise_covar(rn_mode_2, l_cov_syms)*pi_2),
                  l_filter_eqns)                  
l_lhs_cpl = unique(vcat(l_rns[1].syms, vec(l_cov_syms), [:pi_1, :pi_2]))
l_odes_cpl = sym_to_expr(l_sym_cpl, l_lhs_cpl, l_rns[1].params, [])


# Make the ODE function...
@generated function l_iter(du, u, p, t)
    quote
        du .= $(l_odes_iter)
    end
end

@generated function l_cpl(du, u, p, t)
    quote
        du .= $(l_odes_cpl)
    end
end

function w0(t)
    return [1.0, 0.0]
end

function run_noise(iter, u0_cpl, u0_iter, params)
    println("Linear noise  -------------------------")
    println("Solving the direct coupling problem")
    time_cpl = @elapsed sol_cpl = run_single(l_cpl, u0_cpl, tspan, params, [], tplot)
    println("Time taken $time_cpl s")

    println("Solving the weighted ODE problem")
    time_iter = @elapsed sol_iter = run_iterative(l_iter, u0_iter, tspan, params, 
                                                  length(u0_iter), length(l_rns), w0, tplot)
    return sol_cpl, sol_iter
end

function run_noise_cpl(iter, u0_cpl, params)
    time_cpl = @elapsed sol_cpl = run_single(l_cpl, u0_cpl, tspan, params, [], tplot)
    return sol_cpl
end

function run_noise_iter(iter, u0_iter, params)
    time_iter = @elapsed sol_iter = run_iterative(l_iter, u0_iter, tspan, params, 
                                                  length(u0_iter), length(l_rns), w0, tplot)
    return sol_iter[2][2]
end

