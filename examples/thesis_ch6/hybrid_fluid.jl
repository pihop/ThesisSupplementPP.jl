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

f_rns = [rn_mode_1, rn_mode_2]

# Mode switching jump rates.
jump_rates = [(1-d_1)*rs*x_4*scaling]
# Full symbolic equation for iterative.
f_sym_iter = vcat(sum(map((rn,w) -> w*rn.f_symfuncs, f_rns, [w_1, w_2])), jump_rates)
f_lhs_iter = vcat(f_rns[1].syms, [:d_1])
f_odes_iter = sym_to_expr(f_sym_iter, f_lhs_iter, f_rns[1].params, [:w_1, :w_2])
# Full symbolic equations for coupling
f_sym_cpl = (rn_mode_1.f_symfuncs*pi_1 + rn_mode_2.f_symfuncs*pi_2)
f_filter_eqns= [-pi_1*rs*x_4*scaling, pi_1*rs*x_4*scaling]
f_sym_cpl= vcat(f_sym_cpl, f_filter_eqns)
f_lhs_cpl = vcat(f_rns[1].syms, [:pi_1, :pi_2])
f_odes_cpl = sym_to_expr(f_sym_cpl, f_lhs_cpl, f_rns[1].params, [])

# Make the ODE function...
@generated function f_iter(du, u, p, t)
    quote
        du .= $(f_odes_iter)
    end
end

@generated function f_cpl(du, u, p, t)
    quote
        du .= $(f_odes_cpl)
    end
end

function w0(t)
    return [1.0, 0.0]
end

function run_fluid(iter, u0_cpl, u0_iter, params)
    println("Fluid -------------------------------")
    println("Solving the direct coupling problem")
    time_cpl = @elapsed sol_cpl = run_single(f_cpl, u0_cpl, tspan, params, [], tplot)
    println("Time taken $time_cpl s")

    println("Solving the weighted ODE problem")
    time_iter = @elapsed sol_iter = run_iterative(f_iter, u0_iter, tspan, params, 
                                                  length(u0_iter), length(f_rns), w0, tplot)
    return sol_cpl, sol_iter
end

function run_fluid_cpl(iter, u0_cpl,  params)
    time_cpl = @elapsed sol_cpl = run_single(f_cpl, u0_cpl, tspan, params, [], tplot)
    return sol_cpl
end

function run_fluid_iter(iter, u0_iter, params)
    time_iter = @elapsed sol_iter = run_iterative(f_iter, u0_iter, tspan, params, 
                                                  length(u0_iter), length(f_rns), w0, tplot)
    return sol_iter[2][2]
end
