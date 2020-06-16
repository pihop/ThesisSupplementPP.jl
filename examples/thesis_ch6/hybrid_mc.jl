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

m_rns = [rn_mode_1, rn_mode_2]
m_eqns = [make_symbolic_mc(rn,2) for rn in m_rns]
m_syms = m_eqns[1][1]

# Mode switching jump rates.
jump_rates = [(1-d_1)*rs*x_4]
# Full symbolic equation for iterative.
m_sym_iter = vcat(sum(map((rn,w) -> w*rn[2], m_eqns, [w_1, w_2])), jump_rates)
m_lhs_iter = vcat(Symbol.(vec(m_syms)), [:d_1])
m_odes_iter = sym_to_expr(m_sym_iter, m_lhs_iter, m_rns[1].params, [:w_1, :w_2])
# Full symbolic equations for coupling
m_sym_cpl = (m_eqns[1][2]*pi_1 + m_eqns[2][2]*pi_2)
filter_eqns = [-pi_1*rs*x_4, pi_1*rs*x_4]
m_sym_cpl= vcat(m_sym_cpl, filter_eqns)
m_lhs_cpl = vcat(Symbol.(vec(m_syms)), [:pi_1, :pi_2])
m_odes_cpl = sym_to_expr(m_sym_cpl, m_lhs_cpl, m_rns[1].params, [])

# Make the ODE function...
@generated function m_iter(du, u, p, t)
    quote
        du .= $(m_odes_iter)
    end
end

@generated function m_cpl(du, u, p, t)
    quote
        du .= $(m_odes_cpl)
    end
end

function w0(t)
    return [1.0, 0.0]
end

function run_mc(iter, u0_cpl, u0_iter, params)
    println("Moment closure -------------------------")
    println("Solving the direct coupling problem")
    time_cpl = @elapsed sol_cpl = run_single(m_cpl, u0_cpl, tspan, params, [], tplot)
    println("Time taken $time_cpl s")

    println("Solving the weighted ODE problem")
    time_iter = @elapsed sol_iter = run_iterative(m_iter, u0_iter, tspan, params, 
                                                  length(u0_iter), length(m_rns), w0, tplot)
    return sol_cpl, sol_iter
end

function run_mc_cpl(iter, u0_cpl, params)
    time_cpl = @elapsed sol_cpl = run_single(m_cpl, u0_cpl, tspan, params, [], tplot)
    return sol_cpl
end

function run_mc_iter(iter, u0_iter, params)
    time_iter = @elapsed sol_iter = run_iterative(m_iter, u0_iter, tspan, params, 
                                                  length(u0_iter), length(m_rns), w0, tplot)
    return sol_iter[2][2]
end

