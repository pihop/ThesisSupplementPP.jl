using ThesisSupplementPP
using JLD2
using SymEngine

include("../common2.jl")
include("rn.jl")

rs = symbols(:rs)
d_1 = symbols(:d_1)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
x_4 = symbols(:x_4)
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)

scaling = 100
tspan = (0., 10.)
p = (1.0, 0.2)
tplot = 0.0:0.01:10.0
u0 = [0., 1., 0., 0.]

rns = [rn_mode_1, rn_mode_2]
mc_modes = [make_symbolic_mc(rn,2) for rn in rns]

# Mode switching dist.
jump_rates = [(1-d_1)*rs*x_4]

symbolic_iter = sum(map((mc,w) -> w*mc[2], mc_modes, [w_1, w_2]))
syms = mc_modes[1][1]
symbolic_iter = vcat(symbolic_iter, jump_rates)
lhs_iter = vcat(Symbol.(vec(syms)), :d_1)
odes_iter = sym_to_expr(symbolic_iter, lhs_iter, rn_full.params, [:w_1, :w_2])

filtering = [-pi_1*rs*x_4, pi_1*rs*x_4]
lhs_cpl = vcat(Symbol.(vec(syms)), [:pi_1, :pi_2])
symbolic_coupl = mc_modes[1][2]*pi_1 + mc_modes[2][2]*pi_2
symbolic_coupl = vcat(symbolic_coupl, filtering)
odes_coupl = sym_to_expr(symbolic_coupl, lhs_cpl, rn_full.params, [])

jldopen("mc_eqns.jld2", "w") do file
    write(file, "eqns_coupl", odes_coupl)
    write(file, "eqns_coupl_symb", convert.(Expr, symbolic_coupl))
    write(file, "eqns_iter", odes_iter)
    write(file, "eqns_iter_symb", convert.(Expr, symbolic_iter))
end

# Make the ODE function...
@generated function f_mc_iter(du, u, p, t)
    quote
        du .= $(odes_iter)
    end
end

@generated function f_mc_coupl(du, u, p, t)
    quote
        du .= $(odes_coupl)
    end
end

function w0(t)
    return [1.0, 0.0]
end

#Initial conditions
u0_mc_iter = vcat(u0, fill(0.0, length(lhs_iter) - length(u0)))
u0_mc_cpl = vcat(u0, fill(0.0, length(lhs_cpl) - length(u0)))
u0_mc_cpl[end-1] = 1.0 
u0_mc_iter[2] = 100.0
u0_mc_cpl[2] = 100.0

  
println("Solving the direct coupling problem")
@time sol_cpl = run_single(f_mc_coupl, u0_mc_cpl, tspan, p, [], tplot)

println("Solving the weighted ODE problem")
@time hitting_t, sol_iter = run_iterative(f_mc_iter, u0_mc_iter, tspan, p, length(lhs_iter), length(filtering), w0, tplot)

jldopen("data/mc_res.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "params", p)
    write(file, "sol_cpl", sol_cpl)
    write(file, "hitting_t", hitting_t)
    write(file, "sol_iter", sol_iter)
end

