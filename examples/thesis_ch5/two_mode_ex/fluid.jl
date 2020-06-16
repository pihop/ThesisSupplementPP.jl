using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using JLD2
using LinearAlgebra
using SymEngine

include("rn.jl")
include("../common2.jl")

rs = symbols(:rs)
d_1 = symbols(:d_1)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
x_4 = symbols(:x_4)

rns = [rn_mode_1, rn_mode_2]
# Full symbolic equation
symbolic_iter = sum(map((rn,w) -> w*rn.f_symfuncs, rns, [w_1, w_2]))
# Mode switching dist.
jump_rates = [(1-d_1)*rs*x_4*scaling]
symbolic_iter = vcat(symbolic_iter, jump_rates)
odes_iter = sym_to_expr(symbolic_iter, vcat(rns[1].syms, :d_1), rns[1].params, [:w_1, :w_2])

symbolic_coupl = rn_mode_1.f_symfuncs*(1-d_1) + rn_mode_2.f_symfuncs*d_1
symbolic_coupl= vcat(symbolic_coupl, jump_rates)
odes_coupl = sym_to_expr(symbolic_coupl, vcat(rns[1].syms, :d_1), rns[1].params, [])

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
    return [1.0, 0.0]
end
   
println("Solving the direct coupling problem")
@time sol_cpl = run_single(f_coupl, u0, tspan, p, [], tplot)

println("Solving the weighted ODE problem")
@time hitting_t, sol_iter = run_iterative(f_iter, u0, tspan, p, 4, 1, w0, tplot)

jldopen(joinpath(@__DIR__, "data/fluid_res.jld2"), "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "initial", u0)
    write(file, "params", p)
    write(file, "sol_cpl", sol_cpl)
    write(file, "hitting_t", hitting_t)
    write(file, "sol_iter", sol_iter)
end

