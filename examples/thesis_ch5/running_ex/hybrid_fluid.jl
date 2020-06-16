using ThesisSupplementPP
using JLD2
using SymEngine

include("rn.jl")
include("../common2.jl")
include("parametrisation.jl")

rs = symbols(:rs)
d_1 = symbols(:d_1)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
x_4 = symbols(:x_4)
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)

rns = [rn_mode_1, rn_mode_2]
# Full symbolic equation
symbolic_iter = sum(map((rn,w) -> w*rn.f_symfuncs, rns, [w_1, w_2]))
# Mode switching dist.
jump_rates = [(1-d_1)*rs*x_4*scaling]
symbolic_iter = vcat(symbolic_iter, jump_rates)
odes_iter = sym_to_expr(symbolic_iter, vcat(rns[1].syms, :d_1), rns[1].params, [:w_1, :w_2])

filtering = [-pi_1*rs*x_4*scaling, pi_1*rs*x_4*scaling]
symbolic_coupl = rn_mode_1.f_symfuncs*pi_1 + rn_mode_2.f_symfuncs*pi_2
symbolic_coupl= vcat(symbolic_coupl, filtering)
odes_coupl = sym_to_expr(symbolic_coupl, vcat(rns[1].syms, [:pi_1, :pi_2]), rns[1].params, [])

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

# Initial conditions matching the example description. The system starts in mode
# 1.
u0_cpl = vcat(u0_inhom ./ scaling, [1.0, 0.])
println("Solving the direct coupling problem")
@time sol_cpl = run_single(f_coupl, u0_cpl, tspan, p, [], tplot)

# Initial conditions. Probability of jump having happend at time 0 is 0.
u0_iter = vcat(u0_inhom ./ scaling, [0.])
println("Solving the weighted ODE problem")
@time hitting_t, sol_iter = run_iterative(f_iter, u0_iter, tspan, p, 4, 1, w0, tplot)

jldopen("data/fluid_res.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "params", p)
    write(file, "sol_cpl", sol_cpl)
    write(file, "sol_iter", sol_iter)
    write(file, "hitting_t", hitting_t)
end

