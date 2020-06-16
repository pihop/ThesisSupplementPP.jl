using ThesisSupplementPP
using JLD2
using SymEngine

include("../common2.jl")
include("rn.jl")
include("parametrisation.jl")

rs = symbols(:rs)
d_1 = symbols(:d_1)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
x_4 = symbols(:x_4)
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)

rns = [rn_mode_1, rn_mode_2]

# Mode switching dist.
jump_rates = [(1-d_1)*rs*x_4*scaling]

# Full symbolic equation
symbolic_iter = sum(map((rn,w) -> w*rn.f_symfuncs, rns, [w_1, w_2]))
covar_syms = gen_covar_syms(rn_mode_1)
symb_covar = sum(map((rn,w) -> w*gen_noise_covar(rn, covar_syms), rns, [w_1, w_2]))
symbolic_iter = vcat(symbolic_iter, collect(Iterators.flatten(symb_covar)), jump_rates)
lhs_iter = unique(vcat(rns[1].syms, vec(covar_syms), :d_1))
odes_iter = sym_to_expr(symbolic_iter, lhs_iter, rn_full.params, [:w_1, :w_2])


filtering = [-pi_1*rs*x_4*scaling, pi_1*rs*x_4*scaling]

symbolic_coupl = rn_mode_1.f_symfuncs*pi_1 + rn_mode_2.f_symfuncs*pi_2
symbolic_coupl_cov = (gen_noise_covar(rn_mode_1, covar_syms)*pi_1 + 
                      gen_noise_covar(rn_mode_2, covar_syms)*pi_2)
symbolic_coupl = vcat(symbolic_coupl, collect(Iterators.flatten(symbolic_coupl_cov)), filtering)
lhs_cpl = unique(vcat(rns[1].syms, vec(covar_syms), [:pi_1, :pi_2]))
odes_coupl = sym_to_expr(symbolic_coupl, lhs_cpl, rns[1].params, [])

# Make the ODE function...
@generated function f_lne_iter(du, u, p, t)
    quote
        du .= $(odes_iter)
    end
end

@generated function f_lne_coupl(du, u, p, t)
    quote
        du .= $(odes_coupl)
    end
end

function w0(t)
    return [1.0, 0.0]
end

#Initial conditions
u0 = u0_inhom./ scaling 
u0_lne_iter = vcat(u0, fill(0.0, length(lhs_iter) - length(u0)))
u0_lne_cpl = vcat(u0, fill(0.0, length(lhs_cpl) - length(u0)))
u0_lne_cpl[end-1] = 1.0 
  
println("Solving the direct coupling problem")
@time sol_cpl = run_single(f_lne_coupl, u0_lne_cpl, tspan, p, [], tplot)

println("Solving the weighted ODE problem")
@time hitting_t, sol_iter = run_iterative(f_lne_iter, u0_lne_iter, tspan, p, length(lhs_iter), length(filtering), w0, tplot)

jldopen("data/lna_res.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "params", p)
    write(file, "sol_cpl", sol_cpl)
    write(file, "hitting_t", hitting_t)
    write(file, "sol_iter", sol_iter[end])
end

