# Construct and simulate the hybrid linear noise process. Running example of
# Chapter 4.

using DiffEqBiological
using ThesisSupplementPP
using DiffEqBase, StochasticDiffEq
using JLD2
using LinearAlgebra
using LaTeXStrings
using SymEngine

# Reaction networks.
include("rn.jl")

# Parse the arguments for scaling and number of trajectories.
scaling = parse(Int,ARGS[1])
trajectories = parse(Int, ARGS[2]) 

# Time-span of simulation
tspan = (0., 10.)
# Parameters for the model.
p = (1.0, 0.2)

# Setting up the symbols.
rs = symbols(:rs)
d_1 = symbols(:d_1)
w_1 = symbols(:w_1)
w_2 = symbols(:w_2)
x_4 = symbols(:x_4)
pi_1 = symbols(:pi_1)
pi_2 = symbols(:pi_2)

# Two dynamic modes.
rns = [rn_mode_1, rn_mode_2]
# Symbols for elements in the covariance matrix.
covar_syms = gen_covar_syms(rn_mode_1)
# Fluid odes.
odes = [rn.f_symfuncs for rn in rns]
# Diffusion odes.
odes_diff = [gen_diffusion(rn) for rn in rns]
# LHS. 
lhs = rns[1].syms

# Create the expressions.
odes_expr = [sym_to_expr(ode, lhs, rn_full.params, []) for ode in odes]
odes_diff_expr = [sym_to_expr(ode, lhs, rn_full.params, []) for ode in odes_diff]

# Make the ODE function...
@generated function f_lna_hyb(du, u, p, t)
    quote
        if p[2] == 1
            du .= $(odes_expr[1])
        else
            du .= $(odes_expr[2])
        end
    end
end

# Noise function from the expressions. Using the parameter p[2] for mode
# switching. p[1] is reserved for a tuple of model parameters.
@generated function g_lna_hyb(du, u, p, t)
    quote
        if p[2] == 1
            du .= $(odes_diff_expr[1])
        else
            du .= $(odes_diff_expr[2])
        end
        return du
    end
end

# Noise function rescaled and correctly shaped.  
function gn(du, u, p, t)
    du .= scaling^(-0.5)*reshape(g_lna_hyb(vcat(du...), vcat(u...), p, t), size(covar_syms))
end

# Rate of mode switching.
rate_switch(u,p,t) = p[1][2]*u[4]*scaling

function affect_switch!(integrator)
    integrator.p = (p, 0)
end

jump_switch = VariableRateJump(rate_switch,affect_switch!)
u0_lna = zeros(length(lhs))
u0_lna[2] = 1.0

sde_p = SDEProblem(f_lna_hyb, gn, u0_lna, tspan, (p, 1), noise_rate_prototype=zeros(size(covar_syms)))
jump_p = JumpProblem(sde_p, Direct(), jump_switch)

println("The monte carlo simulation")
sim_time = @elapsed sol = solve(EnsembleProblem(jump_p), EulerHeun(), dt = 0.01; dense=false, trajectories=trajectories)
println("Time taken $sim_time")

jldopen("data/hybrid_lna_$scaling.jld2", "w") do file
    write(file, "time", tspan)
    write(file, "pop", scaling)
    write(file, "initial", u0_lna)
    write(file, "params", p)
    write(file, "sim/time", sim_time)
    write(file, "sim/sol", sol)
    write(file, "sim/sum", sum)
end
