using DiffEqBase, DiffEqJump, OrdinaryDiffEq, StochasticDiffEq
using Interpolations
using Plots
using JLD2
using FileIO
using LinearAlgebra
using LaTeXStrings
using Plots.PlotMeasures
using SymEngine
pyplot()

PyPlot.rc("text", usetex="true")
fonts = font(12)

# Simulation output files.
stochastic_sim_file = "data/ssa.jld2"
inhom_stochastic_sim_file = "data/ssa_inhom.jld2"
hybrid_fluid_file = "data/hybrid_fluid.jld2" 
hybrid_lna_file = "data/hybrid_lna.jld2"

# Simulation summaries
sim_sum = load(stochastic_sim_file, "sim_sum")
sim_sum_inhom = load(inhom_stochastic_sim_file, "sim_sum")

# Linear noise approximation
lna_sol_cpl = load("data/lna_res.jld2", "sol_cpl")
lna_sol_iter = load("data/lna_res.jld2", "sol_iter")
lna_hitting_t = load("data/lna_res.jld2", "hitting_t")

# Moment closure-based approximation. This is equivalent to LNA results due to 
# only linear reactions in the reaction network of the model.
#mc_sol_cpl = load("data/mc_res.jld2", "sol_cpl")
#mc_sol_iter = load("data/mc_res.jld2", "sol_iter")
#mc_hitting_t = load("data/mc_res.jld2", "hitting_t")

label_list = ["(1,2)", "(1,1)", "(2,1)", "(2,2)"]

# First set of plots. Time-inhomogeneous simulations. Compared with the full
# stochastic simulation.
idx_toplot = [2,4]
plt_sols = plot(ylabel=latexstring("\\mbox{Population at a given location}"), 
                xlabel=latexstring("\\mbox{Time}"), 
                xlim=(0.0, tplot[end]),
                ylim=(0.0, sum(sim_sum[1])),
                legend=:best, 
                fg_legend=:transparent,
                bg_legend=:transparent,
                size=(700,300),
                xtickfont=fonts, 
                ytickfont=fonts,
                guidefont=fonts, 
                legendfont=fonts)

for (i, idx_sol) in enumerate(idx_toplot)
    # Plotting one standard deviation around the mean.
    mean_sim = getindex.(sim_sum.u.u, idx_sol)
    upper = mean_sim .+ sqrt.(getindex.(sim_sum.v.u, idx_sol))
    lower = max.(0.0, mean_sim .- sqrt.(getindex.(sim_sum.v.u, idx_sol)))

    mean_simi = getindex.(sim_sum_inhom.u.u, idx_sol)
    upperi = mean_sim .+ sqrt.(getindex.(sim_sum_inhom.v.u, idx_sol))
    loweri = max.(0.0, mean_sim .- sqrt.(getindex.(sim_sum_inhom.v.u, idx_sol)))

    plt_sols = scatter!([-1], [-1], color=i, label=latexstring("\\mbox{Location }", string(label_list[idx_sol])))

    plt_sols = plot!(tplot, upper; color=i, lw=1.1, line=(:solid), labels="")
    plt_sols = plot!(tplot, lower; color=i, lw=1.1, line=(:solid), labels="")
#
    plt_sols = plot!(tplot, upperi; color=i, lw=1.1, line=(:dash), labels="")
    plt_sols = plot!(tplot, loweri; color=i, lw=1.1, line=(:dash), labels="")
end
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:solid, 1.5), linealpha=0.5, label="Original mode-switching process")
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:dash, 1.5), linealpha=0.5, label="Time-inhomogeneous approximation ")
savefig(plt_sols, "plots/ssa_inhom.pdf")

# Visual comparison of means from SSA and the iterative and direct coupling for
# the fluid construction.
plt_sols_1 = plot(ylabel="Population density at a given location", 
                xlabel="Time", 
                legend=:right, 
                fg_legend=:transparent,
                bg_legend=:transparent,
                size=(700,300),
                xlims=(0.0, 8.0),
                xtickfont=fonts, 
                ytickfont=fonts,
                guidefont=fonts, 
                legendfont=fonts)

for (i, idx) in enumerate(idx_toplot)
    mean_sim = getindex.(sim_sum.u.u, idx)
    # Fluid results scaled up to population size defined by the scaling
    # variable.
    mean_cpl = scaling*getindex.(lna_sol_cpl(tplot).u, idx)
    # The last solution set in the sol_iter corresponds to the final iteration
    # in the solver.
    mean_iter = scaling*getindex.(lna_sol_iter(tplot).u, idx)
    plt_sols_1 = plot!(tplot, mean_sim; color=i, lw=1.1, line=(0.5, :solid), 
                     labels=string("Location ", label_list[idx], "(SSA)"))
    plt_sols_1 = plot!(tplot, mean_cpl; color=i, lw=1.1, line=(0.9, :dash), 
                     labels=string("Location ", label_list[idx], "(Direct coupling)"))
    plt_sols_1 = plot!(tplot, mean_iter; color=i, lw=1.1, line=(:dashdot), 
                     labels=string("Location ", label_list[idx], "(Iterative)"))
end

savefig(plt_sols_1, "plots/fluid_mean.pdf")

# Hitting time plots.
plt_weight = plot(ylabel=L"Probability $Z(t) = 1$", xlabel="Time", 
                  legend=:right, 
                  fg_legend=:transparent,
                  bg_legend=:transparent,
                  xlims=(0, 5.0), 
                  size=(700,300),
                  xtickfont=fonts, 
                  ytickfont=fonts,
                  guidefont=fonts, 
                  legendfont=fonts)

# Empirical mean values of the mode switching variable. The variable idx 5 in
# the sim_sum tracks whether mode switch has happened or not.
empirical_hitting = LinearInterpolation(sim_sum.t, getindex.(sim_sum.u.u, 5))
plt_weight = plot!(tplot, t->empirical_hitting(t), lw=1.5, line=(0.7, :solid), label="Stochastic simulation")
plt_weight = plot!(tplot, t->lna_sol_cpl(t)[end], lw=1.5, line=(:dash), label="Direct coupling based approximation")
#plt_weight = plot!(tplot, t->lna_sol_iter[1](t)[end], lw=1.5, line=(:dashdot), label="Iterative approximation")
plt_weight = plot!(tplot, t -> lna_hitting_t[1](t), lw=1.5, line=(:dashdot), label="Iterative approximation")
savefig(plt_weight, "plots/hitting_fluid.pdf")

# Relative error for fluid approximation constructions.
# At index rel_idx of the solution
rel_idx = 4

plt_rel = plot(ylabel=L"Relative error for location $(1,1)$ ", 
           xlabel="Time", 
           legend=:right, 
           fg_legend=:transparent,
           bg_legend=:transparent,
           size=(700,300),
           xtickfont=fonts, 
           ytickfont=fonts,
           guidefont=fonts, 
           legendfont=fonts,
           ylims=(0.0, 0.3)
          )

function rel_err(x, y)
    return abs.(x .- y) ./ y
end

# Relative errors.
error_iter = rel_err(getindex.(sim_sum.u.u, rel_idx), 
                     scaling*getindex.(lna_sol_iter(tplot).u, rel_idx))
error_cpl = rel_err(getindex.(sim_sum.u.u, rel_idx), 
                    scaling*getindex.(lna_sol_cpl(tplot).u, rel_idx))
plt_rel = plot!(tplot, error_iter; fillalpha=0.2, lw=1.5, line=(:dashdot), labels="Iterative")
plt_rel = plot!(tplot, error_cpl; fillalpha=0.2, lw=1.5, line=(:dash), labels="Direct coupling")
savefig(plt_rel, "plots/relative_fluid.pdf")


# Linear noise based constructions.
plt_sols = plot(ylabel="Population density at a given location ", 
                xlabel="Time", 
                legend=:right, 
                fg_legend=:transparent,
                bg_legend=:transparent,
                size=(700,300),
                xtickfont=fonts, 
                ytickfont=fonts,
                guidefont=fonts, 
                legendfont=fonts)

# The first element of the tuples is the idx of the mean and the second the
# corresponding variance.
to_plot_pairs = [(2, 9), (4, 14)]

for (i,pair) in enumerate(to_plot_pairs)
    # Appropriate scaling of the variances. Iterative
    iter_var(t) = (scaling^(0.5))sqrt(lna_sol_iter(t)[pair[2]])
    iter_mean(t) = scaling*lna_sol_iter(t)[pair[1]]
    iter_upper(t) = iter_var(t) + iter_mean(t)
    iter_lower(t) = -iter_var(t) + iter_mean(t)

    plt_sols = plot!(tplot, iter_upper; color=i, lw=1.1, line=(:solid), 
                     labels=string("Location ", label_list[pair[1]], "(Iterative)"))
    plt_sols = plot!(tplot, iter_lower; color=i, lw=1.1, line=(:solid), labels="")

    # Appropriate scaling of the variances. Direct coupling
    cpl_var(t) = (scaling^(0.5))*sqrt(lna_sol_cpl(t)[pair[2]])
    cpl_mean(t) = scaling*lna_sol_cpl(t)[pair[1]]
    cpl_upper(t) = cpl_var(t) + cpl_mean(t)
    cpl_lower(t) = -cpl_var(t) + cpl_mean(t)
    plt_sols = plot!(tplot, cpl_upper; color=i, lw=1.1, line=(:dash), 
                     labels=string("Location ", label_list[pair[1]], "(Coupled)"))
    plt_sols = plot!(tplot, cpl_lower; color=i, lw=1.1, line=(:dash), labels="")
end

savefig(plt_sols, "plots/lna.pdf")

#Standard deviations from LNA based approximations compared to SSA.
to_plot = (4, 14)
plt_var = plot(ylabel=latexstring("\\mbox{Standard deviation for location \$(1,1)\$}"), 
               xlabel="Time", 
               xlim=tspan,
               ylim=(0.0, 20.0),
               legend=:best, 
               fg_legend=:transparent,
               bg_legend=:transparent,
               size=(700,300),
               xtickfont=fonts, 
               ytickfont=fonts,
               guidefont=fonts, 
               legendfont=fonts)

var_iter(t) = scaling^(0.5)*real(sqrt.(abs.(lna_sol_iter(t)[to_plot[2]])))    
var_cpl(t) = scaling^(0.5)*real(sqrt.(abs.(lna_sol_cpl(t)[to_plot[2]])))
var_sim  = sqrt.(getindex.(sim_sum.v.u, to_plot[1]))
plt_rel = plot!(tplot, var_sim; fillalpha=0.2, lw=1.1, line=(:solid), labels=latexstring("\\mbox{Stochastic simulation}"))
plt_var = plot!(tplot, var_iter; fillalpha=0.2, lw=1.1, line=(:dashdot), labels=latexstring("\\mbox{Iterative}"))
plt_var = plot!(tplot, var_cpl; fillalpha=0.2, lw=1.1, line=(0.5, :dash), labels=latexstring("\\mbox{Direct coupling}"))
savefig(plt_var, "plots/sd.pdf")
