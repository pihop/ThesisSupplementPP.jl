using ThesisSupplementPP
using StochasticDiffEq
using DiffEqBase
using DiffEqJump
using JLD2
using FileIO
using MAT
using LaTeXStrings
using StatsBase
using Plots
using Interpolations
pyplot()
PyPlot.rc("text", usetex="true")

scaling = 100
# Simulation output files.
stochastic_sim_file = "data/ssa_$scaling.jld2"
hybrid_fluid_file = "data/hybrid_fluid_$scaling.jld2" 
hybrid_lna_file = "data/hybrid_lna_$scaling.jld2"

# Simulation summaries.
sim_sum = load(stochastic_sim_file, "sim_sum")
sim_sum_hf = load(hybrid_fluid_file, "sim/sum")
# Simulation summary broken for SDEs due to intepolation bug in
# StochasticDiffEq.jl. Loading all trajectories for a workaround.
sim_sol_hln = load(hybrid_lna_file, "sim/sol")

# Load in the simulation output file from the CERENA model.
mat_file = matopen("matlab/running_MCM.mat")
mat_sol = read(mat_file, "System_running")["sol"]
close(mat_file)

# Time interval to plot.
tplot = 0.0:0.01:10.0
# Labels for each of the solution indices. Eg. the first component in the
# solution time-series correspondsto location of the model (1,2).
label_list = ["(1,2)", "(1,1)", "(2,1)", "(2,2)"]
# Variable indices to plot.
idxs_to_plot = [2,4]

plt_sols = plot(ylabel=latexstring("\\mbox{Population at a given location}"), 
                xlabel=latexstring("\\mbox{Time}"), 
                xlim=(0.0, tplot[end]),
                ylim=(0.0, sum(sim_sum[1])),
                legend=:best, 
                fg_legend=:transparent,
                bg_legend=:transparent,
                size=(400,300),
                xtickfont=font(10), 
                ytickfont=font(10),
                guidefont=font(10), 
                legendfont=font(10))

# Hybrid fluid plots.
for (i, idx_sol) in enumerate(idxs_to_plot)
    # Get the full simulation mean and quantiles.
    mean_sim = getindex.(sim_sum.u.u, idx_sol)
    upper = mean_sim .+ sqrt.(getindex.(sim_sum.v.u, idx_sol))
    lower = max.(0.0, mean_sim .- sqrt.(getindex.(sim_sum.v.u, idx_sol)))

    # Get the hybrid fluid simulation mean and quantiles.
    # The hybrid fluid trajectories are scaled from [0,1] to [0, scaling] in
    # order to compare them with stochastic simulation.
    mean_sim_h = scaling .* getindex.(sim_sum_hf.u.u, idx_sol)
    upper_h = mean_sim_h .+ scaling*sqrt.(getindex.(sim_sum_hf.v.u, idx_sol))
    lower_h = max.(0.0, mean_sim .- scaling*sqrt.(getindex.(sim_sum_hf.v.u, idx_sol)))
    
    # The following point is outside the plotting region and is used to display
    # the legend.
    plt_sols = scatter!([-1], [-1], color=i, label=latexstring("\\mbox{Location }", string(label_list[idx_sol])))

    # Plot the trajectories.
    plt_sols = plot!(tplot, upper; color=i, lw=1.1, line=(:solid), labels="")
    plt_sols = plot!(tplot, lower; color=i, lw=1.1, line=(:solid), labels="")

    plt_sols = plot!(tplot, upper_h; color=i, lw=1.1, line=(:dash), labels="")
    plt_sols = plot!(tplot, lower_h; color=i, lw=1.1, line=(:dash), labels="")

end

# These points are used for the legend.
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:solid, 1.5), linealpha=0.5, label="Stochastic simulation")
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:dash, 1.5), linealpha=0.5, label="Hybrid fluid approximation")

savefig(plt_sols, "plots/hybrid_fluid_sim.pdf")

# Clear the plot.
plt_sols = plot(ylabel=latexstring("\\mbox{Population at a given location}"), 
                xlabel=latexstring("\\mbox{Time}"), 
                xlim=(0.0, tplot[end]),
                ylim=(0.0, sum(sim_sum[1])),
                legend=:best, 
                fg_legend=:transparent,
                bg_legend=:transparent,
                size=(400,300),
                xtickfont=font(10), 
                ytickfont=font(10),
                guidefont=font(10), 
                legendfont=font(10))

# Hybrid linear noise plot
for (i, idx_sol) in enumerate(idxs_to_plot)
    # Interpolations for SDE solutions structures do not work correctly in the
    # StochasticDiffEq.jl pacakge. The following is a workaround. 
    interpolations = [LinearInterpolation(sim_sol_hln[k].t, getindex.(sim_sol_hln[k].u, idx_sol[1])) for k=1:5000]
    interp(t) = [it(t) for it in interpolations]
    mean_sim_h = [scaling*mean(interp(t)) for t in tplot]
    sd_h = [scaling* sqrt(totalvar(interp(t))) for t in tplot]
    upper_h = mean_sim_h .+ sd_h 
    lower_h = max.(0.0, mean_sim_h .- sd_h)

    # Stochastic simulation of the full process.
    mean_sim = getindex.(sim_sum.u.u, idx_sol)
    upper = mean_sim .+ sqrt.(getindex.(sim_sum.v.u, idx_sol))
    lower = max.(0.0, mean_sim .- sqrt.(getindex.(sim_sum.v.u, idx_sol)))

    plt_sols = scatter!([-1], [-1], color=i, label=latexstring("\\mbox{Location }", string(label_list[idx_sol])))

    plt_sols = plot!(tplot, upper; color=i, lw=1.1, line=(:solid), labels="")
    plt_sols = plot!(tplot, lower; color=i, lw=1.1, line=(:solid), labels="")

    plt_sols = plot!(tplot, upper_h; color=i, lw=1.1, line=(:dash), labels="")
    plt_sols = plot!(tplot, lower_h; color=i, lw=1.1, line=(:dash), labels="")
end

plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:solid, 1.5), linealpha=0.5, label="Stochastic simulation")
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:dash, 1.5), linealpha=0.5, label="Hybrid noise approximation")

savefig(plt_sols, "plots/hybrid_lna_sim.pdf")

# Clear the plot.
plt_sols = plot(ylabel=latexstring("\\mbox{Population at a given location}"), 
                xlabel=latexstring("\\mbox{Time}"), 
                xlim=(0.0, tplot[end]),
                ylim=(0.0, sum(sim_sum[1])),
                legend=:best, 
                fg_legend=:transparent,
                bg_legend=:transparent,
                size=(400,300),
                xtickfont=font(10), 
                ytickfont=font(10),
                guidefont=font(10), 
                legendfont=font(10))

# MCM plots. 
for (i, sim_idx) in enumerate([(2,9), (4,14)])
    # The correspondence of means and variances is completely undocumented in
    # CERENA but can be found by looking at the symbolic expressions it
    # constructs. So the column 2 in the matlab output array turns out to be a
    # mean and column 9 the corresponding variance. Similar with 4 and 14.

    # SSA simulation results.
    mean_sim = getindex.(sim_sum.u.u, sim_idx[1])
    upper = mean_sim .+ sqrt.(getindex.(sim_sum.v.u, sim_idx[1]))
    lower = max.(0.0, mean_sim .- sqrt.(getindex.(sim_sum.v.u, sim_idx[1])))

    # MCM results.
    mean_mcm = mat_sol["y"][:,sim_idx[1]]
    upper_mcm = mean_mcm .+ sqrt.(((mat_sol["y"][:,sim_idx[2]])))
    lower_mcm = mean_mcm .- sqrt.(((mat_sol["y"][:,sim_idx[2]])))

    # Appear off the plot. Used for the legend.
    plt_sols = scatter!([-1], [-1], color=i, label=latexstring("\\mbox{Location }",
                                                               string(label_list[sim_idx[1]])))


    # Plot the results.
    plt_sols = plot!(tplot, upper; color=i, lw=1.1, line=(:solid), labels="")
    plt_sols = plot!(tplot, lower; color=i, lw=1.1, line=(:solid), labels="")
    plt_sols = plot!(tplot[1:1000], upper_mcm; color=i, lw=1.1, line=(:dash), labels="")
    plt_sols = plot!(tplot[1:1000], lower_mcm; color=i, lw=1.1, line=(:dash), labels="")
end
# Appear off the plot. Used for the legend.
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:solid, 1.5), linealpha=0.5, label="Stochastic simulation")
plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:dash, 1.5), linealpha=0.5, label="Conditional moments")

savefig(plt_sols, "plots/mcm_sim.pdf")
