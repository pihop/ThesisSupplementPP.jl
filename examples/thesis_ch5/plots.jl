using Plots
using JLD2
using FileIO
using DiffEqBase, OrdinaryDiffEq
using LaTeXStrings
pyplot()
PyPlot.rc("text", usetex="true")
fonts = font(12)

include("results.jl")

function permute_loc(grid)
    dim = size(grid)[1]
    out = [grid[i, dim+1-j] for i=1:dim, j=1:dim]
    return permutedims(out, [2,1])
end

function permute_loc(grid)
    dim = size(grid)[1]
    out = [grid[i, dim+1-j] for i=1:dim, j=1:dim]
    return permutedims(out, [2,1])
end


# function make_comp_plot(model_dir, test_set, iter, locs) 
#     # Working model
#     include(string(model_dir,"/rn.jl"))
#     vls = convert.(Tuple, permute_loc(m_grid))
#     loc_map = Dict(convert.(Tuple, m_grid) .=> vls)
# 
#     scaling = test_set.scaling
#     str = string(iter, "_", scaling)
#     sim_sum = test_set.sol_sim[iter]
#     sol_cpl = getindex.(test_set.sol_cpl[iter], 1)
#     sol_iter = getindex.(test_set.sol_iter[iter][end], 1)
#     tplot = sim_sum.t
#     plt_sols = plot(ylabel=latexstring("\\mbox{Population at a given location}"), 
#                    xlabel=latexstring("\\mbox{Time}"), 
#                    xlim=tplot,
#                    ylim=(0.0, sum(sim_sum[1])),
#                    legend=:topright, 
#                    fg_legend=:transparent,
#                    bg_legend=:transparent,
#                    size=(400,300),
#                    xtickfont=fonts, 
#                    ytickfont=fonts,
#                    guidefont=fonts, 
#                    legendfont=fonts)
# 
#     for (i, loc) in enumerate(locs)
#         mean_sim = getindex.(sim_sum.u.u, L[loc...])
#         upper = mean_sim .+ sqrt.(getindex.(sim_sum.v.u, L[loc...]))
#         lower = max.(0.0, mean_sim .- sqrt.(getindex.(sim_sum.v.u, L[loc...])))
#         mean_iter = scaling*getindex.(sol_iter, L[loc...])
#         mean_cpl = scaling*getindex.(sol_cpl, L[loc...])
# 
#         plt_sols = scatter!([-1], [-1], color=i, 
#                             label=latexstring("\\mbox{Location }", string(loc_map[loc])))
#         plt_sols = plot!(tplot, upper; color=i, lw=1.1, line=(:solid), labels="")
#         plt_sols = plot!(tplot, lower; color=i, lw=1.1, line=(:solid), labels="")
# 
#         plt_sols = plot!(tplot, mean_iter; color=i, lw=1.1, line=(:dashdot), labels="")
#         plt_sols = plot!(tplot, mean_cpl; color=i, lw=1.1, line=(:dash), labels="")
#     end
#     plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:dashdot, 1.0), label="Iterative")
#     plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:dash, 1.0), label="Direct coupling")
#     plt_sols = plot!(tplot, t -> -1.0, color=:grey, line=(:solid, 1.5), linealpha=0.5, label="SSA")
#     savefig(plt_sols, string(model_dir, "/plots/ssa_$str.pdf"))
# end
# 
# function make_err_plot(mode_dir, test_set, iter, locs)
#     # Working model
#     include(string(model_dir,"/rn.jl"))
#     vls = convert.(Tuple, permute_loc(m_grid))
#     loc_map = Dict(convert.(Tuple, m_grid) .=> vls)
# 
#     scaling = test_set.scaling
#     str = string(iter, "_", scaling)
#     plt_rel = plot(ylabel=latexstring("\\mbox{Relative error for location } (1,1)"), 
#                    xlabel=latexstring("\\mbox{Time}"), 
#                legend=:right, 
#                fg_legend=:transparent,
#                bg_legend=:transparent,
#                size=(700,300),
#                xtickfont=fonts, 
#                ytickfont=fonts,
#                guidefont=fonts, 
#                legendfont=fonts,
#                ylims=(0.0, 1.0)
#               )
#     sim_sum = test_set.sol_sim[iter]
#     sol_cpl = getindex.(test_set.sol_cpl[iter], 1)
#     sol_iter = getindex.(test_set.sol_iter[iter][end], 1)
#     tplot = sim_sum.t
# 
#     for (i, loc) in enumerate(locs)
#         mean_sim = getindex.(sim_sum.u.u, L[loc...])
#         mean_iter = scaling*getindex.(sol_iter, L[loc...])
#         mean_cpl = scaling*getindex.(sol_cpl, L[loc...])
#         plt_rel = plot!(tplot, rel_error_unclean(mean_sim, mean_iter); fillalpha=0.2, lw=1.5, line=(:solid), color=i, labels=latexstring("\\mbox{Iterative}"))
# 
#         plt_rel = plot!(tplot, rel_error_unclean(mean_sim, mean_cpl); fillalpha=0.2, lw=1.5, line=(:dash), color=i, labels=latexstring("\\mbox{Direct coupling}"))
#     end
#     savefig(plt_rel, string(model_dir, "/plots/rel_$str.pdf"))
# end
# 
# function make_prop_plot(model_dir, test_set, iter, locs)
#     # Working model
#     include(string(model_dir,"/rn.jl"))
#     vls = convert.(Tuple, permute_loc(m_grid))
#     loc_map = Dict(convert.(Tuple, m_grid) .=> vls)
# 
#     scaling = test_set.scaling
#     str = string(iter, "_", scaling)
#     xticks = 
#     plt_prop = plot(ylabel="\$\\mbox{Relative error for location } (1,1)\$", 
#                    xlabel="\$\\mbox{Time}\$", 
#                legend=:right, 
#                fg_legend=:transparent,
#                bg_legend=:transparent,
#                size=(700,300),
#                xtickfontsize=10, 
#                ytickfontsize=10,
#                xguidefontsize=10, 
#                legendfontsize=10,
#                ylims=(0.0, 0.2)
#               )
#     sim_sum = test_set.sol_sim[iter]
#     sol_cpl = getindex.(test_set.sol_cpl[iter], 1)
#     sol_iter = getindex.(test_set.sol_iter[iter][end], 1)
#     tplot = sim_sum.t
# 
#     for (i, loc) in enumerate(locs)
#         mean_sim = getindex.(sim_sum.u.u, L[loc...])
#         mean_iter = scaling .* getindex.(sol_iter, L[loc...])
#         mean_cpl = scaling .* getindex.(sol_cpl, L[loc...])
#         plt_prop = plot!(tplot, prop_error(mean_sim, mean_iter, scaling); fillalpha=0.2, lw=1.5, line=(:solid), color=i, labels="\$\\mbox{Iterative}\$")
# 
#         plt_prop = plot!(tplot, prop_error(mean_sim, mean_cpl, scaling); fillalpha=0.2, lw=1.5, line=(:dash), color=i, labels="\$\\mbox{Direct coupling}\$")
#     end
#     savefig(plt_prop, string(model_dir, "/plots/prop_$str.pdf"))
# end
# 
# function make_cerr_plot(model_dir, test_set, iter, locs)
#     # Working model
#     include(string(model_dir,"/rn.jl"))
#     vls = convert.(Tuple, permute_loc(m_grid))
#     loc_map = Dict(convert.(Tuple, m_grid) .=> vls)
# 
#     scaling = test_set.scaling
#     str = string(iter, "_", scaling)
#     plt_cerr = plot(ylabel=latexstring("\\mbox{Error for location } (1,1)"), 
#                    xlabel=latexstring("\\mbox{Time}"), 
#                legend=:right, 
#                fg_legend=:transparent,
#                bg_legend=:transparent,
#                size=(400,300),
#                xtickfont=fonts, 
#                ytickfont=fonts,
#                guidefont=fonts, 
#                legendfont=fonts,
#                ylims=(0.0, 0.5)
#               )
#     sim_sum = test_set.sol_sim[iter]
#     sol_cpl = getindex.(test_set.sol_cpl[iter], 1)
#     sol_iter = getindex.(test_set.sol_iter[iter][end], 1)
#     tplot = sim_sum.t
# 
#     for (i, loc) in enumerate(locs)
#         mean_sim = getindex.(sim_sum.u.u, L[loc...])
#         mean_iter = scaling*getindex.(sol_iter, L[loc...])
#         mean_cpl = scaling*getindex.(sol_cpl, L[loc...])
#         plt_cerr = plot!(tplot, rel_error_unclean(mean_sim, mean_iter); fillalpha=0.2, lw=1.5, line=(:solid), color=1, labels=latexstring("\\mbox{Iterative (relative)}"))
# 
#         plt_cerr = plot!(tplot, rel_error_unclean(mean_sim, mean_cpl); fillalpha=0.2, lw=1.5, line=(:dash), color=1, labels=latexstring("\\mbox{Direct coupling (relative)}"))
# 
#         plt_cerr = plot!(tplot, prop_error(mean_sim, mean_iter, scaling); fillalpha=0.2, lw=1.5, line=(:solid), color=2, labels=latexstring("\\mbox{Iterative (proportional)}"))
# 
#         plt_cerr = plot!(tplot, prop_error(mean_sim, mean_cpl, scaling); fillalpha=0.2, lw=1.5, line=(:dash), color=2, labels=latexstring("\\mbox{Direct coupling (proportional)}"))
#     end
#     savefig(plt_cerr, string(model_dir, "/plots/cerr_$str.pdf"))
# end
# 
# 
function make_heat_plot(model_dir, scaling, idx_sim, nexp)
    # Working model
    include(string(model_dir,"/rn.jl"))
    test_set = process(model_dir, scaling, nexp)

    p_errs_cpl, p_errs_iter = prop_errors(test_set, idx_sim)
    p_cpl_max = maximum.(p_errs_cpl)
    p_iter_max = maximum.(p_errs_iter)
    xs = getindex.(test_set.params, 1)
    ys = getindex.(test_set.params, 2)

    plt_heat_cpl = plot(ylabel=L"r_s", xlabel=L"r_m", zlabel="error", 
                        size=(600,500),
                        xtickfont=font(16), 
                        ytickfont=font(16),
                        ztickfont=font(16),
                        guidefont=font(16), 
                        legendfont=font(16),
                        ylims=(0.0, 1.0),
                        xlims=(0.0, 1.0),
                        zlims=(0.0, 0.4),
                        clims=(0.,0.4))
    plt_heat_cpl = plot!(xs, ys, p_cpl_max; st=:surface, color=:blues,)
    savefig(plt_heat_cpl, string(model_dir, "/plots/heat_cpl_$(string(scaling)).pdf"))

    plt_heat_iter = plot(ylabel=L"r_s", xlabel=L"r_m", zlabel="error", 
                         size=(600,500),
                         xtickfont=font(16), 
                         ytickfont=font(16),
                         ztickfont=font(16), 
                         guidefont=font(16), 
                         legendfont=font(16),
                         ylims=(0.0, 1.0),
                         xlims=(0.0, 1.0),
                         zlims=(0.0, 0.4),
                         clims=(0.,0.4))
    plt_heat_iter = plot!(xs, ys, p_iter_max; st=:surface, color=:blues,)

    savefig(plt_heat_iter, string(model_dir, "/plots/heat_iter_$(string(scaling)).pdf"))
end

function make_heat_plot_var(model_dir, scaling, idx_sim, method, nexp)
    # Working model
    include(string(model_dir,"/rn.jl"))
    test_set = process(model_dir, scaling, nexp, method)

    symbols = get_symbs(model_dir, scaling, 1, method)
    idx_var = get_variance(idx_sim, symbols, length(m_grid), method)
    p_errs_cpl, p_errs_iter = prop_errors_var(test_set, idx_sim, idx_var, method)
    p_cpl_max = maximum.(p_errs_cpl)
    p_iter_max = maximum.(p_errs_iter)

    # Parameter values on x and y axis.
    xs = getindex.(test_set.params, 1)
    ys = getindex.(test_set.params, 2)

    plt_heat_cpl = plot(ylabel=L"r_s", xlabel=L"r_m", zlabel="error", 
                        size=(600,500),
                        xtickfont=font(16), 
                        ytickfont=font(16),
                        ztickfont=font(16),
                        guidefont=font(16), 
                        legendfont=font(16),
                        ylims=(0.0, 1.0),
                        xlims=(0.0, 1.0),
                        zlims=(0.0, 0.4),
                        clims=(0.,0.4))
    plt_heat_cpl = plot!(xs, ys, p_cpl_max; st=:surface, color=:blues,)
    savefig(plt_heat_cpl, string(model_dir, "/plots/var_heat_cpl_$(string(scaling)).pdf"))

    plt_heat_iter = plot(ylabel=L"r_s", xlabel=L"r_m", zlabel="error", 
                         size=(600,500),
                         xtickfont=font(16), 
                         ytickfont=font(16),
                         ztickfont=font(16), 
                         guidefont=font(16), 
                         legendfont=font(16),
                         ylims=(0.0, 1.0),
                         xlims=(0.0, 1.0),
                         zlims=(0.0, 0.4),
                         clims=(0.,0.4))

    plt_heat_iter = plot!(xs, ys, p_iter_max; st=:surface, color=:blues,)
    savefig(plt_heat_iter, string(model_dir, "/plots/var_heat_iter_$(string(scaling)).pdf"))
end

function make_var_plot(model_dir, scaling, iter, locs, method, nexp; ylim=(0.0, 1.0), kwargs...)
    # Model file
    include(string(model_dir,"/rn.jl"))

    test_set = process(model_dir, scaling, nexp, method)
    tspan = test_set.tspan
#    println("Paramters of the simulation were ", test_set.params[iter])

    sim_sum = test_set.sol_sim[iter]
    sol_cpl = getindex.(test_set.sol_cpl[iter], 1)
    sol_iter = getindex.(test_set.sol_iter[iter][end], 1)
    tplot = sim_sum.t
    plt_vars = plot(ylabel=latexstring("\\mbox{Standard deviation for a given location}"), 
                    xlabel=latexstring("\\mbox{Time}"), 
                    xlim=tspan,
                    ylim=ylim,
                    legend=:topright, 
                    fg_legend=:transparent,
                    bg_legend=:transparent,
                    size=(400,300),
                    xtickfont=fonts, 
                    ytickfont=fonts,
                    guidefont=fonts, 
                    legendfont=fonts)
    coeff = 1.0
    if method == "lna"
        coeff = scaling^(-0.5)
    elseif method == "mc"
        coeff = scaling^(-1) 
    end

    for (i, loc) in enumerate(locs)
        symbols = get_symbs(model_dir, scaling, iter, method)
        var_idx = get_variance(L[loc...], symbols, length(m_grid), method)
        # Correct scaling of the results.
        var_iter = coeff .* sqrt.(abs.(getindex.(sol_iter, var_idx)))
        var_cpl =  coeff .* sqrt.(abs.(getindex.(sol_cpl, var_idx)))
        var_sim =  sqrt.(getindex.(sim_sum.v.u, L[loc...])) / scaling

        plt_vars = plot!(tplot, var_iter; color=i, lw=1.1, line=(:dashdot), labels="Iterative")
        plt_vars = plot!(tplot, var_cpl; color=i, lw=1.1, line=(0.5, :dash), labels="Direct coupling")
        plt_vars = plot!(tplot, var_sim; color=i, lw=1.1, line=(:solid), labels="Stochastic Simulation")
    end
    str = string(method, "_", iter, "_", scaling)
    savefig(plt_vars, string(model_dir, "/plots/var_$str.pdf"))
    open(string(model_dir, "/plots/var_$str.txt"), "w") do io
        write(io, string(test_set.params[iter]))
    end
end

