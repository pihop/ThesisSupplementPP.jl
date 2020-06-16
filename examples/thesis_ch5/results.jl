using ThesisSupplementPP
using DiffEqBase, OrdinaryDiffEq
using LinearAlgebra
using StatsBase
using JLD2
using FileIO

struct Compare
    iters
    scaling
    params
    # Time to generate the simulation trajectories only.
    sim_time
    # Includes the construction of simulation summary.
    sim_time_total
    sol_sim
    cpl_time
    iter_time
    sol_cpl
    sol_iter
    tspan
    sampling
end

function process(model_dir, scaling, nexp, approx="fluid")
    # nexp -- number of experiments
    sim_time = []
    sim_time_total = []
    sol_sim = []
    cpl_time = []
    iter_time = []
    sol_cpl = []
    sol_iter = []
    params = []
    sampling = 0
    tspan = 0 
    for i in 1:nexp 
        str = string(i, "_", Int(scaling))
        path = joinpath(@__DIR__, model_dir, "data/", "$(approx)_res_$str.jld2")
        jldopen(path, "r") do file
            # first is needed for the 1000 runs experiment...
            push!(cpl_time, file["cpl/time"])
            push!(iter_time, file["iter/time"])
            push!(sol_cpl, file["cpl/sol"])
            push!(sol_iter, file["iter/sol"])
            tspan = file["tspan"]
        end
        path = joinpath(@__DIR__, model_dir, "data/ssa_res_$str.jld2")
        jldopen(path, "r") do file
            push!(sim_time, file["sim_time"])
            push!(sim_time_total, file["total_time"])
            push!(params, file["params"])
            push!(sol_sim, file["sim_sum"])
            sampling = file["sampling"]
        end
    end
    return Compare(nexp, scaling, params, sim_time, sim_time_total, 
                        sol_sim, cpl_time, iter_time, sol_cpl, sol_iter, tspan, sampling)
end

function get_symbs(model_dir, scaling, iter, approx="fluid")
    str = string(iter, "_", Int(scaling))
    path = joinpath(@__DIR__, model_dir, "data/", "$(approx)_res_$str.jld2")
    file = jldopen(path, "r")
    out = file["cpl/lhs"]
    close(file)
    return out 
end

function get_variance(idx, symbs, len, method)
    if method == "lna"
        var = Symbol("c_$(idx)$(idx)")
        return findall(x->x==var, symbs)[1]
    elseif method == "mc"
        array = zeros(Int, len)
        array[idx] = 2 
        str = string(array...)
        var = Symbol("M_$str")
        return findall(x->x==var, symbs)[1]
    else
        return 0
    end
end

function rel_error(series_ref, series_test)
    return filter(x -> (!isnan(x) && !isinf(x)), abs.(series_ref - series_test) ./ series_ref)
end

function rel_error_unclean(series_ref, series_test)
    return abs.(series_ref - series_test) ./ series_ref
end

function prop_error(series_ref, series_test, scaling)
    return abs.(series_ref - series_test) ./ scaling
end

function rel_errors(test_set, idx)
    rels_cpl = []
    rels_iter = []
    sca = test_set.scaling
    for i in 1:test_set.iters
        s_iter = sca .* getindex.(getindex.(test_set.sol_iter[i][end],1), idx)
        s_cpl =  sca .* getindex.(getindex.(test_set.sol_cpl[i],1), idx)
        s_sim = getindex.(test_set.sol_sim[i].u.u, idx)
        rel_iter = rel_error_unclean(s_sim, s_iter) 
        rel_cpl = rel_error_unclean(s_sim, s_cpl) 
        push!(rels_cpl, rel_cpl)
        push!(rels_iter, rel_iter)
    end
    return rels_cpl, rels_iter
end

function prop_errors(test_set, idx)
    props_cpl = []
    props_iter = []
    sca = test_set.scaling
    for i in 1:test_set.iters
        s_sim = getindex.(test_set.sol_sim[i].u.u, idx)
        s_iter = sca .* getindex.(getindex.(test_set.sol_iter[i][end],1), idx)
        s_cpl =  sca .* getindex.(getindex.(test_set.sol_cpl[i],1), idx)

        prop_iter = prop_error(s_sim, s_iter, sca) 
        prop_cpl = prop_error(s_sim, s_cpl, sca) 
        push!(props_cpl, prop_cpl)
        push!(props_iter, prop_iter)
    end
    return props_cpl, props_iter
end

function prop_errors_var(test_set, idx1, idx2, method)
    props_cpl = []
    props_iter = []
    sca = test_set.scaling

    coeff = 1
    if method == "lna"
        coeff = sca^(-0.5)
    elseif method == "mc"
        coeff = 1.0
    end

    for i in 1:test_set.iters
        sol_cpl = getindex.(test_set.sol_cpl[i], 1)
        sol_iter = getindex.(test_set.sol_iter[i][end], 1)

        s_sim = sqrt.(getindex.(test_set.sol_sim[i].v.u, idx1))
        s_iter = coeff .* real(sqrt.(abs.(getindex.(sol_iter, idx2))))
        s_cpl =  coeff .* real(sqrt.(abs.(getindex.(sol_cpl, idx2))))

        prop_iter = prop_error(s_sim, s_iter, sca) 
        prop_cpl = prop_error(s_sim, s_cpl, sca) 
        push!(props_cpl, prop_cpl)
        push!(props_iter, prop_iter)
    end
    return props_cpl, props_iter
end

function timings(test_set)
    # First result removed because of JIT compilation.
    mean_cpl = mean(test_set.cpl_time[2:end])
    mean_iter = mean(test_set.iter_time[2:end])
    mean_sim = mean(test_set.sim_time[2:end])
    mean_sim_total = mean(test_set.sim_time_total[2:end])
    println("Mean timings for parameters population size $(test_set.scaling)")
    println("Coupling $mean_cpl \t Iterative $mean_iter \t Simulation $mean_sim ($mean_sim_total)")
#    return mean(test_set.cpl_time[2:end]), mean(test_set.iter_time[2:end]), mean(test_set.sim_time[2:end])
    return mean_cpl, mean_iter, mean_sim
end

function all_timings(dir, nexp)
    # Corresponds to Table 5.1 from the thesis.
    pops = [100, 200, 300, 500, 1000]
    for p in pops
        fluid = process(dir, p, nexp, "fluid")
        lna = process(dir, p, nexp, "lna")
        mc = process(dir, p, nexp, "mc")
        println("Timings for fluid")
        timings(fluid)
        println("Timings for LNA")
        timings(lna)
        println("Timings for MC")
        timings(mc)
        println("\n")
    end
end

function all_errors_means(dir, idx, nexp)
    # Corresponds to Table 5.2 in the thesis.
    pops = [100, 200, 300, 500, 1000]
    for p in pops
        test_set = process(dir, p, nexp, "fluid")
        scaling = test_set.scaling
        # Maximum proportional errors
        prop_errs_cpl, prop_errs_iter = prop_errors(test_set, idx)
        prop_cpl_max = maximum.(prop_errs_cpl)
        prop_iter_max = maximum.(prop_errs_iter)

        # Mean and standard deviation of maximum errors.
        prop_mean_cpl = mean(prop_cpl_max)*100
        prop_mean_iter = mean(prop_iter_max)*100
        prop_sd_cpl = StatsBase.std(prop_cpl_max)*100
        prop_sd_iter = StatsBase.std(prop_iter_max)*100
        # Maximum of the maximum errors.
        max_cpl = maximum(prop_cpl_max)*100
        max_iter = maximum(prop_iter_max)*100

        # Print the results
        println("Population size \t mean (sd) \t\t\t\t maximum")
        println("Direct coupling")
        println("$scaling \t\t $prop_mean_cpl ($prop_sd_cpl) % \t $max_cpl")
        println("Iterative")
        println("$scaling \t\t $prop_mean_iter ($prop_sd_iter) % \t $max_iter")
        println()
    end
end
#
function all_errors_vars(model_dir, idx, method, nexp)
    # In order to find the correct variances from the stored solutions it is
    # useful to have the model file at hand.
    include(string(model_dir,"/rn.jl"))

    pops = [100, 200, 300, 500, 1000]

    for p in pops
        test_set = process(model_dir, p, nexp, method)
        scaling = test_set.scaling
        symbols = get_symbs(model_dir, scaling, 1, method)
        # Find the index of the variance.
        idx_var = get_variance(idx, symbols, length(m_grid), method)

        prop_errs_cpl, prop_errs_iter = prop_errors_var(test_set, idx, idx_var, method)
        prop_cpl_max = maximum.(prop_errs_cpl)
        prop_iter_max = maximum.(prop_errs_iter)

        # Mean and standard deviation of maximum errors.
        prop_mean_cpl = mean(prop_cpl_max)*100
        prop_mean_iter = mean(prop_iter_max)*100
        prop_sd_cpl = StatsBase.std(prop_cpl_max)*100
        prop_sd_iter = StatsBase.std(prop_iter_max)*100
        # Maximum of the maximum errors.
        max_cpl = maximum(prop_cpl_max)*100
        max_iter = maximum(prop_iter_max)*100

        # Print the results
        println("Population size \t mean (sd) \t\t\t\t maximum")
        println("Direct coupling")
        println("$scaling \t\t $prop_mean_cpl ($prop_sd_cpl) % \t $max_cpl")
        println("Iterative")
        println("$scaling \t\t $prop_mean_iter ($prop_sd_iter) % \t $max_iter")
        println()
    end
end
