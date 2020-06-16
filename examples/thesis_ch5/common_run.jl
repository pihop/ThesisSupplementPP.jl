using JLD2

function run_fluid(iter, u0_coupl, u0_iter, params)
    println("Fluid")
    println("Solving the direct coupling problem")
    time_cpl = @elapsed sol_cpl = run_single(f_coupl, u0_coupl, tspan, params, [], tplot)
    println("Time taken $time_cpl s")

    println("Solving the weighted ODE problem")
    time_iter = @elapsed sol_iter = run_iterative(f_iter, 
                                                  u0_iter, 
                                                  tspan, 
                                                  params, 
                                                  length(u0_iter), 
                                                  length(rns), 
                                                  w0, 
                                                  tplot)
    println("Time taken $time_iter s")

    str = string(iter, "_", Int(scaling))

    hitting_t = sol_iter[1]
    iter_save = [tuples(sol_iter[2][i]) for i in 1:length(rns)]
    cpl_save = tuples(sol_cpl)

    jldopen(joinpath("data/fluid_res_$str.jld2"), "w") do file
        file["tspan"] = tspan
        file["tplot"] = tplot
        file["pop"] = scaling
        file["params"] = params
        file["cpl/lhs"] = lhs_fluid_cpl
        file["cpl/initial"] = u0_coupl
        file["cpl/sol"] = cpl_save 
        file["cpl/time"] = time_cpl
        file["iter/hitting_t"] = hitting_t 
        file["iter/sol"] = iter_save
        file["iter/time"] = time_iter
        file["iter/initial"] = u0_iter
        file["iter/lhs"] = lhs_fluid_iter
    end
end

function run_lna(iter, u0_coupl, u0_iter, params)
    println("LNA")
    println("Solving the direct coupling problem")
    time_cpl = @elapsed sol_cpl = run_single(f_lna_coupl, u0_coupl, tspan, params, [], tplot)
    println("Time taken $time_cpl s")

    println("Solving the weighted ODE problem")
    time_iter = @elapsed sol_iter = run_iterative(f_lna_iter, 
                                                  u0_iter, 
                                                  tspan, 
                                                  params, 
                                                  length(u0_iter), 
                                                  length(rns), 
                                                  w0, 
                                                  tplot)
    println("Time taken $time_iter s")

    str = string(iter, "_", Int(scaling))
    path = joinpath("data/lna_cpl_$str.jld2")

    hitting_t = sol_iter[1]
    iter_save = [tuples(sol_iter[2][i]) for i in 1:length(rns)]
    cpl_save = tuples(sol_cpl)

    jldopen(joinpath("data/lna_res_$str.jld2"), "w") do file
        file["tspan"] = tspan
        file["tplot"] = tplot
        file["pop"] = scaling
        file["params"] = params
        file["cpl/initial"] = u0_coupl
        file["cpl/sol"] = cpl_save 
        file["cpl/time"] = time_cpl
        file["cpl/lhs"] = lhs_lna_coupl
        file["iter/hitting_t"] = hitting_t 
        file["iter/sol"] = iter_save
        file["iter/time"] = time_iter
        file["iter/initial"] = u0_iter
        file["iter/lhs"] = lhs_lna_iter 
    end
end

function run_mc(iter, u0_coupl, u0_iter, params)
    println("MC")
    println("Solving the direct coupling problem")
    time_cpl = @elapsed sol_cpl = run_single(f_mc_coupl, u0_coupl, tspan, params, [], tplot)
    println("Time taken $time_cpl s")

    println("Solving the weighted ODE problem")
    time_iter = @elapsed sol_iter = run_iterative(f_mc_iter, 
                                                  u0_iter, 
                                                  tspan, 
                                                  params, 
                                                  length(u0_iter), 
                                                  length(rns), 
                                                  w0, 
                                                  tplot)
    println("Time taken $time_iter s")

    str = string(iter, "_", Int(scaling))
    path = joinpath("data/mc_cpl_$str.jld2")

    hitting_t = sol_iter[1]
    iter_save = [tuples(sol_iter[2][i]) for i in 1:length(rns)]
    cpl_save = tuples(sol_cpl)

    jldopen(joinpath("data/mc_res_$str.jld2"), "w") do file
        file["tspan"] = tspan
        file["tplot"] = tplot
        file["pop"] = scaling
        file["params"] = params
        file["cpl/initial"] = u0_coupl
        file["cpl/sol"] = cpl_save 
        file["cpl/time"] = time_cpl
        file["cpl/lhs"] = lhs_mc_coupl
        file["iter/hitting_t"] = hitting_t 
        file["iter/sol"] = iter_save
        file["iter/time"] = time_iter
        file["iter/initial"] = u0_iter
        file["iter/lhs"] = lhs_mc_iter 
    end
end

