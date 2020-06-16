struct State
    u :: Array{Float64}
    t :: Float64
end

struct Moment{N}
    n_vector :: NTuple{N, Int}
    symbol :: Basic
end

function order(m::Moment)
    return sum(m.n_vector)
end

function is_mixed(m::Moment)
    return !(order(m) in m.n_vector)
end

function replace_names(expr, old_names, new_names)
    mapping = Dict(zip(old_names, new_names))
    MacroTools.postwalk( x -> x in old_names ? x= mapping[x] : x, expr)
end

function make_symbolic_mc(rn, max_order::Int)
    # Generate counters for central and raw moments.
    raw_ctr, central_ctr = gen_counters(max_order, rn.syms)
    # This should be a matrix with each row a reaction and each column
    # corresponding to element counter.  
    # (-2, 2 ) * ( <a_2(x)> ) \approx (-2, 2 ) * ( taylor series summands )
    # (1 , -1)   ( <a_1(x)> )         (1 , -1)   ( taylor series summands )
    mfk = gen_moment_expressions(rn, central_ctr, raw_ctr, max_order) 
    mfk = close_normal(mfk[1], mfk[2], central_ctr, raw_ctr, max_order)
    lhs = generate_problem_left_hand_side(central_ctr, raw_ctr, max_order)
    return lhs, mfk
end

function gen_moment_expressions(rn, central_ctr, raw_ctr, max_order)
    stoich_mat = get_stoich_mat(rn) 
    dmu_dt = gen_dmu_dt(rn, stoich_mat, central_ctr)
    println("Making central expressions")
    @time central_mnt_exprs = make_central_exprs(central_ctr, 
                                                 raw_ctr, 
                                                 dmu_dt,
                                                 max_order,
                                                 rn,
                                                 stoich_mat)
    
    println("Raw to central")
    @time central_from_raw_exprs = raw_to_central(central_ctr, 
                                                  convert.(Basic, rn.syms), 
                                                  raw_ctr)

    println("Substitusions")
    @warn "This makes use of a solver function that in general is not correct.
    For the examples in the thesis it is ok though."
    
    @time map!((x) -> substitute_raw_with_central(x, 
                                                  central_from_raw_exprs, 
                                                  central_ctr, 
                                                  raw_ctr), central_mnt_exprs, central_mnt_exprs) 

#   @time central_mnt_exprs = map((x) -> substitute_raw_with_central(x, 
#                                                  central_from_raw_exprs, 
#                                                  central_ctr, 
#                                                  raw_ctr), central_mnt_exprs) 
#
    
    println("Kinetic expressions")
    @time mfk = generate_mass_fluctuation_kinetics(central_mnt_exprs, 
                                             dmu_dt, 
                                             central_ctr)
    
    return expand.(mfk), central_from_raw_exprs
end

function sym_to_expr(exprs, vars, params, params_t)
    new_vars = [:(u[$i]) for i in 1:length(vars)]
    new_params = [:(p[1][$i]) for i in 1:length(params)]
    new_params_t = [:(p[2](t)[$i]) for i in 1:length(params_t)]
    exprs = [DiffEqBiological.replace_names(convert(Expr, expr), vars, new_vars) for expr in exprs]
    exprs = [DiffEqBiological.replace_names(convert(Expr, expr), params, new_params) for expr in exprs]
    exprs = [DiffEqBiological.replace_names(convert(Expr, expr), params_t, new_params_t) for expr in exprs]
    return Expr(:tuple, exprs...)
end

function get_stoich(rn, idx)
    v = State(zeros(length(rn.syms)), 0.)
    (rn.jumps[idx]).affect!(v)
    return v.u
end

function get_prop(rn, idx)
    return convert(Basic, rn.jump_rate_expr[idx])
end

function get_stoich_mat(rn)
    stoich_mat = [get_stoich(rn, i) for i in 1:length(rn.jump_rate_expr)]
    return reshape(collect(Iterators.flatten(stoich_mat)), :, length(rn.jump_rate_expr))
end

function gen_diffusion(rn)
    changes = [(get_stoich(rn, idx), get_prop(rn,idx)) for idx in 1:length(rn.jumps)]
    return sum([(vec[1]) * (transpose(vec[1])) * vec[2] for vec in changes])
end

function gen_noise_covar(rn, lhs)
    jac = jacobian(rn.f_symfuncs, rn.syms)
    lhs = convert.(Basic, lhs)
    return unique(expand.(jac*lhs + lhs*transpose(jac) + gen_diffusion(rn)))
end

function lower_first(i, j)
    if i <= j
        return Symbol("c_$i$j")
    else 
        return Symbol("c_$j$i")
    end
end

function gen_covar_syms(rn)
    lhs = convert.(Basic, rn.syms)
    return [lower_first(i, j) for i in 1:length(lhs), j in 1:length(lhs)]
end

function jacobian(exprs, vars)
    jac = [[diff(expr, var) for var in vars] for expr in exprs]
    return transpose(hcat(jac...))
end
