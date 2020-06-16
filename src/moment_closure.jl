"""

For reference see
"A general moment expansion method for stochastic kinetic models" by A Ale, et. al

This is a very partial Julia reimplementation of parts of the Means package for
python (https://github.com/theosysbio/means).

"""

@inline function partial_derivative(expression::Basic, var_symbs::Array{Basic, 1}, ctr_entry::Tuple)
    """
    Take partial derivatives. ctr_entry = [1,0,1] maps to the partial derivative
    of expression with respect to x_1 followed by x_3
    """
#    @time diff(expression, Tuple(var_symbs), Tuple([convert.(Basic,ctr_entry)]...))
    return diff(expression, Tuple(var_symbs), Tuple([convert.(Basic,ctr_entry)]...))
end

function one_over_factorial(counter_entry)
    # (n_1, ..., n_m) -> 1/(n_1!*...*n_m!)
    return Basic(1 / mapreduce(x -> factorial(x), *, counter_entry))
end

function vec_binomial(vec_n, vec_k)
    # Product of binomials
    return Basic(mapreduce(x -> binomial(x[1],x[2]), * , zip(vec_n, vec_k)))
end

function make_alpha(vars, vec_n, vec_k)
    return Basic(mapreduce(x -> x[1]^(x[2] - x[3]), *, zip(vars, vec_n, vec_k)))
end

function sign(vec_n, vec_k)
    return Basic(mapreduce(x -> (-1)^(x[1] - x[2]), *, zip(vec_n, vec_k)))
end

function solve(expr, var)
    # TODO: This is silly but for the examples in the thesis turns out to be ok.
    # However, an actual symbolic solver (SymPy.solve for example) is needed in
    # general. 
    # display(expr)
    res = -(expr - var)
    return res 
end

function gen_counters(max_order, var_symbs; central_symb_prefix="M_", raw_symb_prefix="x_")
    """
    Generate the moment descriptors up to (max_order + 1) moments.
    Nore that first order central moments (E[x_i - mu_i]) = 0.
    """
    n_moments = max_order + 1

    central_ctr = [Moment{length(var_symbs)}(ntuple(x -> 0, Val{length(var_symbs)}()), Basic(1))]
    raw_ctr = [Moment{length(var_symbs)}(ntuple(x -> 0, Val{length(var_symbs)}()), Basic(1))]
    
    # Descriptors for first order moment equations. 
    identity = Matrix(1I, length(var_symbs), length(var_symbs)) 
    units = collect(eachrow(identity))
    raw_ctr = vcat(raw_ctr, 
                   [Moment{length(var_symbs)}(tuple(identity[:,i]...), symb) 
                                 for (i, symb) in enumerate(var_symbs)])
    # Descriptors for higher order moment equations.
    higher_order_descriptors = vcat([sum.(collect(with_replacement_combinations(units,i))) for i in 2:n_moments]...)
    higher_order_descriptors = [ntuple(i->x[i], Val{length(var_symbs)}()) for x in higher_order_descriptors]

    higher_order_symbols = vec([symbols(string(raw_symb_prefix, collect(c for c in count)...)) 
                      for count in higher_order_descriptors])

    append!(raw_ctr, [Moment{length(var_symbs)}(higher_order_descriptors[i], symb) 
                        for (i, symb) in enumerate(higher_order_symbols)])

    central_ctr_descriptors = filter(x -> sum(x) > 1, higher_order_descriptors)

    central_ctr_symbols = [symbols(string(central_symb_prefix, 
                                          collect(c for c in count)...))
                               for count in central_ctr_descriptors]

    central_ctr = vcat(central_ctr, 
                       [Moment{length(var_symbs)}(central_ctr_descriptors[i], symb) 
                        for (i, symb) in enumerate(central_ctr_symbols)])

    return raw_ctr, central_ctr
end

function gen_dmu_dt(rn, stoich_mat, central_ctr)
    """
    Returns the matrix corresponding to truncated Taylor expansion of the first
    order raw moments.
    Rows correspond to the propensities and columns to the central moments up to
    order (max_order + 1).
    """
    syms = convert.(Basic, rn.syms)
    rate_exprs = convert.(Basic, rn.jump_rate_expr)
    derivatives = [partial_derivative(expr, syms, c.n_vector) 
                   for (expr, c) in Iterators.product(rate_exprs, central_ctr)]
    factorial_terms = transpose(reduce(hcat,
                                       fill([one_over_factorial(c.n_vector) for c in central_ctr],
                                            length(rate_exprs))))
    return stoich_mat * (factorial_terms .* derivatives)
end

function central_expr(central_ctr, central_itr, raw_ctr, dmu_dt, rn, stoich_mat)
    """
    The time derivative of a central moment given by central_itr in terms of
    raw and central moments.
    """ 
    syms = convert.(Basic, rn.syms)
    central_vec = central_itr.n_vector
    taylor_expr = Array{Basic}[]
    for raw in raw_ctr
        raw_vec = raw.n_vector
        coeff = vec_binomial(central_vec, raw_vec) * sign(central_vec, raw_vec) 

        alpha = make_alpha(syms, central_vec, raw_vec) 

        dalpha_dt = mapreduce(x -> (x[2]-x[3])/x[1] * alpha .* x[4], 
                              +, 
                              zip(syms, central_vec, raw_vec, eachrow(dmu_dt)))

        e_ctr = filter(x -> all(raw_vec .>= x.n_vector) && order(x) > 0, raw_ctr)

        #slow
        dbeta_dt = get_dbeta_dt(rn, stoich_mat, central_ctr, raw_vec, e_ctr)

        if length(e_ctr) == 0
            beta = Basic(1)
        else
            beta = raw.symbol
        end

        push!(taylor_expr, coeff .* (alpha .* dbeta_dt + beta .* dalpha_dt))
    end

    return expand.(sum(taylor_expr))
end

function make_central_exprs(central_ctr, raw_ctr, dmu_dt, max_order, rn, stoich_mat)
    """
    Collect the time-derivatives. 
    """
    central_moments = Array{Basic}[]
    
    @inbounds for central_itr in central_ctr
        if (order(central_itr) == 0  || order(central_itr) > max_order)
            continue
        end
        # The other terms are going to be 0 in the expansion.
        raw_lower = filter(x -> all(central_itr.n_vector .>= x.n_vector), raw_ctr)
        push!(central_moments, central_expr(central_ctr, central_itr, raw_lower, dmu_dt, rn, stoich_mat))
    end
    return central_moments
end


function get_dbeta_dt(rn, stoich_mat, central_ctr, raw_vec, e_ctr)
#    println("Dbeta is slow at the moment")
    syms = convert.(Basic, rn.syms)
    rate_exprs = convert.(Basic, rn.jump_rate_expr)

    if length(e_ctr) == 0
        return transpose(fill(Basic(0), 1, length(central_ctr)))
    end
    f_of_x_vec = [make_f_of_x(syms, raw_vec, ek.n_vector,  reac)
                      for (reac, ek) in Iterators.product(rate_exprs, e_ctr)]

    f_expectation_vec = transpose(hcat([make_f_expectation(f, syms, central_ctr)
                                        for f in f_of_x_vec]...))

    s_pow_e_vec = reduce(hcat, [make_s_pow_e(stoich_mat, react_idx, ek.n_vector) 
                                for (react_idx, ek) in Iterators.product(1:length(rate_exprs), e_ctr)])

    k_choose_e_vec = reduce(vcat, [fill([vec_binomial(raw_vec, ek.n_vector)], length(rate_exprs)) for ek in e_ctr])

    s_times_ke = k_choose_e_vec .* transpose(s_pow_e_vec)

    return mapreduce(x -> x[1] .* x[2], +, zip(eachrow(f_expectation_vec), s_times_ke))
end

function make_f_of_x(var_symbs, k_vec, e_vec, propensity)
    return mapreduce(x -> x[2]^(k_vec[x[1]] - e_vec[x[1]]), *, enumerate(var_symbs))* propensity
end

function make_f_expectation(expr, var_symbs, central_ctr)
    return  map(x -> 
                partial_derivative(expr, var_symbs, x.n_vector) * one_over_factorial(x.n_vector),
                central_ctr)
end

function make_s_pow_e(stoichiometry, react_idx, e_vec)
    return mapreduce(x -> stoichiometry[x[1], react_idx]^x[2], *, enumerate(e_vec))
end

function raw_to_central(central_ctr, var_symbs, raw_ctr)
    central_in_terms_of_raw = Basic[]
    for central_iter in central_ctr
        if order(central_iter) == 0
            continue
        end
        central_vec = central_iter.n_vector
        raw_lower = [raw for raw in raw_ctr if all(central_iter.n_vector .>= raw.n_vector)]
        n_choose_k_vec = map(x -> vec_binomial(central_vec, x.n_vector), raw_lower)
        minus_one_pow_n_minus_k = [sign(central_vec, k_vec.n_vector)  for k_vec in raw_lower] 
        alpha_vec = [make_alpha(var_symbs, central_vec, k_vec.n_vector) for k_vec in raw_lower]

        beta = [raw.symbol for raw in raw_lower]
        product = [(n * m * a * b) for (n, m, a, b) in zip(n_choose_k_vec, minus_one_pow_n_minus_k, alpha_vec, beta)]

        push!(central_in_terms_of_raw, sum(product))
    end
    return central_in_terms_of_raw
end


function substitute_raw_with_central(central_moments_exprs, central_from_raw_exprs, central_ctr, raw_ctr)
    # TODO: Substitutions are slow...
    positive_raw_mnt_symbs = [raw.symbol for raw in raw_ctr if order(raw) > 1]
    central_symbols = [central.symbol for central in central_ctr if order(central) > 1]

    eq_to_solve = [cfr - cs for (cs, cfr) in zip(central_symbols, central_from_raw_exprs)]
    solved_xs = [solve(eq,raw) for (eq, raw) in zip(eq_to_solve, positive_raw_mnt_symbs)]

    substitution_pairs = Dict(zip(positive_raw_mnt_symbs, solved_xs))
    filter!(k -> k.second != Any[], substitution_pairs) 
    map!(x -> subs(x, substitution_pairs), central_moments_exprs, central_moments_exprs)

    return central_moments_exprs
end

function  generate_mass_fluctuation_kinetics(central_moments, dmu_over_dt, central_ctr)
    central_moments_symbols = [n.symbol for n in central_ctr]
    mfk = [e for e in dmu_over_dt * central_moments_symbols] 
    return vcat(mfk, [(transpose(central_moments_symbols) * cm) for cm in central_moments]...)
end

function generate_problem_left_hand_side(central_ctr, raw_ctr, max_order)
    prob_moments_over_dt = [k.symbol for k in raw_ctr if order(k) == 1]
    append!(prob_moments_over_dt, [n.symbol for n in central_ctr if max_order >= order(n) > 1])
    return prob_moments_over_dt
end
