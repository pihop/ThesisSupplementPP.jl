"""
For reference see
"A general moment expansion method forstochastic kinetic models" by A Ale, et. al

This is a very partial Julia reimplementation of parts of the Means package for
python (https://github.com/theosysbio/means).
In particular,
https://github.com/theosysbio/means/blob/master/src/means/approximation/mea/closure_normal.py

"""

function get_covariance_symbol(q_counter, sp1_idx, sp2_idx, multivariate=true)
    # The diagonal positions in the matrix are the variances
    if sp1_idx == sp2_idx
        return [q.symbol for q in q_counter if q.n_vector[sp1_idx] == 2 && order(q) == 2][1]
    # In multivariate cases, return covariances
    elseif multivariate
        return [q.symbol for q in q_counter if q.n_vector[sp1_idx] == 1 && q.n_vector[sp2_idx] == 1 && order(q) == 2][1]
    # In univariate cases, covariances are 0s
    else
        return Basic(0)
    end
end


function compute_one_closed_central_moment(moment, covariance_matrix)
    # If moment order is odd, higher order moments equals 0
    if order(moment) % 2 != 0
        return Basic(0)
    end

    # index of species
    idx = [i for i in 1:length(moment.n_vector)]

    # repeat the index of a species as many time as its value in counter
    list_for_partition = reduce(vcat, collect([Iterators.repeat([i], k) for (i,k) in zip(idx, moment.n_vector)]))

    if order(moment) == 2
        return covariance_matrix[list_for_partition[1], list_for_partition[2]]
    else
        error("TODO: Only setting higher order moments to 0 has been implemented.")
    end
end


function compute_closed_central_moments(central_from_raw_exprs, n_counter, k_counter)
    n_species = length([0 for pm in k_counter if order(pm) == 1])

    covariance_matrix = [get_covariance_symbol(n_counter, x, y) for (x,y) 
                         in Iterators.product(1:n_species, 1:n_species)]

    positive_n_counter = [n for n in n_counter if order(n) > 1]
    out_mat = [compute_one_closed_central_moment(n, covariance_matrix) for n in positive_n_counter]

    return out_mat
end

function close_normal(mfk, central_from_raw_exprs, n_counter, k_counter, max_order)
    # we obtain expressions for central moments in terms of variances/covariances

    closed_central_moments = compute_closed_central_moments(central_from_raw_exprs, n_counter, k_counter)

    # set mixed central moment to zero iff univariate
    closed_central_moments = set_mixed_moments_to_zero(true, closed_central_moments, n_counter)

    positive_n_counter = [n for n in n_counter if order(n) > 0]
    substitutions_pairs = [(n.symbol, ccm) for (n,ccm) in
                           zip(positive_n_counter, closed_central_moments) if order(n) > max_order]
    new_mfk = [subs(expr, substitutions_pairs...) for expr in mfk]
    return expand.(new_mfk)
end

function set_mixed_moments_to_zero(multivariate, closed_central_moments, n_counter)
    positive_n_counter = [n for n in n_counter if order(n) > 1]
    if multivariate
        return closed_central_moments
    else
        return [mixed_help(n, ccm) for (n,ccm) in zip(positive_n_counter, closed_central_moments)]
    end
end

function mixed_help(counter, ccm)
    if is_mixed(counter)
        return 0
    else 
        return ccm
    end
end
