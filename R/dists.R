get_MC_dist = function(laet, iterations){
    if (is.null(iterations)){
        cat("error: MC simulation without a specified number of iterations")
        return()
    }
    # Run Monte Carlo simulation n times to generate empirical distribution
    # of P(Z = j | X = m, Y = k) for each j

    # x_convs[j] gives unconditional distribution of sum_{i > j} X_i for each j, same for y_convs
    x_convs = iterative_paired_convolution(laet$X, laet$m + 1)
    y_convs = iterative_paired_convolution(laet$Y, laet$k + 1)

    result = rep(0, laet$N + 1)

    for (i in 1:iterations){
        sampled_X = sample_from_probs_cond_on_total(laet$X, laet$m, x_convs)
        sampled_Y = sample_from_probs_cond_on_total(laet$Y, laet$k, y_convs)
        Z = sum(sampled_X * sampled_Y)
        result[Z + 1] = result[Z + 1] + 1
    }

    result / iterations
}

get_exact_dist = function(laet) {
    # establish dimensions of 3D array
    m_dim = laet$m + 1
    k_dim = laet$k + 1
    z_dim_without_shortcut = min(laet$m, laet$k) + 1
    z_dim = min(z_dim_without_shortcut, max(laet$o)+1)

    dist = array(0, c(m_dim, k_dim, z_dim))
    dist[1,1,1] = 1 # set observed value of empty sum to 0 with probability 1
    n = 0

    unconditional_X_probs = successive_paired_convolution(laet$X, max_len = laet$m + 1)
    unconditional_Y_probs = successive_paired_convolution(laet$Y, max_len = laet$k + 1)

    for(p_x in laet$X){
        n = n + 1
        p_y = laet$Y[n]
        if(n %% 10 == 0)
            print(n)

        m = min(m_dim, n + 1)
        k = min(k_dim, n + 1)
        z = min(z_dim, n + 1)

        # set up some indexes we will use a few times
        a0 = 1:m
        a1 = 1:(m - 1)
        b0 = 1:k
        b1 = 1:(k - 1)
        c0 = 1:z
        c1 = 1:(z - 1)

        new_dist = dist * (1 - p_x) * (1 - p_y)
        new_dist[a1+1,b0,c0] = new_dist[a1+1,b0,c0] + dist[a1,b0,c0] * p_x * (1 - p_y)
        new_dist[a0,b1+1,c0] = new_dist[a0,b1+1,c0] + dist[a0,b1,c0] * (1 - p_x) * p_y
        new_dist[a1+1,b1+1,c1+1] = new_dist[a1+1,b1+1,c1+1] + dist[a1,b1,c1] * p_x * p_y

        dist = new_dist
    }

    d = dist[laet$m+1,laet$k+1,] / unconditional_X_probs[laet$m + 1] / unconditional_Y_probs[laet$k + 1]

    if(z_dim < z_dim_without_shortcut){
        # we have taken a shortcut and skipped higher observed values.
        # just assign any missing mass to the end of the pmf so it doesn't interfere
        d = c(d, rep(0, z_dim_without_shortcut - z_dim - 1), 1 - sum(d))
    }
    return(d)
}

get_binom_dist = function(laet) successive_paired_convolution(laet$X*laet$Y)

get_saddlepoint_dist = function(laet) {
    hess_fn = hess_K(laet$X, laet$Y)
    
    base_grad_fn = grad_K(laet$m, laet$k, 0.5, laet$X, laet$Y)
    restricted_base_grad_fn = function(v) base_grad_fn(c(v, 0))[1:2]
    base_soln = pracma::broyden(restricted_base_grad_fn, c(0, 0), hess_fn(c(0, 0, 0))[1:2,1:2])$zero
    det_of_hess_at_base_soln = det(hess_fn(c(base_soln, 0))[1:2,1:2])

    K_at_base_soln = K(c(base_soln, 0), laet$X, laet$Y)
    
    min_possible_value = max(0, laet$m + laet$k - laet$N)
    max_possible_value = min(laet$m, laet$k, max(laet$o)+1)
    cdf = matrix(rep(NaN, laet$N*2), nrow=2) # populate results into here

    result = sapply(min_possible_value:max_possible_value,
        saddlepoint_cdf,
        laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln
    )

    if(min_possible_value == 0)
        cdf[1,1:max_possible_value] = result[1,2:length(result[1,])]
    else
        cdf[1,min_possible_value:max_possible_value] = result[1,]

    cdf[2,1+(min_possible_value:max_possible_value)] = result[2,]

    # saddlepoint struggles with the obvious endpoints
    cdf[1,max_possible_value+1] = 1.0
    cdf[2, 1] = 1.0
    return(cdf)
}