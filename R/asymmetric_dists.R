get_exact_asymmetric_dist = function(laet){
    probs_X_equals_one = laet$Y[laet$X == 1] # Y-label probabilities for i s.t. X_i = 1
    probs_X_equals_zero = laet$Y[laet$X == 0] # Y-label probabilities for i s.t. X_i = 0
    dist_X_equals_one = successive_paired_convolution(probs_X_equals_one)
    dist_X_equals_zero = successive_paired_convolution(probs_X_equals_zero)
    dist_all = capped_convolution(dist_X_equals_one, dist_X_equals_zero)
    pmf = rep(0, laet$N + 1)
    for (j in max(0, laet$k + laet$m - laet$N):min(laet$m, laet$k)){
        pmf[j + 1] = dist_X_equals_one[j + 1] * dist_X_equals_zero[laet$k - j + 1] / dist_all[laet$k + 1]
    }

    return(pmf)
}

get_saddlepoint_asymmetric_dist = function(laet){
    hess_fn = hess_asymmetric(laet$X, laet$Y)

    base_grad_fn = grad_asymmetric(laet$X, laet$Y, laet$k, 0.5)
    restricted_base_grad_fn = function(v) base_grad_fn(c(v, 0))[1:1]
    base_soln = pracma::broyden(restricted_base_grad_fn, c(0), matrix(hess_fn(c(0, 0))[1:1,1:1], nrow=1))$zero
    det_of_hess_at_base_soln = hess_fn(c(base_soln, 0))[1,1]#det(matrix(hess_fn(c(base_soln, 0))[1:1,1:1], nrow=1))

    K_at_base_soln = K_asymmetric(laet$X, laet$Y)(c(base_soln, 0))
    min_possible_value = max(0, laet$m + laet$k - laet$N)
    max_possible_value = min(laet$m, laet$k)
    cdf = matrix(rep(NaN, 2*laet$N), nrow=2) # populate results into here

    result = sapply(min_possible_value:max_possible_value,
        asymmetric_saddlepoint_cdf,
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