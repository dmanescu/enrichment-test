K = function(v, p_x, p_y){
    # returns the value of the cumulant generating function of
    # (X ~ Ber(p_x), Y ~ Ber(p_y), Z = X * Y). This is denoted
    # K(s_0, s_1, t), and is evaluated at v = (s_0, s_1, t)
    ev = exp(v)
    es_0 = ev[1]
    es_1 = ev[2]
    et = ev[3]
    MGF_terms = get_MGF_terms(p_x, p_y, es_0, es_1, et)
    return(sum(log(apply(MGF_terms, 2, sum))))
}

grad_K = function(m, k, j, p_x, p_y){
    # function to calculate grad K(s_0, s_1, t)
    # the result is then shifted by (m, k, j - 0.5)
    # so that when this function is plugged into a root solver
    # the solver returns the values (s_0, s_1, t) for which
    # grad K(s_0, s_1, t) = (m, k, j - 0.5)
    # which is the required form for the saddlepoint solver
    function(v, x = p_x, y = p_y){
        ev = exp(v)
        es_0 = ev[1]
        es_1 = ev[2]
        et = ev[3]
        MGF_terms = get_MGF_terms(x, y, es_0, es_1, et)
        return(get_grad(MGF_terms) - c(m, k, j - 0.5))
    }
}

hess_K = function(p_x, p_y){
    # function to calculate hess K(s_0, s_1, t) which can
    # be fed into the root solver and is also used in the
    # saddlepoint p-value formula
    function(v){
        ev = exp(v)
        es_0 = ev[1]
        es_1 = ev[2]
        et = ev[3]
        MGF_terms = get_MGF_terms(p_x, p_y, es_0, es_1, et)
        result = matrix(0, nrow=3, ncol=3)
        result[1,] = get_hess_first_row(MGF_terms)
        result[2,] = get_hess_second_row(MGF_terms)
        result[3,] = get_hess_third_row(MGF_terms)
        return(result)
    }
}


get_MGF_terms = function(p_x, p_y, es_0, es_1, et){
    # returns a 4xN matrix L containing all the summands of
    # the joint MGF/CGF of (X, Y, Z). Recall that the MGF will be given by
    # M(s_0, s_1, t) = \sum_i A[i] + B[i] + C[i] + D[i]
    # and
    # K(s_0, s_1, t) = \sum_i log(A[i] + B[i] + C[i] + D[i])
    # where A, B, C and D are the four arrays which form the rows
    # of this matrix:
    result = matrix(0, nrow=4, ncol=length(p_x))
    result[1,] = p_x * p_y * es_0 * es_1 * et
    result[2,] = p_x * (1 - p_y) * es_0
    result[3,] = (1 - p_x) * p_y * es_1
    result[4,] = (1 - p_x) * (1 - p_y)
    result
}

get_hess_first_row = function(MGF){
    # given all the summands evaluated by get_MGF_terms at the point (s_0, s_1, t)
    # this function returns the first row of the hessian matrix at that point.
    A = MGF[1,]
    B = MGF[2,]
    C = MGF[3,]
    D = MGF[4,]
    total = (A + B + C + D) ^ 2
    c(sum((A + B) * (C + D) / total),
      sum((A*D - B*C) / total),
      sum(A*(C + D) / total))
}

get_hess_second_row = function(MGF){
    # given all the summands evaluated by get_MGF_terms at the point (s_0, s_1, t)
    # this function returns the second row of the hessian matrix at that point.
    A = MGF[1,]
    B = MGF[2,]
    C = MGF[3,]
    D = MGF[4,]
    total = (A + B + C + D) ^ 2
    c(sum((A*D - B*C) / total),
      sum((A + C) * (B + D) / total),
      sum(A * (B + D) / total))
}

get_hess_third_row = function(MGF){
    # given all the summands evaluated by get_MGF_terms at the point (s_0, s_1, t)
    # this function returns the third row of the hessian matrix at that point.
    A = MGF[1,]
    B = MGF[2,]
    C = MGF[3,]
    D = MGF[4,]
    total = (A + B + C + D) ^ 2
    c(sum(A * (C + D) / total),
      sum(A * (B + D) / total),
      sum(A * (B + C + D) / total))
}

get_grad = function(MGF){
    # given all the summands evaluated by get_MGF_terms at the point (s_0, s_1, t)
    # this function returns the grad of the CGF evaluated at that point
    A = MGF[1,]
    B = MGF[2,]
    C = MGF[3,]
    D = MGF[4,]
    total = A + B + C + D
    result = matrix(0, nrow=3, ncol=length(A))
    c(sum((A + B) / total),
      sum((A + C) / total),
      sum(A / total))
}

saddlepoint_cdf = function(j, laet, hess_fn,
                            base_soln, det_of_hess_at_base_soln, K_at_base_soln){
    # function to return the saddlepoint cdf approximation at a single point j.
    # handles the case of a singularity at j by averaging the results 0.01 either side of j.

    grad_fn = grad_K(laet$m, laet$k, j, laet$X, laet$Y)

    soln = tryCatch({pracma::broyden(grad_fn, c(1, 1, 1), hess_fn(c(1, 1, 1)))},
                    error=function(err) NaN)

    if(!is.list(soln)){
        if (j != floor(j)) {
            # we were already at what we thought was a removable singularity, let's return NaN
            return(c(NaN, NaN))
        }
        else
            # (removable) singularity just perturb
            return(saddlepoint_cdf(j + 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln) * 0.5 +
                   saddlepoint_cdf(j - 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln) * 0.5)
    }

    soln = soln$zero
    sqrt_term_u2 = det(hess_fn(soln)) / det_of_hess_at_base_soln
    sqrt_term_w2 = (K_at_base_soln +
                    laet$m * (soln[1] - base_soln[1]) +
                    laet$k * (soln[2] - base_soln[2]) +
                    soln[3] * (j - 0.5) -
                    K(soln, laet$X, laet$Y))

    if (min(sqrt_term_u2, sqrt_term_w2) < 1e-8){
        if (j != floor(j))
            # we were already at what we thought was a singularity, let's return NaN
            return(c(NaN, NaN))
        else
            # (removable) singularity just perturb
            return(saddlepoint_cdf(j + 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln) * 0.5 +
                   saddlepoint_cdf(j - 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln) * 0.5)
    }

    u_2 = 2 * sinh(soln[3] / 2.0) * sqrt(sqrt_term_u2)
    w_2 = sign(soln[3]) * sqrt(2.0) * sqrt(sqrt_term_w2)

    if(is.nan(u_2)){
        if (j != floor(j))
            # we were already at what we thought was a singularity, let's return NaN
            return(c(NaN, NaN))
        else
            # (removable) singularity just perturb
            return(saddlepoint_cdf(j + 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln) * 0.5 +
                   saddlepoint_cdf(j - 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln) * 0.5)
    }

    return(c(pnorm(w_2) + dnorm(w_2)*(1/w_2 - 1/u_2),
             pnorm(-w_2) - dnorm(w_2)*(1/w_2 - 1/u_2)
             ))
}