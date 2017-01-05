hess_asymmetric = function(X, Y){
  function(v){
    ev = exp(v)
    er = ev[1]
    et = ev[2]
    ert = exp(v[1] * v[2])
    ert_or_er = ert * X + er * (1 - X)

    d2K_dr2 = sum(ert_or_er * (1 - Y) / (ert_or_er * Y + 1 - Y)**2)
    d2K_drdt = sum(ert_or_er * (1 - Y) * X / (ert_or_er * Y + 1 - Y)**2)
    d2K_dt2 = sum(ert * (1 - Y) * X / (ert_or_er * Y + 1 - Y)**2)
    matrix(c(
        d2K_dr2, d2K_drdt,
        d2K_drdt, d2K_dt2), nrow=2)
  }
}

grad_asymmetric = function(X, Y, k, j){
  # returns grad function shifted by (k, j - 0.5)
  offset = -c(k, j - 0.5)
  function(v){
    ev = exp(v)
    er = ev[1]
    et = ev[2]
    ert = exp(v[1] * v[2])
    ert_or_er = ert * X + er * (1 - X)
    
    c(sum(ert_or_er * Y / (ert_or_er * Y + 1 - Y)),
             sum(ert * Y * X / (ert_or_er * Y + 1 - Y))) + offset
  }
}

K_asymmetric = function(X, Y) {
  function(v){
    ev = exp(v)
    er = ev[1]
    et = ev[2]
    ert = exp(v[1] * v[2])
    ert_or_er = ert * X + er * (1 - X)
    sum(log(ert_or_er * Y + (1 - Y)))
  }  
}

asymmetric_saddlepoint_cdf = function(j, laet, hess_fn,
                                     base_soln, det_of_hess_at_base_soln, K_at_base_soln){
    grad_fn = grad_asymmetric(laet$X, laet$Y, laet$k, j)
    soln = tryCatch({pracma::broyden(grad_fn, c(1, 1), hess_fn(c(1, 1)))}, error=function(err) NaN)
    if (!is.list(soln)){
        if (j != floor(j)){
            return(c(NaN, NaN))
            }
        else
            return(asymmetric_saddlepoint_cdf(j + 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln)*0.5
                   + asymmetric_saddlepoint_cdf(j - 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln)*0.5)
    }

    soln = soln$zero
    sqrt_term_det = det(hess_fn(soln)) / det_of_hess_at_base_soln
    sqrt_term_w = (K_at_base_soln +
                    laet$k * (soln[1] - base_soln[1]) +
                    (j - 0.5) * soln[2] -
                    K_asymmetric(laet$X, laet$Y)(soln))

    if (min(sqrt_term_det, sqrt_term_w) < 1e-8){
        if (j != floor(j))
            return(c(NaN, NaN))
        else
            return(asymmetric_saddlepoint_cdf(j + 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln)*0.5
                   + asymmetric_saddlepoint_cdf(j - 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln)*0.5)
    }

    u = 2 * sinh(soln[2]/2.0) * sqrt(sqrt_term_det)
    w = sign(soln[2]) * sqrt(2 * sqrt_term_w)
    if (is.nan(u)){
        if (j != floor(j))
            return(c(NaN, NaN))
        else
            return(asymmetric_saddlepoint_cdf(j + 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln)*0.5
                   + asymmetric_saddlepoint_cdf(j - 1e-2, laet, hess_fn, base_soln, det_of_hess_at_base_soln, K_at_base_soln)*0.5)
    }

    return(c(pnorm(w) + dnorm(w)*(1/w - 1/u),
             pnorm(-w) - dnorm(w)*(1/w - 1/u)
             ))
}