get_moments_fast = function(laet) {
  # this, the less smart deconvolution function, simply computes
  # E(sum(X_i | i != j)) = E(sum(X_i)) - p^X_j = sum(x_probs) - p^X_j
  # and
  # V(sum(X_i | i != j)) = V(sum(X_i)) - p^X_j (1 - p^X_j)
  #    = sum(x_probs * (1 - x_probs)) - p^X_j (1 - p^X_j)
  # and then assumes sum(X_i | i != j) is normal.

  N = laet$N
  m = laet$m
  k = laet$k
  X = laet$X
  Y = laet$Y

  x_total_mean = sum(X)
  x_total_var = sum(X * (1 - X))
  y_total_mean = sum(Y)
  y_total_var = sum(Y * (1 - Y))
  
  denom = dnorm(m, x_total_mean, sqrt(x_total_var)) * dnorm(k, y_total_mean, sqrt(y_total_var))
  
  mean = sum(X * Y *
            dnorm(m - 1, mean = x_total_mean - X, sd = sqrt(x_total_var - X * (1 - X))) *
            dnorm(k - 1, mean = y_total_mean - Y, sd = sqrt(y_total_var - Y * (1 - Y)))
            )

  second_moment = 0

  for (i in 1:N) {
    second_moment_terms = (X * X[i] * Y * Y[i] *
                    dnorm(m - 2, x_total_mean - X[i] - X, sqrt(x_total_var - X[i] * (1 - X[i]) - X * (1 - X))) *
                    dnorm(k - 2, y_total_mean - Y[i] - Y, sqrt(y_total_var - Y[i] * (1 - Y[i]) - Y * (1 - Y))))
    # add all but the diagonal term to the running total
    second_moment_terms[i] = 0.0 # this covers the case where that term is NaN
    second_moment = second_moment + sum(second_moment_terms)
  }

  second_moment = second_moment + mean
  second_moment = second_moment / denom

  mean = mean / denom
  var = second_moment - mean**2
  cat("mean =", mean, "sd =", sqrt(var), "\n")
  list(mean=mean, sd=sqrt(var))
}

get_moments_slow = function(laet) {
  # this, the less smart deconvolution function, simply computes
  # E(sum(X_i | i != j)) = E(sum(X_i)) - p^X_j = sum(x_probs) - p^X_j
  # and
  # V(sum(X_i | i != j)) = V(sum(X_i)) - p^X_j (1 - p^X_j)
  #    = sum(x_probs * (1 - x_probs)) - p^X_j (1 - p^X_j)
  # and then assumes sum(X_i | i != j) is normal.

  N = laet$N
  m = laet$m
  k = laet$k
  X = laet$X
  Y = laet$Y

  # decide on lengths for FFT
  n_X = min(N, floor(4*m/3) + 3)
  n_Y = min(N, floor(4*k/3) + 5)

  if (n_X %% 2 == 0)
    n_X = n_X + 1
  if (n_Y %% 2 == 0)
    n_Y = n_Y + 1
 
  safe_fft = function(x, inverse=FALSE, len=NaN){
    if (!is.nan(len) && (length(x) < len))
            x = c(x, rep(0, len - length(x)))

    #if(length(x) == 2 && all(x == 0.5))
    #    1.0
    #else
    f = fft(x, inverse)
    if(inverse)
        return(f / length(f))
    else
        return(f)
  }

  # get exact unconditional distribution of sum(X_j) and sum(Y_j)
  dist_X = successive_paired_convolution(X, N+1)#m + 2)
  dist_Y = successive_paired_convolution(Y, N+1)#k + 2)

  denom = dist_X[m+1] * dist_Y[k+1] # equals P(sum(X_j) = m, sum(Y_j) = k)

  #fft_X = c(safe_fft(dist_X), rep(0, n_X - length(dist_X)))
  fft_X = safe_fft(dist_X[1:n_X], len=n_X)
  #fft_Y = c(safe_fft(dist_Y), rep(0, n_Y - length(dist_Y)))
  fft_Y = safe_fft(dist_Y[1:n_Y], len=n_Y)

  mean = sum(na.omit(X * Y *
            sapply(X, function(p) {Re(safe_fft(fft_X / safe_fft(c(1 - p, p), len=n_X), inverse=TRUE)[m])}) *
            sapply(Y, function(p) {Re(safe_fft(fft_Y / safe_fft(c(1 - p, p), len=n_Y), inverse=TRUE)[k])})
            ))

  second_moment = 0

  for (i in 1:N) {
    fft_X_i = safe_fft(c(1 - X[i], X[i]), len=n_X)
    fft_Y_i = safe_fft(c(1 - Y[i], Y[i]), len=n_Y)
    second_moment_terms = (X * X[i] * Y * Y[i] *
                    sapply(X, function(p) {Re(safe_fft(fft_X / fft_X_i / safe_fft(c(1 - p, p), len=n_X), inverse=TRUE)[m - 1])}) *
                    sapply(Y, function(p) {Re(safe_fft(fft_Y / fft_Y_i / safe_fft(c(1 - p, p), len=n_Y), inverse=TRUE)[k - 1])})
                    )

    # add all but the diagonal term to the running total
    second_moment_terms[i] = 0.0
    second_moment = second_moment + sum(na.omit(second_moment_terms))
  }

  second_moment = second_moment + mean
  second_moment = second_moment / denom
  mean = mean / denom
  var = second_moment - mean**2
  cat("mean =", mean, "sd =", sqrt(var), "\n")
  list(mean=mean, sd=sqrt(var))
}