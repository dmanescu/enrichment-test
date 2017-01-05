swap_side = function(side){
  side = c(side)
  lt = side == "lt"
  gt = side == "gt"
  side[lt] = "gt"
  side[gt] = "lt"
  return(side)
}

successive_paired_convolution = function(L, max_len=0) {
  # given a list L = (p_1, p_2, ...) of probabilities, return the convolution of all the pairs (1 - p_i, p_i) together.
  # if p_i are probabilities the result is the pmf of the sum of length(L) bernoulli RVs each with probability p_i.
  # for optimisation purposes you can specify you only need the pmf up to an observed value max_len.
  # note that as vectors are 1-indexed, results[1] is P(sum_j X_j) = 0, not P(sum_j X_j) = 1, and so on.

  results = c(1)
    
  for (p in L) {
    results = capped_convolution(results, c(1 - p, p), max_len=max_len)
  }

  results
}

iterative_paired_convolution = function(L, max_len=0) {
  # same as successive_paired_convolution but stores the intermediate results in the rows of a matrix.
  # the i-th row in the matrix (1 <= i <= length(L) + 1) stores the distribution of
  # sum_{j >= i - 1} X_j
  # where X_j are independent Bernoulli(p_j) RVs (p_j being the j-th element of L). Hence
  # results[i, k] = P(sum_{j >= i - 1} X_j = k - 1)
  # in particular the first row results[0,] matches the output of successive_paired_convolution

  if (max_len == 0)
    max_len = length(L) + 1

  i = length(L) + 1
  intermediate = c(1)
  results = matrix(0, nrow=length(L) + 1, ncol=max_len)
  results[i, 1] = 1

  for (p in rev(L)) {
    i = i - 1
    results[i,] = capped_convolution(results[i+1,], c(1 - p, p), max_len=max_len)
  }

  results
}

capped_convolution = function(X, Y, max_len=0){
    # simple implementation of a convolution with a cap on the length of the output (for efficiency purposes)

    if (max_len == 0)
        max_len = length(X) + length(Y) - 1

    results = rep(0, max_len)
    i = 1

    for (v_x in X){
        if (i > max_len)
            break

        slice_Y = Y[1:min(length(Y), max_len - i + 1)]
        results_end = (i + length(slice_Y) - 1)
        results[i:results_end] = results[i:results_end] + v_x * slice_Y

        i = i + 1
    }

    results
}

# return vector of sampled Bernoulli(P[i]) RVs
sample_from_probs_uncond = function(P) runif(length(P)) < prob

# return vector of sampled X_i ~ Bernoulli(P[j]) RVs conditioned on sum(X_j) = total
# more efficient than sampling with rejection.
# expects as a third argument the precomputed distributions of sum_{j >= i}(X_j) for each i
# (which can be found using the iterative_paired_convolution function)
sample_from_probs_cond_on_total = function(probs, total, uncond_dists){
    j = 0
    result = rep(0, length(probs))
    r = runif(length(probs))

    for(prob in probs) {
        j = j + 1

        if (total == 0) {
            # we have already seen all our allowed 1's
            break
        }

        # computes P(X_j = 1 | X_1, X_2, ..., X_{j - 1}, sum_{i > j} X_j)
        p = prob *  uncond_dists[j + 1, total] / uncond_dists[j, total + 1]
        sampled_j = r[j] < p
        result[j] = sampled_j
        total = total - sampled_j
    }

    result
}

library(zoo) # used for vector fill-forward implementation

# implementation of a simple wrapper that abstracts away the difficulty of thinking about
# PMF <-> CDF conversions and CDF alternatives (handles two-sided, gt, lt alternatives)
# each of the length-aware enrichment tests implements a distribution function which is packaged
# into a dist object. The resulting dist object is just a list,
#    d = list(pmf=function(x){...},
#             cdf=function(q){...})
# which allows the length-aware enrichment test code to query
# the underlying implementation for either a pmf or a cdf.
# cdf_functions are expected to return 2xN matrices
#    - the first row gives the values P(Z <= z) for each z
#    - the second row gives the values 1 - P(Z <= z)
# this is useful for precision reasons
dist = function(laet, pmf_function = NULL, cdf_function = NULL){
    new_dist = list(laet = laet)

    if (!is.null(pmf_function)){
        # use the pmf function provided
        new_dist$pmf = function(x, side){
            pmf_function(laet)[x + 1]
        }
    }
    else {
        # convert the given cdf_function to generate a pmf function
        new_dist$pmf = function(x, side){
            result = clean_cdf(cdf_function(laet))
            cdf_lt = result[1,]
            cdf_gt = result[2,]
            pmf = diff(c(0, cdf_lt))
            mode = which(pmf == max(pmf))[1]
            better_pmfs_for_upper_tail = -diff(c(cdf_gt[mode:length(cdf_gt)], 0))
            pmf[mode:length(pmf)] = better_pmfs_for_upper_tail
            pmf[x + 1]
        }
    }

    if (!is.null(cdf_function)){
        # use the cdf function provided
        new_dist$cdf = function(q, side){
            result = clean_cdf(cdf_function(laet))
            if(side == "lt"){
                return(result[1,q+1])
            } else if(side == "gt"){
                return(result[2,q+1])
            } else {
                # take CDFs with lt and gt alternatives respectively.
                # For each observed value x on either tail,
                # find the observed value x' on the other tail which maximises 0 <= P(X = x') <= P(X = x),
                # that is, the most likely outcome x' which is "at least as surprising" as x but lies on the other tail.
                cdf_lt = result[1,]
                cdf_gt = result[2,]
                pmf = diff(c(0, cdf_lt))
                mode = which(pmf == max(pmf))[1]
                better_pmfs_for_upper_tail = -diff(c(cdf_gt[mode:length(cdf_gt)], 0))
                pmf[mode:length(pmf)] = better_pmfs_for_upper_tail
                two_sided_cdf_pvalue(pmf, cdf_lt, cdf_gt, q)
            }
        }
    }
    else {
        # convert the given pmf_function which generates the whole pmf to a cdf function
        new_dist$cdf = function(q, side){
            pmf = pmf_function(laet)
            if (side == "lt")
              cumsum(pmf)[q + 1]
            else if (side == "gt")
              cumsum(rev(pmf))[length(pmf) - q]
            else {
              cdf_lt = cumsum(pmf)
              cdf_gt = rev(cumsum(rev(pmf)))
              two_sided_cdf_pvalue(pmf, cdf_lt, cdf_gt, q)
            }
        }
    }

    new_dist
}

two_sided_cdf_pvalue = function(pmf, cdf_lt, cdf_gt, q){
    # simple function designed to handle the process of extracting two-sided
    # p-values. For our purposes the p-value for outcome x is defined as the sum of probabilities
    # of all outcomes which are "at least as surprising" as x. We also assume that all pmfs
    # we will handle are unimodal, i.e. pmfs are increasing up to the mode and then decreasing
    # from that point onwards.
    mode = which(pmf == max(pmf))[1]
    pmf_order = order(pmf)

    other_tail_prob = rep(0, length(cdf_lt))
    tail_prob = rep(0, length(cdf_lt))

    # search upper tail for values to match lower tail
    pmf_order_with_nans_below_mode = as.vector(pmf_order)
    pmf_order_with_nans_below_mode[pmf_order <= mode] = NaN
    matched_indexes_on_upper_tail = zoo::na.locf(pmf_order_with_nans_below_mode, na.rm=F)[pmf_order <= mode]
    other_tail_prob[1:mode] = zoo::na.fill(cdf_gt[matched_indexes_on_upper_tail], 0)
    tail_prob[1:mode] = cdf_lt[1:mode]

    # search lower tail for values to match upper tail
    pmf_order_with_nans_above_mode = as.vector(pmf_order)
    pmf_order_with_nans_above_mode[pmf_order >= mode] = NaN
    matched_indexes_on_lower_tail = zoo::na.locf(pmf_order_with_nans_above_mode, na.rm=F)[pmf_order >= mode]
    other_tail_prob[length(tail_prob):mode] = zoo::na.fill(cdf_lt[matched_indexes_on_lower_tail], 0)
    tail_prob[mode:length(tail_prob)] = cdf_gt[mode:length(tail_prob)]
    return((tail_prob + other_tail_prob)[q + 1])
}

clean_cdf = function(result){
    # clean up NaNs, etc.
    result[1,] = zoo::na.locf(c(0, result[1,], 1), na.rm=F, fromLast=F)[2:(length(result[1,]) + 1)]
    result[2,] = zoo::na.locf(c(1, result[2,], 0), na.rm=F, fromLast=T)[2:(length(result[1,]) + 1)]
    return(result)
}

normal_dist = function(laet, moments_function){
    # convert a moments_function = function(laet) list(mean=x, sd=y) into a regular dist.
    # this is used for the two normal approximations
    dist(laet,
        pmf_function = function(laet) {
            x = 0:min(laet$m, laet$k)
            moments = moments_function(laet)
            pnorm(x + 0.5, moments$mean, moments$sd) - pnorm(x - 0.5, moments$mean, moments$sd)
        },
        cdf_function = function(laet) {
            moments = moments_function(laet)
            q = 0:min(laet$m, laet$k)
            t(matrix(c(
               pnorm(q + 0.5, moments$mean, moments$sd) - pnorm(-0.5, moments$mean, moments$sd),
               pnorm(q - 0.5, moments$mean, moments$sd, lower.tail = FALSE) - pnorm(min(laet$m, laet$k) + 0.5, moments$mean, moments$sd, lower.tail = FALSE)
               ), ncol=2))
        })
}

base_laet_config = function(X, Y, m, k, o) {
  # distills the basic inputs of the length aware enrichment test into a simple list to be passed around internally
  X[X < 0] = 0
  X[X > 1] = 1
  
  Y[Y < 0] = 0
  Y[Y > 1] = 1
  
  list(X=X, Y=Y, N=length(X), m=m, k=k, o=o)
}