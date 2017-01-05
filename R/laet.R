#TODO: asymmetric saddlepoint uses approximated hessian function
library('pracma')
library("zoo")
asymmetric_methods = c("exact", "saddlepoint")

#' A Length-Aware Enrichment Test
#'
#' This function implements the symmetric length-aware enrichment test, a generalisation
#' of Fisher's exact test. Fisher's exact test invokes the hypergeometric distribution
#' to describe the size of the overlap between N objects of which m have property X with equal likelihood,
#' and k of the N objects have property Y with equal likelihood. In this symmetric, length-aware generalisation
#' we permit the N objects to have independent but not identical probabilities of having properties X & Y.
#' This function also implements an asymmetric generalisation of Fisher's exact test which may be more appropriate
#' under some circumstances. In it we assume it is known which m of the N objects have property X, but each
#' object has a certain probability of having property Y.
#'
#' The binomial method (method="binom") calculates the distribution of the overlap given only
#' the marginal probabilities x_probs and y_probs. In this situation the probabilities that each
#' of the N objects has both properties X and Y are given by x_probs * y_probs, resulting in a simple binomial distribution.
#' The exact method uses dynamic programming to calculate the distribution of the overlap exactly, conditioned on the
#' totals m and k, and is fairly computationally intensive.
#' The normal method calculates the first and second moments of the same distribution to generate a normal distribution.
#' The moments are derived using FFTs and provide neither a significant speedup nor particularly good accuracy.
#' The faster normal method (method="fast_normal") approximates the moments and provides a more meaningful speedup
#' in exchange for only a small reduction in accuracy over the normal approximation which derives the moments exactly.
#' A Monte Carlo method is available and provides a good balance between running time and accuracy, but cannot be used
#' if full precision is required or the p-values in question are too small.
#' Finally a saddlepoint approximation is available. Its implementation is considerably more involved than the normal
#' approximation but gives much greater accuracy even down to the lowest p-values and its running time is excellent.
#'
#' For asymmetric tests (test="x-cond" and test="y-cond") we assume those objects with property X (resp. Y) are known,
#' and the probabilities that each of the objects has property Y are given by y_probs (resp. X, x_probs).
#' As such it is assumed that the vector of x_probs (resp. y_probs) consists only of the values 1 and 0
#' corresponding to objects which have/do not have property X (resp. Y), and therefore this vector sums to m (resp. k) by design.
#' The exact method for the asymmetric tests is more performant than for the symmetric test, but a saddlepoint
#' approximation is also provided.
#'
#' The results of the test come in a variety of forms. Given kind="pmf", the test will return a vector giving
#' the probability that each observed value occurs, i.e. a probability density function. With kind="cdf" the
#' output will be the probability that a result "at least as surprising" as each observed value occurs.
#' The definition of "at least as surprising" differs based on the side given. For side="lt" and side="gt"
#' we "at least as small" or "at least as large" are selected. With side="two-sided" we define the p-value to
#' be the probability of an event will occur whose probability mass is less than or equal that of the observed value.
#'
#' @param observed observed value or vector of observed values
#' @param m number of the N objects which have property X
#' @param k number of the N objects which have property Y
#' @param x_probs vector of length N, holding the probabilities that each object has property X
#' @param y_probs vector of length N, holding the probabilities that each object has property Y
#' @param test one of "symmetric", "x-cond" (in which case valid x_probs are 1 and 0) and "y-cond"
#' @param side one of "gt", "lt" and "two-sided" (in which we define the p-value as the probability of at least as unlikely an outcome)
#' @param kind one of "cdf" and "pmf" (parameter side is optional if "pmf" is given)
#' @param method method(s)/approximation(s) to use. Valid methods are MC, fast_normal, normal, binom, exact and saddlepoint. Only the latter two are available for asymmetric tests.
#' @param MC.iterations (optional) number of iterations to use if Monte Carlo simulation is selected
#' @export
#' @examples
#' > laet_out = laet(0:3, m=3, k=5, x_probs=rep(0.2, 8), y_probs=rep(0.3, 8), test="symmetric",
#'        side="lt", kind="cdf", method=c("exact", "saddlepoint"))
#' > print(laet_out$results$saddlepoint)
#' > print(laet_out$results$exact)
laet = function(observed=NULL, m=NULL, k=NULL, x_probs=NULL, y_probs=NULL,
                test = c("symmetric", "x-cond", "y-cond"),
                side = c("two-sided", "gt", "lt"),
                kind = c("cdf", "pmf"),
                method = c("exact", "MC", "fast_normal", "normal", "binom", "saddlepoint"),
                MC.iterations = NULL) {
  if(any(c(kind) == "cdf") && length(c(side)) != 1){
    stop("if requesting cdf, the side must be specified")
  }

  stopifnot(is.element(side, c('two-sided', 'gt', 'lt')))
  if((min(observed) < max(0, m + k - length(x_probs))) || max(observed) > min(m, k)){
    stop("some observed values are not possible (must be between max(0, m + k - N) and min(m, k))")
  }

  laet_obj = base_laet_config(x_probs, y_probs, m, k, observed)

  if(laet_obj$k > laet_obj$N / 2){
    laet_obj$k = laet_obj$N - laet_obj$k
    laet_obj$Y = 1 - laet_obj$Y
    laet_obj$o = laet_obj$m - laet_obj$o
    side = swap_side(side)
  }

  if(laet_obj$m > laet_obj$N / 2){
    laet_obj$m = laet_obj$N - laet_obj$m
    laet_obj$X = 1 - laet_obj$X
    laet_obj$o = laet_obj$k - laet_obj$o
    side = swap_side(side)
  }

  o = laet_obj$o

  # if necessary prevent any calculations from taking calculation shortcuts based on the observed values,
  #  as two-sided calculation requires full cdfs to calculate any single p-value
  if(any(c(kind) == "cdf") && any(c(side) == "two-sided"))
    laet_obj$o = 0:laet_obj$N

  if(length(c(test)) != 1)
    stop("test must be exactly one of symmetric, x-cond, y-cond")

  if(test == "symmetric") {
    binom = dist(laet_obj, pmf_function = get_binom_dist)
    fast_normal = normal_dist(laet_obj, get_moments_fast)
    slow_normal = normal_dist(laet_obj, get_moments_slow)
    MC = dist(laet_obj, pmf_function = function(laet) get_MC_dist(laet, MC.iterations))
    saddlepoint = dist(laet_obj, cdf_function=get_saddlepoint_dist)
    exact = dist(laet_obj, pmf_function=get_exact_dist)
  
    tests = list(binom=binom, fast_normal=fast_normal, normal=slow_normal, MC=MC, saddlepoint=saddlepoint, exact=exact)
  } else {
    if(length(setdiff(method, asymmetric_methods)) > 0)
      stop(paste("valid methods for asymmetric tests are", toString(asymmetric_methods)))

    if(test == "y-cond")
      # interchange X<->Y and m<->k at this level so the underlying implementation can assume x-cond
      laet_obj_to_pass_in = list(N=laet_obj$N, m=laet_obj$k, k=laet_obj$m, X=laet_obj$Y, Y=laet_obj$X, o=laet_obj$o)
    else
      laet_obj_to_pass_in = laet_obj

    if(length(setdiff(laet_obj_to_pass_in$X, c(0, 1))) > 0)
      stop("for asymmetric tests, one set of label probabilities must be identically one or zero")
    stopifnot(sum(laet_obj_to_pass_in$X) == laet_obj_to_pass_in$m)

    exact = dist(laet_obj_to_pass_in, pmf_function=get_exact_asymmetric_dist)
    saddlepoint = dist(laet_obj_to_pass_in, cdf_function=get_saddlepoint_asymmetric_dist)

    tests = list(exact=exact, saddlepoint=saddlepoint)
  }

  results = list()
  summary = paste("Probability of overlap",
                  toString(observed, 20),
                  "from",
                  toString(laet_obj$N),
                  "objects of which",
                  toString(laet_obj$m),
                  "have property X, and",
                  toString(laet_obj$k),
                  "have property Y",
                  sep=" ")
  for(m in c(method)){
    dist_for_test = tests[[m]]
    if(length(kind) == 2) {
      results[[m]] = list(pmf=dist_for_test$pmf(o, side),
                             cdf=dist_for_test$cdf(o, side))
    } else {
      results[[m]] = dist_for_test[[kind]](o, side)
    }

    summary = paste(summary,
                    paste("\t",
                          m,
                          ": ",
                          toString(signif(results[[m]], digits=2), width=65),
                          sep=""),
                    sep='\n')
  }
  
  list(summary=summary, results=results)
}