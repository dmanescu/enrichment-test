library('sets')
#' Estimate the label probabilities of two pairs of binary labels as a function of length
#'
#' Given two binary (0 or 1) features, X and Y, of a set of individuals sequenced by some numerical attribute (length, say),
#' this function performs loess smoothing on those observed features in order to the estimate the probability
#' that each feature is 1, expressed as a function of the numerical attribute.
#'
#' The function requires four vectors, namely that of the lengths for which X = 0, X = 1, Y = 0 and Y = 1.
#' Not all four need to be provided, as the missing ones can be inferred provided enough information:
#' the function can accept a vector of all lengths, meaning labels with X = 0 (resp. Y = 0) can then be inferred
#' from the list of those with X = 1 (resp. Y = 1) and vice versa.
#'
#' If set to 'increasing' (resp. 'decreasing'), the monotonicity argument enforces that the resulting probabilities be monotone
#' by setting each label probability equal to the max (resp. min) of all estimated label probabilities
#' on a shorter length individual.
#' 
#' @param lengths_all NULL, or a vector of all lengths 
#' @param lengths_X_0 NULL, or a vector of lengths for which X is observed to be 0
#' @param lengths_X_1 NULL, or a vector of lengths for which X is observed to be 1
#' @param lengths_Y_0 NULL, or a vector of lengths for which Y is observed to be 0
#' @param lengths_Y_1 NULL, or a vector of lengths for which X is observed to be 1
#' @param monotonicity NULL, 'increasing' or 'decreasing'.
#' @export
#' @examples
#' > label_probs = loess_smoothed_marginals(1:10, lengths_X_0=c(4, 7, 8), lengths_Y_1=8:10, monotonicity='increasing')
#' > label_probs$X
#' > label_probs$Y
#' > label_probs$length
loess_smoothed_marginals = function(lengths_all = NULL,
                               lengths_X_0 = NULL, lengths_X_1 = NULL,
                               lengths_Y_0 = NULL, lengths_Y_1 = NULL,
                               monotonicity=NULL){
  stopifnot(sum(c(is.null(lengths_X_0), is.null(lengths_X_1), is.null(lengths_all))) <= 1)
  stopifnot(sum(c(is.null(lengths_Y_0), is.null(lengths_Y_1), is.null(lengths_all))) <= 1)

  X = loess_smooth(lengths_all, lengths_X_0, lengths_X_1, monotonicity)
  Y = loess_smooth(lengths_all, lengths_Y_0, lengths_Y_1, monotonicity)

  stopifnot(X$length == Y$length)

  list(X=X$loess, Y=Y$loess, length=X$length)
}

#' Estimate the probabilities of a binary label as a function of a numeric attribute
#'
#' Given a binary feature (0 or 1) of a set of individuals sequenced by some numerical attribute (length, say),
#' this function performs loess smoothing on the observed feature values in order to the estimate the probability
#' that each feature is 1, expressed as a function of the numerical attribute.
#'
#' If set to 'increasing' (resp. 'decreasing'), the monotonicity argument enforces that the resulting probabilities be monotone
#' by setting each label probability equal to the max (resp. min) of all estimated label probabilities
#' on a shorter length individual.
#' 
#' @param lengths_all NULL, or a vector of all lengths 
#' @param lengths_0 NULL, or a vector of lengths for which the attribute is observed to be 0
#' @param lengths_1 NULL, or a vector of lengths for which the attribute is observed to be 1
#' @param monotonicity NULL, 'increasing' or 'decreasing'.
#' @export
#' @examples
#' > probs = loess_smooth(1:10, lengths_0=c(6, 7, 10), monotonicity='increasing')
#' > label_probs$length
loess_smooth = function(lengths_all = NULL, lengths_0 = NULL, lengths_1 = NULL, monotonicity = NULL){
  stopifnot(sum(c(is.null(lengths_0), is.null(lengths_1), is.null(lengths_all))) <= 1)

  if (is.null(lengths_0))
    lengths_0 = setminus(lengths_all, lengths_1)
  else if (is.null(lengths_1))
    lengths_1 = setminus(lengths_all, lengths_0)

  lengths = c(lengths_0, lengths_1)
  labels = c(rep(0, length(lengths_0)), rep(1, length(lengths_1)))
  order_length = order(lengths)
  lengths_sorted = lengths[order_length]
  labels_sorted = labels[order_length]
  all = data.frame(lengths = lengths_sorted, labels = labels_sorted)
  all.loess = loess(labels ~ lengths, all)$fitted
  all.loess[all.loess > 1] = 1
  all.loess[all.loess < 0] = 0
  if(monotonicity == 'increasing'){
    all.loess = enforce_monotone(all.loess, max)
  } else if(monotonicity == 'decreasing'){
    all.loess = enforce_monotone(all.loess, min)
  }
  return (list(length = lengths_sorted, loess = all.loess, label = labels_sorted))
}

setminus = function(a, b){
  g = sets::gset_difference(sets::as.gset(a), sets::as.gset(b))
  return(rep(unlist(g), times=sets::gset_memberships(g)))
}

enforce_monotone = function(v, f=max){
    v2 = c(v)
    current_max = v2[1]
    i = 1
    for(x in v){
        current_max = f(current_max, x)
        v2[i] = current_max
        i = i + 1
    }
    return(v2)
}


"setminus = function(a, b){
  # the usual setdiff in R doesn't (to my knowledge) account for multiplicity,
  # this function does so by using run length encoding
  encoded_a = rle(a)
  encoded_b = rle(b)
  x = encoded_a$values %in% encoded_b$values
  y = encoded_b$values %in% encoded_a$values
  print(sum(x))
  print(sum(y))
  print(length(x))
  print(length(y))
  indices = encoded_a$values %in% encoded_b$values
  encoded_a$lengths[indices] = encoded_a$lengths[indices] - encoded_b$lengths[encoded_b$values %in% encoded_a$values]
  encoded_a$lengths[encoded_a$lengths < 0] = 0
  inverse.rle(encoded_a)
}"