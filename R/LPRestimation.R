#' Estimate LPR
#'
#' This function estimates the LPR on a data set.
#' The user has the option of providing two data sets (one to use as training
#' and the other as a test set). Otherwise, if only one data set is provided,
#' the function will estimate LPR values for the provided data set using
#' leave one out CV.
#'
#' Users have the option of using bagging to produce more robust LPR estimates,
#' and can choose whether to fit the precision function using a local
#' polynomial kernel smoother or a spline.
#'
#' @param scores Numeric vector containing classifier scores for each observation
#' @param labels A 0-1 vector containing the true status of each observation
#' @param newdata An optional numeric vector containing scores for a new data
#'                set, e.g. a test set. If omitted, then LPRs for \code{scores}
#'                are estimated by leave-one-out cross-validation.
#' @param num_bags Integer indicating the number of bagging steps. If 0,
#'                 no bagging is done.
#' @param method Method used to estimate the precision function
#' @param ... Additional parameters to pass to the function for estimating the
#'            precision (e.g., additional arguments to fitKernelLPR).
#'
#' @export
#'
estimateLPR <- function(scores, labels,
                        newdata = NULL, num_bags,
                        method = c("kernel", "spline"), ...){
  
  # set estimation method for precision function (kernel or spline)
  method = match.arg(method, choices = c("kernel", "spline"))
  
  if(!is.null(newdata)){
    # train on scores, labels and predict on newdata
    output <- predictNewLPR(train_scores = scores, train_labels = labels,
                            test_scores = newdata,
                            num_bags = num_bags,
                            method = method, ...)
    
  } else{
    # use leave-one-out cross validation to estimate LPR of each score
    output <- sapply(scores, function(x){
      cat(sprintf("...bagging obs. %s", match(x, scores)), "\n")
      
      predictNewLPR(train_scores = scores[scores != x],
                    train_labels = labels[scores != x],
                    test_scores = x,
                    num_bags = num_bags,
                    method = method, ...)})
    
    #     # foreach parallelized version
    #     output <- foreach(x = scores, .combine = c) %dopar% {
    #       predictNewLPR(train_scores = scores[scores != x],
    #                     train_labels = labels[scores != x],
    #                     test_scores = x,
    #                     num_bags = num_bags,
    #                     method = method, ...)
    #     }
  }
  
  cat("\n all done! \n")
  
  return(output)
}


#' Estimate LPR on a new data set using bagging
#'
#' This function estimates the LPR given training data and a new test set.
#' Users have the option of using bagging to produce more robust LPR estimates,
#' and can choose whether to fit the precision function using a local
#' polynomial kernel smoother or a spline.
#'
#' @param train_scores Numeric vector containing classifier scores for each observation
#' @param train_labels A 0-1 vector containing the true status of each observation
#' @param test_scores Numeric vector containing the scores for the test set.
#' @param num_bags Integer indicating the number of bagging steps. If 0,
#'                 no bagging is done.
#' @param method Method used to estimate the precision function
#' @param ... Additional parameters to pass to the function for estimating the
#'            precision function (e.g., fitKernelLPR)
#' @param verbose If TRUE, prints out at most five progress messages
#'                during LPR estimation using bagging,
#'                since this can be slow when the
#'                number of bags or observations is large.
#'
predictNewLPR <- function(train_scores, train_labels,
                          test_scores,
                          num_bags,
                          method = c("kernel", "spline"),
                          verbose = FALSE,
                          ...){
  
  # set estimation method for precision function (kernel or spline)
  method = match.arg(method, choices = c("kernel", "spline"))
  FUN <- fitSplineLPR
  if(method == "kernel") FUN = fitKernelLPR
  
  if(num_bags > 0){
    # estimate test sample LPRs using bagged training samples
    bag_idx <- replicate(num_bags, createBag(train_labels))
    bag_lprs <- apply(bag_idx, 2, function(x){
      fxn_args <- list(train_scores = train_scores[x],
                       train_labels = train_labels[x],
                       test_scores = test_scores, ...)
      do.call(FUN, fxn_args)[["LPR_test"]]})
    
    #     foreach parallelized version
    #     bag_lprs <- foreach(i = seq_len(num_bags), .combine = cbind, .verbose = verbose) %dopar% {
    #       bag_idx <- sample.int(n_train, replace = TRUE)
    #
    #       fxn_args <- list(train_scores = train_scores[bag_idx],
    #                        train_labels = train_labels[bag_idx],
    #                        test_scores = test_scores, ...)
    #       do.call(FUN, fxn_args)[["LPR_test"]]
    #     }
    
    # take the median of the bagged LPR estimates
    if(is.matrix(bag_lprs)){
      lpr_test <- apply(bag_lprs, 1, median, na.rm=TRUE)
    } else{
      lpr_test <- median(bag_lprs, na.rm=TRUE)
    }
    
  } else{
    # fit the precision function without bagging
    fxn_args <- list(train_scores = train_scores,
                     train_labels = train_labels,
                     test_scores = test_scores, ...)
    lpr_test <- do.call(FUN, fxn_args)[["LPR_test"]]
  }
  
  return(lpr_test)
}

createBag <- function(labels){
  n <- length(labels)
  idx <- sample.int(n, replace=TRUE)
  # if any group is not represented, resample
  if( any(table(labels[idx]) == 0) ){
    createBag(labels)
  } else{
    return(idx)
  }
}


#' Estimate the precision function
#'
#' Given a vector of classifier scores and true labels,
#' \code{estimatePrecision} estimates the precision as a function of the
#' lower tail probabilities (cdf) of the scores.
#'
#' @param scores A vector of classifier scores.
#' @param labels True labels for each of the instances in \code{scores}.
#' @param plot If TRUE, generates a plot of the precision function.
#'
#' @return A list with three elements.
#' \describe{
#'  \item{u}{Empirical CDF evaluated at each score.}
#'  \item{prec}{Empirical precision estimates for each score.}
#'  \item{size}{Sample size used to estimate each precision value.}
#'  }
#'
estimatePrecision <- function(scores, labels, plot=FALSE){
  # This function empirically estimates the precision function
  # Here, the precision is a function of u, the empirical percentiles
  # (actually proportions rather than percents here, but same concept)

  n <- length(scores)
  # find empirical percentiles
  u <- convertToPercentile(scores)
  # find number of values used to estimate precision for each percentile
  size <- sapply(scores, function(i) sum( scores > i ) );

  # compute precision for each percentile
  prec = sapply(scores, function(i) sum( (scores > i) & (labels > 0) ))/size

  # by this method, the last term is undefined since no scores
  # can be larger than the maximum.
  # to address this, inherit the last precision value for convenience
  na_ind <- is.na(prec)
  prec[na_ind] <- prec[!na_ind][which.max(u[!na_ind])]

  if(plot){
    par(mfrow=c(2,1));
    plot(u, prec,
         xlab="score percentile: P(S < s)",
         ylab="precision: P(Q = 1 | S > s)");
    plot(u, (1-u)*prec,
         xlab="score percentile: P(S < s)",
         ylab="joint probability: P(Q = 1, S > s)")
  }

  # return the score percentiles, precision, and
  # number of values used to estimate the precision for each percentile
  return(list(u = u, prec = prec, size = size))
}


#' Transform classifier scores to percentiles
#'
#' Given a vector of classifier scores \code{convertToPercentile}
#' returns the empirical CDF evaluated at each of the scores.
#'
#' @param scores A vector of classifier scores.
#'
#' @return A vector of percentiles (as decimals rather than percentages), the
#' same length as \code{scores}.
#'


convertToPercentile <- function(scores){
  sapply(scores, function(i) mean(scores <= i))
}


#' Create cross-validation folds (from caret package)
#'
#' This is the caret package's "createFolds" function, copied here
#' because this is the only function used from the caret package,
#' so it seemed excessive to list caret as a dependency.
#'
#' @param y outcome vector, i.e. the indicator labels for each instance
#' @param k number of folds to create
#' @param list logical; should folds be returned as a list or matrix?
#'
caretFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE){
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

#' Estimate LPR using local polynomial smoothing
#'
#' Given a vector of classifier scores and true labels,
#' \code{fitSplineLPR} fits a local quadratic smoother to
#' the precision function on training data, and provides predictions
#' on a test set of classifier scores.
#'
#' @param train_scores A vector of classifier scores to use as the training set
#' @param test_scores A vector of classifier scores to obtain LPRs for. If
#'                    this is set to NULL, then the function only returns the
#'                    best smoothing parameter chosen from an nfold CV on
#'                    the training data.
#' @param train_labels True labels for the instances in the training set
#' @param param Smoothing parameter for the spline fit; if NULL, estimate
#'              using \code{nfold} cross validation. If a non-NULL is
#'              provided and test_scores = NULL, then this function
#'              trivially returns the provided value of param.
#' @param nfold Use \code{nfold} cross validation to estimate param,
#'              if not given. By default, 5-fold CV is done.
#'
#' @return A list with three elements.
#' \describe{
#'  \item{u}{Empirical CDF evaluated at each score.}
#'  \item{prec}{Empirical precision estimates for each score.}
#'  \item{size}{Sample size used to estimate each precision value.}
#'  }
#'

fitKernelLPR <- function(train_scores, train_labels,
                         test_scores = NULL, param = NULL,
                         weights = NULL, nfold = 5){
  
  # estimate the precision function on training data
  prec_list = estimatePrecision(train_scores, train_labels)
  n <- length(prec_list[[1]])
  
  # use weights proportional to sample size used in estimating cdf, u
  if(is.null(weights)) weights = 1/(n - prec_list[[3]])
  # another popular choice is equal weight per observation,
  # weights <- rep(1, n)
  
  # use cross-validation to find the best spline smoothing parameters
  if(is.null(param)){
    folds <- caretFolds(y = as.factor(train_labels), k = nfold)
    param <- optimize(kernelloss,
                      interval=c(0,10),
                      x = prec_list[[1]],
                      y = prec_list[[2]],
                      folds = folds,
                      w = weights)$min
  }
  
  if(!is.null(test_scores)){
    # estimate cdf (i.e., u, our x values) for the test cases
    u_test <- sapply(test_scores, function(x) mean(train_scores < x))
    
    # fit the kernel smoother on the entire training data set
    # and predict on the test set
    kernel_fit = locpol::locPolSmootherC(x = prec_list[[1]],
                                 y = prec_list[[2]],
                                 xeval = u_test,
                                 bw = param,
                                 deg = 2,
                                 kernel = EpaK,
                                 weig = weights)
    
    # estimate LPR per equation 2.4 in Jiang et al (2013):
    # LPR = G_k(u) - (1-u) G'k(u)
    LPR_test <- kernel_fit$beta0 - (1- u_test)*kernel_fit$beta1
    # censor the estimated LPR to be between 0 and 1
    LPR_test <- ifelse(LPR_test < 0, 0,
                       ifelse(LPR_test > 1, 1, LPR_test))
    
    return(list(LPR_test = LPR_test, param = param))
  }
  
  return(param)
}

#' Compute mean squared error loss, calculated per fold (for cv)
#' @param param bandwidth parameter for kernel smoother
#' @param x a vector giving the values of the predictor variable
#' @param y responses
#' @param folds a list of vectors, each giving the indices for each CV fold
kernelloss <- function(param, x, y, folds, w){
  error2 <- mapply(function(i, x, y, param, w){
    z <- locPolSmootherC(x = x[-i], y = y[-i],
                         xeval = x[i],
                         bw = param, deg = 2,
                         kernel = EpaK,
                         weig = w[-i])
    pred_y <- z$beta0 - (1 - x[i])*z$beta1
    mean((y[i] - pred_y)^2)},
    i = folds,
    MoreArgs = list(x = x, y = y, param = param, w = w))
  return(sum(error2))
}

#' Estimate LPR using a spline fit
#'
#' Given a vector of classifier scores and true labels,
#' \code{fitSplineLPR} fits a spline to the precision
#' function on training data, and provides predictions
#' on a test set of classifier scores.
#'
#' @param train_scores A vector of classifier scores to use as the training set
#' @param test_scores A vector of classifier scores to obtain LPRs for. If
#'                    this is set to NULL, then the function only returns the
#'                    best smoothing parameter chosen from an nfold CV on
#'                    the training data.
#' @param train_labels True labels for the instances in the training set
#' @param param Smoothing parameter for the spline fit; if NULL, estimate
#'              using \code{nfold} cross validation. If a non-NULL is
#'              provided and test_scores = NULL, then this function
#'              trivially returns the provided value of param.
#' @param nfold Use \code{nfold} cross validation to estimate param,
#'              if not given. By default, 5-fold CV is done.
#'
#' @return A list with three elements.
#' \describe{
#'  \item{u}{Empirical CDF evaluated at each score.}
#'  \item{prec}{Empirical precision estimates for each score.}
#'  \item{size}{Sample size used to estimate each precision value.}
#'  }
#'
fitSplineLPR <- function(train_scores, train_labels,
                         test_scores = NULL, param = NULL,
                         nfold = 5){
  # estimate the precision function
  prec_list = estimatePrecision(train_scores, train_labels)
  n <- length(prec_list[[1]])
  weights = 1/(n - prec_list[[3]]) # weights for spline fit
  
  # use cross-validation to find the best spline smoothing parameters
  if(is.null(param)){
    folds <- caretFolds(y = as.factor(train_labels), k = nfold)
    param <- optimize(splineloss,
                      interval=c(0,1),
                      x = prec_list[[1]],
                      y = prec_list[[2]],
                      folds = folds,
                      w = weights)$min
  }
  
  if(!is.null(test_scores)){
    # fit the spline model on the entire training data set
    spline_fit = smooth.spline(x = prec_list[[1]],
                               y = prec_list[[2]],
                               spar = param,
                               w = weights,
                               control.spar = list(tol=1e-06))
    
    # estimate precision for the test cases
    u_test <- sapply(test_scores, function(x) mean(train_scores < x))
    prec = predict(object = spline_fit, x = u_test)$y
    dprec = predict(object = spline_fit, x = u_test, deriv = 1)$y
    LPR_test = prec - (1 - u_test)*dprec
    
    # censor the estimated LPR to be between 0 and 1
    LPR_test <- ifelse(LPR_test < 0, 0,
                       ifelse(LPR_test > 1, 1, LPR_test))
    return(list(LPR_test = LPR_test, param = param))
  }
  
  return(param)
}

#' Define mean squared error loss function, calculated per fold (for cv)
#' @param param smoothing parameter, typically in (0, 1]
#' @param x a vector giving the values of the predictor variable
#' @param y responses
#' @param folds a list of vectors, each giving the indices for each CV fold
#'
splineloss=function(param, x, y, folds, w=NULL){
  if(is.null(w)) w <- rep(1, length(x))
  error2 = mapply(function(i, y, param, x, w){
    z = smooth.spline(x[-i], y[-i], w = w[-i],
                      spar = param, cv = NA,
                      control.spar = list(tol=1e-06))
    fit = predict(object = z, x[i])
    mean((y[i] - fit$y)^2)},
    i = folds,
    MoreArgs = list(y = y, param = param, x = x, w))
  return(sum(error2));
}

