% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mLPRestimation.R
\name{estimatemLPR}
\alias{estimatemLPR}
\title{Estimate LPR}
\usage{
estimatemLPR(
  S_train,
  S_test,
  Y_train,
  graph_structure,
  X_train = NULL,
  X_test = NULL,
  cdensity_method = "empirical",
  pr_s_y_module = "LPR-",
  pr_y_module = "hierarchical-",
  approx_mlpr = "full",
  num_cores = 1,
  saving_name = "",
  ...
)
}
\arguments{
\item{S_train}{Numeric matrix containing classifier scores for each training observation,
of shape=n_train x K, where n_train is the number of training objects and
K is the number of class.}

\item{S_test}{Numeric matrix containing classifier scores for each testing observation,
of shape=n_test x K, where n_test is the number of testing objects and
K is the number of class.}

\item{Y_train}{Binary matrix containing labels for each training observation, of the shape
as \code{S_train}.}

\item{graph_structure}{A list that indicates the parents of each node.}

\item{X_train}{Numeric matrix containing the features for training data.}

\item{X_test}{Numeric matrix containing the features for testing data.}

\item{cdensity_method}{Methods of estimating the conditional probablity of \code{Pr(Y_i|Y_pa(i))},
including 'empirical', 'logistic' (requires \code{X_train} and \cdoe{X_test}),
'svm' (requires \code{X_train} and \cdoe{X_test}).}

\item{pr_s_y_module}{A string indicating which method and model to estimate \code{Pr(S|Y) or Pr(Y|S)}. The method
and the model is concatenated by '-'. Methods include 'LPR', 'Gaussian', 'thresholding',
models include 'logistic', 'svm'.}

\item{pr_y_module}{A string indicating which method and model to estimate \code{Pr(Y)}. The method
and the model is concatenated by '-'. Methods include 'hierarchical', 'empirical',
models include 'logistic', 'svm'.}

\item{approx_mlpr}{Three versions of mLPR estimations: 'full', 'nbh', 'indpt'.}

\item{num_cores}{Number of cores for parallel computing.}

\item{saving_name}{The prefix name used to save intermediate objects.}
}
\value{
A numetric matrix containing mLPR values for testing objects, of shape=n_test x K.
}
\description{
This function estimates the mLPR on a data set.
The user needs to provide two data sets (one to use as training
and the other as a test set).
}
