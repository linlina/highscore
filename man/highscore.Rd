\name{highscore}
\alias{highscore}
\title{Regularized score matching method}
\usage{
highscore(x, model, lambda = 0, centered = TRUE, tol = 1e-6, maxit = 1000, ...)
}
\arguments{
  \item{x}{The n x p data matrix, that is, n observations and p variables}
  \item{model}{Either "gaussian", "nonnegative" (for non-negative Gaussians), "conditional" (for conditional Gaussian)}

  \item{lambda}{The l1 regularization parameter}

  \item{centered}{Either TRUE or FALSE: if TRUE then (1) for Gaussian and non-negative Gaussian models, we estimate the mean/shift parameter mu and (2) for conditional Gaussian models, we estimate the matrix B for x_ix_j cross terms}

  \item{tol}{Convergence threshold}

  \item{maxit}{Maximum number of iterations to run}

  \item{...}{ignored}
}
\description{
Implements the regularized score matching method from Lin, Drton and Shojaie (2015) http://arxiv.org/abs/1507.00433
}
\examples{
library(mvtnorm)
n <- 500
p <- 1000
x <- rmvnorm(n, rep(0, p), diag(1,p))

highscore(x, "gaussian", lambda = 0.5, centered = TRUE, tol = 1e-06, maxit = 1000)
}
