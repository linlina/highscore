\name{hello}
\alias{hello}
\title{Hello, World!}
\usage{
highscore(x, model, lambda = 0, centered = TRUE, tol = 1e-6, maxit = 1000, ...)
}
\arguments{
  \item{x}{The n x p data matrix, that is, n observations and p variables}
  \item{model}{Either "gaussian", "nonnegative" (for non-negative Gaussians), "conditional" (for conditional Gaussian)}

  \item{lambda}{The l1 regularization parameter}

  \item{centered}{Either TRUE or FALSE: (1) For Gaussian and non-negative Gaussian, this means we estimate the mean/shift parameter (2) For conditional Gaussian, this means we estimate the matrix B for x_ix_j cross terms}

  \item{tol}{Convergence threshold}

  \item{maxit}{Maximum number of iterations before
  termination}

  \item{save.iterates}{Returns iterates if TRUE}

  \item{...}{ignored}
}
\description{

}
\examples{
hello()
}
