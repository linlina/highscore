
highscore <- function(x, model, lambda = 0, centered = TRUE, tol = 1e-6, maxit = 1000){
  n <- nrow(x)
  p <- ncol(x)

  if (!model %in% c("gaussian", "nonnegative", "conditional")){
    stop("We haven't implemented the model you've picked!")
  }

  orig_lambda <- lambda

  if (model == "gaussian"){
    if (centered){
      if (n >= p) {
        S <- crossprod(x, x)/n
      } else {
        S <- diag(colSums(x * x)/n)
      }
    } else {
      if (n >= p) {
        S <- cov(x)
      } else {
        stdx <- apply(x, 2, scale, scale = FALSE, center = TRUE)
        if (n >= p) {
          S <- crossprod(stdx, stdx)/n
        } else {
          S <- diag(colSums(stdx * stdx)/n)
        }
      }
    }
    K <- diag(1, p)
    if (n >= p){
      temp <- .C("Score1", p = as.integer(p), S = as.double(S), lambda = as.double(lambda), K = as.double(K), tol = as.double(tol), maxit = as.integer(maxit), PACKAGE = "highscore")
    } else {
      temp <- .C("Score2", n = as.integer(n), p = as.integer(p), X = as.double(x), S = as.double(S), lambda = as.double(lambda), K = as.double(K), tol = as.double(tol), maxit = as.integer(maxit), PACKAGE = "highscore")
    }
    temp$K <- matrix(temp$K, p, p)
    if (centered){
      output <- list(K = temp$K, lambda = orig_lambda, iter = temp$maxit, tol = tol)
    } else {
      mu <- colMeans(x)
      output <- list(K = temp$K, mu = mu, lambda = orig_lambda, iter = temp$maxit, tol = tol)
    }
  } else if (model == "nonnegative") {
    K <- diag(1, p)
    if (centered){
      S1 <- crossprod(x, x)/n
      S2 <- matrix(0, nrow = p, ncol = p*p)
      for (k in 1:p){
        S2[ , ((k-1)*p + 1):(k*p)] <- 1/n*crossprod(sweep(x, MARGIN = 1, x[,k]^2, `*`), x)
      }
      temp <- .C("ScoreNN1", p = as.integer(p), S1 = as.double(S1), S2 = as.double(S2), lambda = as.double(lambda), K = as.double(K), tol = as.double(tol), maxit = as.integer(maxit), package = "highscore")
      temp$K <- matrix(temp$K, p, p)
      output <- list(K = temp$K, lambda = orig_lambda, iter = temp$maxit, tol = tol)
    } else {
      S1 <- crossprod(x, x)/n
      S2 <- crossprod(x^2, x)/n
      S3 <- matrix(0, nrow = p, ncol = p*p)
      for (k in 1:p){
        S3[ , ((k-1)*p + 1):(k*p)] <- 1/n*crossprod(sweep(x, MARGIN = 1, x[,k]^2, `*`), x)
      }
      xbar <- colSums(x)
      b <- rep(0, p)
      temp <- .C("ScoreNN2", p = as.integer(p), Xbar = as.double(xbar), S1 = as.double(S1), S2 = as.double(S2), S3 = as.double(S3), lambda = as.double(lambda), K = as.double(K), b = as.double(b), tol = as.double(tol), maxit = as.integer(maxit), package = "highscore")
      temp$K <- matrix(temp$K, p, p)
      output <- list(K = temp$K, mu = mu, lambda = orig_lambda, iter = temp$maxit, tol = tol)
    }
  } else if (model == "conditional"){
    if (centered){
      S <- crossprod(x, x)/n
      SS <- crossprod(x^2, x)/n
      SD <- crossprod(x^2, x^2)/n
      SDQ <- crossprod(x^4, x^2)/n
      ST <- matrix(0, nrow = p, ncol = p*p)
      for (k in 1:p){
        ST[ , ((k-1)*p + 1):(k*p)] <- 1/n*crossprod(sweep(x^2, MARGIN = 1, x[,k]^2, `*`), x^2)
      }
      xbar <- colMeans(x)
      A <- matrix(0, p, p)
      d <- rep(1, p)
      e <- rep(0, p)
      temp <- .C("Score_cond", p = as.integer(p), Xbar = as.double(xbar), S = as.double(S), SS = as.double(SS), SD = as.double(SD), SDQ = as.double(SDQ), ST = as.double(ST), lambda = as.double(lambda), A = as.double(A), d = as.double(d), e = as.double(e), tol = as.double(tol), maxit = as.integer(maxit), package = "highscore")
      temp$A <- matrix(temp$A, p, p)
      output <- list(A = temp$A, d = temp$d, e = temp$e, lambda = orig_lambda, iter = temp$maxit, tol = tol)
    } else {
      S <- crossprod(x, x)/n
      SS <- crossprod(x^2, x)/n
      SD <- crossprod(x^2, x^2)/n
      SDQ <- crossprod(x^4, x^2)/n
      ST <- matrix(0, nrow = p, ncol = p*p)
      SSS <- matrix(0, nrow = p, ncol = p*p)
      for (k in 1:p){
        ST[ , ((k-1)*p + 1):(k*p)] <- 1/n*crossprod(sweep(x^2, MARGIN = 1, x[,k]^2, `*`), x^2)
        SSS[ , ((k-1)*p + 1):(k*p)] <- 1/n*crossprod(sweep(x, MARGIN = 1, x[,k]^2, `*`), x)
      }
      xbar <- colMeans(x)
      A <- matrix(0, p, p)
      B <- matrix(0, p, p)
      d <- rep(1, p)
      e <- rep(0, p)
      temp <- .C("Score_cond2", p = as.integer(p), Xbar = as.double(xbar), S = as.double(S), SS = as.double(SS), SSS = as.double(SSS), SD = as.double(SD), SDQ = as.double(SDQ), ST = as.double(ST), lambda = as.double(lambda), A = as.double(A), B = as.double(B), d = as.double(d), e = as.double(e), tol = as.double(tol), maxit = as.integer(maxit), package = "highscore")
      temp$A <- matrix(temp$A, p, p)
      output <- list(A = temp$A, B= temp$B, d= temp$d, e = temp$e, lambda = orig_lambda, iter = temp$maxit, tol = tol)
    }
  } else {
    cat("Something went wrong, oh no", "\n")
    output <- 0
  }
  return(output)
}
