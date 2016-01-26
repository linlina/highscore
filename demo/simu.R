library(igraph)
library(mvtnorm)
library(tmvtnorm)


cat("Gaussian demo \n")
side <- 10
p <- side^2
gr <- graph.lattice(c(side, side))
Adj <- as.matrix(get.adjacency(gr))
iCor <- matrix(0, p, p)
iCor[which(Adj > 0)] <- runif(length(which(Adj > 0)), 0.5, 1)
for (i in 1:p){
  if (sum(abs(iCor[i,])) != 0) iCor[i, ] <- iCor[i, ]/sum(abs(iCor[i,]) * 1.5)
}
iCor <- (iCor + t(iCor))/2
iCor[which(Adj > 0)] <- iCor[which(Adj > 0)]
diag(iCor) <- 1
Sigma <- solve(iCor)
Sigma <- cov2cor(Sigma)
Prec <- solve(Sigma)
Prec[which(abs(Prec) < 1e-05)] <- 0

n <- 1000
x <- rmvnorm(n, rep(0, p), Sigma)
set.seed(1)
test1 <- highscore(x, "gaussian", lambda = 0.5, centered = TRUE, tol = 1e-06, maxit = 1000)


cat("non-negative Gaussian demo \n")
gr <- erdos.renyi.game(side^2, 0.1)
Adj <- as.matrix(get.adjacency(gr))
Cor <- matrix(0, p, p)
iCor[which(Adj > 0)] <- runif(length(which(Adj > 0)), 0.5, 1)
signs <- sample(c(-1, 1), length(which(Adj > 0)), replace = TRUE)
for (i in 1:p){
  if (sum(abs(iCor[i,])) != 0) iCor[i, ] <- iCor[i, ]/sum(abs(iCor[i,]) * 1.5)
}
iCor <- (iCor + t(iCor))/2
diag(iCor) <- 1
Sigma <- solve(iCor)
Sigma <- cov2cor(Sigma)
Prec <- solve(Sigma)
Prec[which(abs(Prec) < 1e-05)] <- 0

n <- 1000
set.seed(1)
x <- rtmvnorm(n, mean = rep(0, nrow(Sigma)), sigma = Sigma, lower = rep(0, nrow(Sigma)), upper = rep(Inf, nrow(Sigma)), algorithm = "gibbs", burn.in.samples = 100, thinning = 10)
test2 <- highscore(x, "nonnegative", lambda = 0.5, centered = TRUE, tol = 1e-06, maxit = 1000)


cat("conditional Gaussian demo \n")

cat("function for simulating data \n")

gibbs.cond <- function(n, A, B, d, e, burnin, thin) {
  p <- ncol(A)
  samplesize <- iter <- 0
  x <- rep(0, p)
  x.track <- NULL

  while (samplesize < n) {
    for (i in 1:p) {
      v <- -2 * (2 * crossprod(A[i, -i], x[-i]^2) + d[i])
      #print(i)
      #print(v)
      m <- -1/2 * (B[i, -i] * x[-i] + e[i])
      x[i] <- rnorm(1, m, sqrt(v))
      #print(x[i])
    }
    iter <- iter + 1
    if ((iter%%thin == 0) && (iter > burnin)) {
      x.track <- rbind(x.track, x)
      samplesize <- samplesize + 1
    }
  }
  r <- list(x = x.track, n = n, p = p, A = A, B = B, d = d, e = e)
  return(r)
}


p <- 50
gr <- erdos.renyi.game(p, 0.04)
Adj <- as.matrix(get.adjacency(gr))
d <- rep(-1, p)
A <- Adj * -1/25 #Adj * -1/50
B <- matrix(0, nrow = p, ncol = p)
e <- rep(8/50, p)#rep(8/50, p)

n <- 1000
set.seed(1)
x <- gibbs.cond(n, A, B, d, e, 100, 10)$x
test3 <- highscore(x, "conditional", lambda = 0.2, centered = TRUE, tol = 1e-06, maxit = 1000)

