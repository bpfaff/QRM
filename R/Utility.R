##
Pconstruct <- function(theta){
  n <- length(theta)
  d <- (1 + sqrt(1 + 8 * n)) / 2
  A <- matrix(0, nrow = d, ncol = d)
  A[lower.tri(A)] <- theta
  diag(A) <- 1
  Q <- A %*% t(A)
  P <- cov2cor(Q)
  P
}
##
Pdeconstruct <- function(P){
  A <- t(chol(P))
  Adiag <- diag(diag(A))
  Astar <- solve(Adiag) %*% A
  Astar[lower.tri(Astar)]
}
## Empirical (cumulative) distribution function
edf <- function(v, adjust = FALSE){
  original <- v
  v <- sort(v)
  vv <- cumsum(!duplicated(v))
  repeats <- tapply(v, v, length)
  add <- rep(cumsum(repeats - 1.), repeats)
  df <- (vv + add)/(length(vv) + as.numeric(adjust))
  as.numeric(df[rank(original)])
}
##
ESnorm <- function(p, mu = 0, sd = 1){
  mu + sd * dnorm(qnorm(p)) / (1 - p)
}
##
ESst <- function(p, mu = 0, sd = 1, df, scale = FALSE){
  ES <- (dt(qt(p, df), df)/(1 - p)) * ((df + (qt(p, df))^2) / (df - 1))
  if (scale) ES <- ES * sqrt((df - 2) / df)
  mu + sd * ES
}
