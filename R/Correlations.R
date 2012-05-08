## Function for computing equal correlations
equicorr <- function(d, rho){
  d <- as.integer(d)
  if(abs(rho) > 1)
    stop("Correlation is greater than unity in absolute value")
  if(rho < ( - (d - 1.0)^(-1.0)))
    stop(paste("rho must be at least",  - (d - 1.0)^(-1.0)))
  out <- matrix(rho, nrow = d, ncol = d)
  diag(out) <- 1
  out
}
## Make a matrix positive definite
eigenmeth <- function(mat, delta = 0.001){
  decomp <- eigen(mat)
  Lambda <- decomp$values
  Lambda[Lambda < 0] <- delta
  Gamma <- decomp$vectors
  newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
  D <- 1/sqrt(diag(newmat))
  diag(D) %*% newmat %*% diag(D)
}
## Spearman rank correlations
Spearman <- function(data, ...){
  out <- cor(data, method = "spearman", ...)
  out
}
## Kendall rank correlations
Kendall <- function(data, ...){
  out <- cor(data, method = "kendall", ...)
  out
}
