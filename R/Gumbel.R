## Gumbel Distribution
pGumbel <- function(q, mu = 0, sigma = 1){
  if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!")
  q <- (q - mu) / sigma
  exp(-exp(-q))
}
qGumbel <- function(p, mu = 0, sigma = 1){
  if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!")
  if(length(p[(p > 0) & (p < 1)]) < length(p)) stop("p parameter does not represent probabilities")
  mu + sigma * (-log(-log(p)))
}
dGumbel <- function(x, mu = 0, sigma = 1, log = FALSE){
  if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!")
  q <- (x - mu) / sigma
  out <- -q - exp(-q) - log(sigma)
  if(!(log)) out <- exp(out)
  out
}
rGumbel <- function(n, mu = 0, sigma = 1){
  if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!")
  U <- runif(n)
  qGumbel(U, mu, sigma)
}
