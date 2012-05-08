## GPD Distribution
pGPD <- function(q, xi, beta = 1){
  x <- q / beta
  if(xi == 0){
    out <- pexp(x)
  } else {
    out <- 1 -(1+xi*x)^(-1/xi)
  }
  out[x < 0] <- 0
  if(xi < 0) out[x > -1 / xi] <- 1
  out 
}
qGPD <- function(p, xi, beta = 1){
  if(xi == 0){
    out <- qexp(p)
  } else {
    inner <- (1 / xi) * ((1 - p)^(-xi) - 1)
    out <- pmax(inner, 0)
    if (xi < 0) out <- pmin(inner, 1 / abs(xi))
  }
  beta * out
}
rGPD <- function(n, xi, beta = 1){
    U <- runif(n)
    qGPD(U, xi, beta)
}
dGPD <- function(x, xi, beta = 1, log = FALSE){
  xx <- x / beta
  if(xi == 0){
    out <- log(dexp(xx)) - log(beta)
  } else {
    out <- rep(-Inf, length(x))
    cond <- (xx > 0)
    if (xi < 0) cond <- cond & (xx < 1 / abs(xi))
    out[cond] <- (-1 / xi-1) * log(1 + xi * xx[cond]) - log(beta)
    }
    if(!log) out <- exp(out)
    out
}
