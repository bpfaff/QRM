## Random variates of Gauss Copula
rcopula.gauss <- function(n, Sigma){
  d <- dim(Sigma)[1]
  diagvals <- diag(Sigma)
  if(!(all(diagvals == 1))) stop("Sigma should be correlation matrix")
  mnorm <- rmnorm(n, Sigma = Sigma, mu = 0)
  matrix(pnorm(mnorm), ncol = d)
}
## Density of Gauss Copula
dcopula.gauss <- function(Udata, Sigma, log = FALSE){
  d <- dim(Udata)[2]
  Qdata <- apply(Udata, 2, qnorm)
  out <- dmnorm(Qdata, rep(0, d), Sigma, log = TRUE) - apply(log(apply(Qdata, 2, dnorm)), 1, sum)
  if(!(log)) out <- exp(out)
  out
}
## Fitting of Gauss Copula
fit.gausscopula <- function(Udata, ...){
  negloglik <- function(theta, data){
    Sigma <- Pconstruct(theta)
    -sum(dcopula.gauss(data, Sigma, log = TRUE))
  }
  theta <- Pdeconstruct(Spearman(Udata))
  fit <- nlminb(theta, negloglik, data = Udata, ...)
  theta <- fit$par
  Sigma <- Pconstruct(theta)
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  list(P = Sigma, converged = conv, ll.max = -fit$objective, fit = fit)
}
