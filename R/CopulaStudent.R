## Random variates of Student's t Copula
rcopula.t <- function(n, df, Sigma){
  d <- dim(Sigma)[1]
  diagvals <- diag(Sigma)
  if(!(all(diagvals == 1))) stop("Sigma should be correlation matrix")
  tmp <- rmt(n, df, Sigma, mu = 0)
  matrix(pt(tmp, df), ncol = d)
}
## Density of Student's t copula
dcopula.t <- function(Udata, df, Sigma, log = FALSE){
  d <- dim(Udata)[2]
  Qdata <- apply(Udata, 2, qt, df = df)
  out <- dmt(Qdata, df, rep(0, d), Sigma, log = TRUE) - apply(log(apply(Qdata, 2, dt, df = df)), 1, sum)
  if(!(log)) out <- exp(out)
  out
}
## Fitting
fit.tcopula <- function(Udata, method = c("all", "Kendall", "Spearman"), startdf = 5, ...){
  method <- match.arg(method)
  if(method == "all"){
    negloglik1 <- function(theta, data){
      nu <- theta[1]
      P <- Pconstruct(theta[-1])
      -sum(dcopula.t(data, nu, P, log = TRUE))
    }
  } else {
    negloglik2 <- function(theta, data, P){
      -sum(dcopula.t(data, theta, P, log = TRUE))
    }
  }
  if(method == "Kendall"){
    Rtau <- Kendall(Udata)
    P <- sin(pi * Rtau / 2)
  } else {
    P <- Spearman(Udata)
  }
  if(min(eigen(P)$values) < 0) stop("Non psd covariance matrix")
  if(method == "all"){
    theta <- c(startdf, Pdeconstruct(P))
    fit <- nlminb(theta, negloglik1, data = Udata, ...)
    P <- Pconstruct(fit$par[-1])
  } else {
    fit <- nlminb(startdf, negloglik2, data = Udata, P = P, ...)
  }
  nu <- fit$par[1]
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  list(P = P, nu = nu, converged = conv, ll.max = -fit$objective, fit = fit)
}


