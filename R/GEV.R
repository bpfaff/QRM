## GEV Distribution
pGEV <- function(q, xi, mu = 0, sigma = 1){
  x <- (q - mu) / sigma
  if (xi==0){
    return(out <- pGumbel(q, mu, sigma))
  } else {
    out <- exp( - (1 + xi * x)^(-1 / xi))
  }
  if (xi > 0) out[1 + xi * x  <= 0] <- 0
  if (xi < 0) out[1 + xi * x <=0] <- 1
  out
}
qGEV <- function(p, xi, mu = 0, sigma = 1){
  if(xi==0){
    return(out <- qGumbel(p, mu, sigma))
  } else {
    inner <- ((-log(p))^(-xi) - 1) / xi
    if (xi > 0) out <- pmax(inner, -1 / xi)
    if (xi < 0) out <- pmin(inner, 1 / abs(xi))
  }
  mu + sigma * out
}
dGEV <- function(x, xi, mu = 0, sigma = 1, log = FALSE){
  xx <- (x - mu) / sigma
  if(xi==0){
    return(out <- dGumbel(x, mu, sigma, log))
  } else {
    out <- rep(-Inf,length(x))
    out[1 + xi * xx > 0] <- (-1 / xi - 1) * log(1 + xi * xx[1 + xi * xx > 0]) - (1 + xi * xx[1 + xi * xx > 0])^(-1 / xi) - log(sigma)
  }
  if (!log) out <- exp(out)
  out
}
rGEV <- function(n, xi, mu = 0, sigma = 1){
  U <- runif(n)
  qGEV(U, xi, mu, sigma)
}
fit.GEV <- function(maxima, ...){
  sigma0 <- sqrt((6. * var(maxima))/pi)
  mu0 <- mean(maxima) - 0.57722 * sigma0
  xi0 <- 0.1
  theta <- c(xi0, mu0, sigma0)
  negloglik <- function(theta, maxvalue){
    -sum(dGEV(maxvalue, theta[1], theta[2], abs(theta[3]), log = TRUE))
  }
  optimfit <- optim(theta, fn = negloglik, maxvalue = maxima, ...)
  par.ests <- optimfit$par
  ifelse(optimfit$convergence == 0, converged <- TRUE, converged <- FALSE)
  par.ests[3] <- abs(par.ests[3])
  fisher <- hessian(negloglik, par.ests, maxvalue=maxima)
  varcov <- solve(fisher)
  par.ses <- sqrt(diag(varcov))
  out <- list(par.ests = par.ests, par.ses = par.ses, varcov = varcov, converged = converged, llmax = -negloglik(par.ests,maxima))
  names(out$par.ests) <- c("xi", "mu", "sigma")
  names(out$par.ses) <- c("xi", "mu", "sigma")
  out
}
