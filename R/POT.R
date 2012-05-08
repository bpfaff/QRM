fit.GPD <- function(data, threshold = NA, nextremes = NA, type = c("ml", "pwm"), information = c("observed", "expected"), optfunc = c("optim", "nlminb"), ...){
  type <- match.arg(type)
  optfunc <- match.arg(optfunc)
  information <- match.arg(information)
  if(is.na(nextremes) & is.na(threshold))
    stop("Enter either a threshold or the number of upper extremes")
  if(!is.na(nextremes) & !is.na(threshold))
    stop("Enter EITHER a threshold or the number of upper extremes")
  if(is.timeSeries(data)) data <- series(data) 
  data <- as.numeric(data)
  n <- length(data)
  if(!is.na(nextremes)){
    threshold <- findthreshold(data, nextremes)
  }
  exceedances <- data[data > threshold]
  excess <- exceedances - threshold
  Nu <- length(excess)
  if(type == "ml"){
    xbar <- mean(excess)
    a0 <- xbar
    gamma <- -0.35
    delta <- 0.0
    pvec <- ((1.:Nu) + delta)/(Nu + delta)
    a1 <- mean(sort(excess) * (1. - pvec))
    xi <- 2. - a0/(a0 - 2. * a1)
    beta <- (2. * a0 * a1)/(a0 - 2. * a1)
    par.ests <- c(xi,beta)
    negloglik <- function(theta, ydata){
      -sum(dGPD(ydata, theta[1], abs(theta[2]), log = TRUE))
    }
    deriv <- function(theta, ydata){
      xi <- theta[1]
      beta <- theta[2]
      term1 <- sum(ydata / (beta + xi * ydata))
      term2 <- sum(log(1 + xi * ydata / beta))
      d1 <- -term2 * xi^(-2) + (1 + 1 / xi) * term1
      d2 <- (length(ydata) - (xi + 1) * term1) / beta
      c(d1, d2)
    }
    if(optfunc == "optim"){
      fit <- optim(par.ests, fn = negloglik, gr = deriv, ydata = excess, ...)
    }
    if(optfunc == "nlminb"){
      fit <- nlminb(start = par.ests, objective = negloglik, gradient = deriv, ydata = excess, ...)
    }
    par.ests <- fit$par
    par.ests[2] <- abs(par.ests[2])
    ifelse(fit$convergence == 0, converged <- TRUE, converged <- FALSE)
    ll.max <- -negloglik(fit$par, ydata = excess)
    if(information == "observed"){
      fisher <- hessian(negloglik, fit$par, ydata = excess)
      varcov <- solve(fisher)
    }
    if(information == "expected"){
      one <- (1. + par.ests[1.])^2. / Nu
      two <- (2. * (1. + par.ests[1.]) * par.ests[2.]^2.) / Nu
      cov <-  - ((1. + par.ests[1.]) * par.ests[2.]) / Nu
      varcov <- matrix(c(one, cov, cov, two), 2.)
    }
  }
  if(type == "pwm"){
    xbar <- mean(excess)
    a0 <- xbar
    gamma <- -0.35
    delta <- 0.0
    pvec <- ((1.:Nu) + delta)/(Nu + delta)
    a1 <- mean(sort(excess) * (1. - pvec))
    xi <- 2. - a0/(a0 - 2. * a1)
    beta <- (2. * a0 * a1)/(a0 - 2. * a1)
    par.ests <- c(xi, beta)
    denom <- Nu * (1. - 2. * xi) * (3. - 2. * xi)
    if(xi > 0.5){
      denom <- NA
      warning("Asymptotic standard errors not available for PWM Type when xi > 0.5")
    }
    one <- (1. - xi) * (1. - xi + 2. * xi^2.) * (2. - xi)^2.
    two <- (7. - 18. * xi + 11. * xi^2. - 2. * xi^3.) * beta^2.
    cov <- beta * (2. - xi) * (2. - 6. * xi + 7. * xi^2. - 2. *xi^3.)
    varcov <- matrix(c(one, cov, cov, two), 2.)/denom
    information <- "expected"
    converged <- NA
    ll.max <- NA
  }
  par.ses <- sqrt(diag(varcov))
  p.less.thresh <- 1. - Nu / n
  out <- list(n = length(data), data = exceedances, threshold = threshold,
              p.less.thresh = p.less.thresh, n.exceed = Nu, type = type,
              par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              information = information, converged = converged, ll.max = ll.max)
  names(out$par.ests) <- c("xi", "beta")
  names(out$par.ses) <- c("xi", "beta")
  out
}
## POT: Threshold
findthreshold <- function(data, ne){
  if(is.timeSeries(data)) data <- as.vector(series(data))
  if(!is.vector(data)) stop("data input to findthreshold() must be a vector or timeSeries with only one data column")
  if(all(length(data) < ne)) stop("data length less than ne (number of exceedances")
  data <- rev(sort(as.numeric(data)))
  thresholds <- unique(data)
  indices <- match(data[ne], thresholds)
  indices <- pmin(indices + 1., length(thresholds))
  thresholds[indices]
}
## Tail plot
plotTail <- function(object, extend = 2, fineness = 1000, ...){
  data <- as.numeric(object$data)
  threshold <- object$threshold
  xi <- object$par.ests[names(object$par.ests) == "xi"]
  beta <- object$par.ests[names(object$par.ests) == "beta"]
  xpoints <- sort(data)
  ypoints <- ppoints(sort(data))
  xmax <- max(xpoints) * extend
  prob <- object$p.less.thresh
  ypoints <- (1- prob) * (1 - ypoints)
  x <- threshold + qGPD((0:(fineness - 1)) / fineness, xi, beta)
  x <- pmin(x,xmax)
  y <- pGPD(x - threshold, xi, beta)
  y <- (1 - prob) * (1 - y)
  plot(xpoints, ypoints, xlim = range(threshold, xmax), ylim = range(ypoints, y), xlab = "x (on log scale)", ylab = "1-F(x) (on log scale)", log = "xy", ...)
  lines(x, y)
  return(invisible(list(xpoints = xpoints, ypoints = ypoints, x = x, y = y)))
}
## Risk Measures
showRM <- function(object, alpha, RM = c("VaR", "ES"), extend=2, ci.p = 0.95, like.num = 50., ...){
  threshold <- object$threshold
  par.ests <- object$par.ests
  xihat <- par.ests[names(par.ests) == "xi"]
  betahat <- par.ests[names(par.ests) == "beta"]
  p.less.thresh <- object$p.less.thresh
  a <- (1 - alpha) / (1 - p.less.thresh)
  quant <- threshold + betahat * (a^(-xihat) - 1) / xihat
  es <- quant / (1 - xihat) + (betahat - xihat * threshold) / (1 - xihat)
  point.est <- switch(RM, VaR = quant, ES = es)
  plotTail(object, extend = 2)
  abline(v = point.est)
  xmax <- max(object$data) * extend
  parloglik <- function(theta, excessesIn, xpiIn, aIn, uIn, RMIn){
    xi <- theta
    if(RMIn == "VaR") beta <- xi * (xpiIn - uIn) / (aIn^(-xi) - 1)
    if(RMIn == "ES")  beta <- ((1 - xi) * (xpiIn - uIn))/(((aIn^( - xi) - 1) / xi) + 1)
    if(beta <= 0){
      out <- 1.0e17
    } else {
      out <- -sum(dGPD(excessesIn, xi, beta, log = TRUE))
    }
    out
  }
  parmax <- NULL
  start <- switch(RM, VaR = threshold, ES = quant)
  xp <- exp(seq(from = log(start), to = log(xmax), length = like.num))
  for(i in 1.:length(xp)){
    optimfit2 <- optim(xihat, parloglik, excessesIn=(object$data - threshold), xpiIn=xp[i], 
                       aIn=a, uIn=threshold, RMIn=RM, ...)
    parmax <- rbind(parmax, -parloglik(optimfit2$par, excessesIn = (object$data - threshold), xpiIn = xp[i], aIn = a, uIn = threshold, RMIn = RM))
  }
  overallmax <-  object$ll.max
  crit <- overallmax - qchisq(0.999, 1) / 2.
  cond <- parmax > crit
  xp <- xp[cond]
  parmax <- parmax[cond]
  par(new = TRUE)
  plot(xp, parmax, type = "n", xlab = "", ylab = "", axes = FALSE,
       ylim = range(overallmax, crit), xlim = range(threshold, xmax), log = "x")
  axis(4, at = overallmax - qchisq(c(0.95, 0.99), 1.)/2., labels = c("95", "99"), tick = TRUE)
  aalpha <- qchisq(ci.p, 1.)
  abline(h = overallmax - aalpha / 2, lty = 2, col = 2)
  cond <- !is.na(xp) & !is.na(parmax)
  smth <- spline(xp[cond], parmax[cond], n = 200.)
  lines(smth, lty = 2., col = 2.)
  ci <- smth$x[smth$y > overallmax - aalpha/2.]
  out <- c(min(ci), point.est, max(ci))
  names(out) <- c("Lower CI", "Estimate", "Upper CI")
  out
}
## ME plot
MEplot <- function(data, omit = 3., labels = TRUE, ...){
  if(is.timeSeries(data)) data <- series(data)   
  data <- as.numeric(data)
  n <- length(data)
  myrank <- function(x, na.last = TRUE){
    ranks <- sort.list(sort.list(x, na.last = na.last))
    if(is.na(na.last)) x <- x[!is.na(x)]
    for(i in unique(x[duplicated(x)])){
      which <- x == i & !is.na(x)
      ranks[which] <- max(ranks[which])
    }
    ranks
  }
  data <- sort(data)
  n.excess <- unique(floor(length(data) - myrank(data)))
  points <- unique(data)
  nl <- length(points)
  n.excess <- n.excess[ - nl]
  points <- points[ - nl]
  excess <- cumsum(rev(data))[n.excess] - n.excess * points
  y <- excess/n.excess
  plot(points[1.:(nl - omit)], y[1.:(nl - omit)], xlab = "", ylab = "",...)
  if(labels) title(xlab = "Threshold", ylab = "Mean Excess")
  return(invisible(list(x = points[1.:(nl - omit)], y = y[1.:(nl - omit)])))
}
## Risk Measures
RiskMeasures <- function(out, p){
  u <- out$threshold
  par.ests <- out$par.ests
  xihat <- par.ests[names(par.ests) == "xi"]
  betahat <- par.ests[names(par.ests) == "beta"]
  p.less.thresh <- out$p.less.thresh
  lambda <- 1. / (1. - p.less.thresh)
  quant <- function(pp, xi, beta, u, lambda){
    a <- lambda * (1. - pp)
    u + (beta * (a^( - xi) - 1.))/xi
  }
  short <- function(pp, xi, beta, u, lambda){
    a <- lambda * (1. - pp)
    q <- u + (beta * (a^( - xi) - 1.))/xi
    (q * (1. + (beta - xi * u)/q))/(1. - xi)
  }
  q <- quant(p, xihat, betahat, u, lambda)
  es <- short(p, xihat, betahat, u, lambda)
  cbind(p, quantile = q, sfall = es)
}
## GPD shape varies with threshold or number of extremes. 
xiplot <- function(data, models = 30., start = 15., end = 500., reverse = TRUE,
                   ci = 0.95, auto.scale = TRUE, labels = TRUE, table = FALSE, ...){
  if(is.timeSeries(data)) data <- series(data)
  data <- as.numeric(data)
  qq <- 0.
  if(ci) qq <- qnorm(1. - (1. - ci)/2.)
  x <- trunc(seq(from = min(end, length(data)), to = start, length = models))
  gpd.dummy <- function(nex, data){
    out <- fit.GPD(data = data, nextremes = nex, information = "expected")
    c(out$threshold, out$par.ests[1.], out$par.ses[1.])
  }
  mat <- apply(as.matrix(x), 1., gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("threshold", "shape", "se", "exceedances"), NULL)
  thresh <- mat[1,  ]
  y <- mat[2,  ]
  yrange <- range(y)
  if(ci){
    u <- y + mat[3.,  ] * qq
    l <- y - mat[3.,  ] * qq
    yrange <- range(y, u, l)
  }
  index <- x
  if(reverse) index <-  - x
  if(auto.scale){
    plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  } else {
    plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  }
  axis(1, at = index, labels = paste(x), tick = FALSE)
  axis(2)
  axis(3, at = index, labels = paste(format(signif(thresh, 3.))), tick = FALSE)
  box()
  if(ci){
    lines(index, u, lty = 2., col = 2.)
    lines(index, l, lty = 2., col = 2.)
  }
  if(labels){
    labely <- "Shape (xi)"
    if(ci) labely <- paste(labely, " (CI, p = ", ci, ")", sep = "")
    title(xlab = "Exceedances", ylab = labely)
    mtext("Threshold", side = 3., line = 3.)
  }
  if(table) print(mat)
  return(invisible(list(x = index, y = y, upper = u, lower = l)))
}
## Hill Plot
hillPlot <- function (data, option = c("alpha", "xi", "quantile"), start = 15, 
    end = NA, reverse = FALSE, p = NA, ci = 0.95, auto.scale = TRUE, labels = TRUE, ...){
  if(is.timeSeries(data)) data <- as.vector(series(data))
  data <- as.numeric(data)
  ordered <- rev(sort(data))
  ordered <- ordered[ordered > 0]
  n <- length(ordered)
  option <- match.arg(option)
  if((option == "quantile") && (is.na(p)))
    stop("Input a value for the probability p")
  if ((option == "quantile") && (p < 1 - start/n)){
    cat("Graph may look strange !! \n\n")
    cat(paste("Suggestion 1: Increase `p' above", format(signif(1 - start/n, 5)), "\n"))
    cat(paste("Suggestion 2: Increase `start' above ", ceiling(length(data) * (1 - p)), "\n"))
  }
  k <- 1:n
  loggs <- logb(ordered)
  avesumlog <- cumsum(loggs)/(1:n)
  xihat <- c(NA, (avesumlog - loggs)[2:n])
  alphahat <- 1/xihat
  y <- switch(option, alpha = alphahat, xi = xihat, quantile =
              ordered * ((n * (1 - p))/k)^(-1/alphahat))
  ses <- y/sqrt(k)
  if(is.na(end)) end <- n
  x <- trunc(seq(from = min(end, length(data)), to = start))
  y <- y[x]
  ylabel <- option
  yrange <- range(y)
  if(ci && (option != "quantile")){
    qq <- qnorm(1 - (1 - ci)/2)
    u <- y + ses[x] * qq
    l <- y - ses[x] * qq
    ylabel <- paste(ylabel, " (CI, p =", ci, ")", sep = "")
    yrange <- range(u, l)
  }
  if(option == "quantile") ylabel <- paste("Quantile, p =", p)
  index <- x
  if(reverse) index <- -x
  if(auto.scale){
    plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  } else {
    plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
  }
  axis(1, at = index, labels = paste(x), tick = FALSE)
  axis(2)
  threshold <- findthreshold(data, x)
  axis(3, at = index, labels = paste(format(signif(threshold, 3))), tick = FALSE)
  box()
  if(ci && (option != "quantile")){
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if(labels){
    title(xlab = "Order Statistics", ylab = ylabel)
    mtext("Threshold", side = 3, line = 3)
  }
  return(invisible(list(x = index, y = y)))
}
## QQ-Plot for GPD
plotFittedGPDvsEmpiricalExcesses <- function(data, threshold = NA, nextremes = NA){
  if(is.na(nextremes) & is.na(threshold))
    stop("Enter either a threshold or the number of upper extremes")
  if(is.timeSeries(data)) data <- series(data)
  mod <- fit.GPD(data, threshold, nextremes)
  if(!is.na(nextremes)) threshold <- findthreshold(as.vector(data), nextremes)
  pECDF <- edf(mod$data)
  maxVal <- as.numeric(max(data))
  quantVector <- seq(threshold, maxVal, 0.25)
  pG <- pGPD(quantVector-threshold, mod$par.ests["xi"], mod$par.ests["beta"]) 
  split.screen(c(1, 1))
  plot(quantVector, pG, type = "l", log = "x", xlab = "x (log scale)", ylab = "Fu(x-u)")
  screen(1, new = FALSE)
  plot(mod$data, pECDF, log = "x", pch = 19, xlab = "x (log scale)", ylab= "Fu(x-u)", col = "blue")
  close.screen(all.screens = TRUE)
}



