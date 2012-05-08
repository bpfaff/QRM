##
## Functions for plotting
##
## Bivariate density plot
BiDensPlot <- function(func, xpts =  c(-2, 2), ypts = c(-2, 2), npts = 50, type = c("persp", "contour"), ...){
  type <- match.arg(type)
  npts <- as.integer(npts)
  x <- seq(from = xpts[1], to = xpts[2], length.out = npts)
  y <- seq(from = ypts[1], to = ypts[2], length.out = npts)
  xval <- rep(x, length(y))
  yval <- rep(y, rep(length(x), length(y)))
  data <- cbind(xval, yval)
  z <- eval(func(data, ...))
  z <- matrix(z, nrow = length(x), ncol = length(y), byrow = TRUE)
  switch(type, persp = persp(x, y, z), contour = contour(x, y, z))
  return(invisible(list(x = x, y = y, z = z)))
}
##
QQplot <- function(x, a = 0.5, reference = c("normal", "exp", "student"), ...){
  n <- length(x)
  reference <- match.arg(reference)
  plot.points <- ppoints(n, a)
  func <- switch(reference,
                 normal = qnorm,
                 exp = qexp,
                 student = qt)
  xp <- func(plot.points,...)
  y <- sort(x)
  plot(xp, y, xlab = paste("Theoretical", reference), ylab = "Empirical")
  return(invisible(list(x = x, y = y)))
}
