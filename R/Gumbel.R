## Copyright (C) 2013 Marius Hofert, Berhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


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
