## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
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


### Gumbel Distribution ########################################################

## distribution function
pGumbel <- function(q, mu = 0, sigma = 1)
{
  stopifnot(sigma > 0)
  exp(-exp(-( (q-mu)/sigma )))
}

## quantile function
qGumbel <- function(p, mu = 0, sigma = 1)
{
  stopifnot(0 < p, p < 1, sigma > 0)
  mu + sigma * (-log(-log(p)))
}

## density
dGumbel <- function(x, mu = 0, sigma = 1, log = FALSE)
{
  stopifnot(sigma > 0)
  q <- (x - mu) / sigma
  out <- -q - exp(-q) - log(sigma)
  if(!log) out <- exp(out)
  out
}

## random number generation
rGumbel <- function(n, mu = 0, sigma = 1)
{
  stopifnot(sigma > 0)
  qGumbel(runif(n), mu, sigma)
}
