#' Creates a grain size distribution
#'
#' Creates a grain size distribution, given a set of grain sizes, D50, and
#' sp.

#' @param D50 Median grain size (mm).
#' @param sp Spread of distribution (default is 1).
#' @param ds Vector of grain sizes to map the gsd to (mm).
#' @param plot Should the GSD be plotted (default is `TRUE``).
#'
#' @return The size fraction for each grain size.
#'
#' @export

gsd_maker <- function(D50, sp, ds, plot = TRUE){
  par(mfrow = c(1,1), mar = c(4, 4, 0.5, 0.5), oma = rep(0, 4))
  gsd <- stochasim::sim_gsd(D50 = D50, sp = sp, Map = plot)
  cdf <- stats::approx(gsd$size_class, gsd$cdf, xout = ds)$y
  cdf[which(is.na(cdf[1:3]))] <- 0
  cdf[which(is.na(cdf))] <- 1
  cdf[length(cdf)] <- 1
  if (plot){
    lines(ds, cdf, col = "red")
    points(ds, cdf, pch = 16, col = "red")
  }
  ps <- diff(c(0,cdf))

  return(ps)
}

#' Calculates grain size statistic from given distribution
#'
#' Calculates Dx from a grain size distribtion, where `x` is the fraction
#' of the grain size distribution finer than the calculated grain size (e.g. `x
#' = 0.5` for D50).
#'
#' @param ps Grain size fractions
#' @param Ds Grain sizes (mm)
#' @param x Percetile to find (e.g. 0.5 for D50)
#'
#' @return Size of Dx
#'
calc_Dx <- function(ps, Ds, x){

  #convert to cdf
  cdf <- cumsum(as.numeric(ps))

  #convert Ds to log2
  Ds_log <- log2(Ds)

  #Find Dx in log space
  Dx_log <- stats::approx(cdf, Ds_log, xout = x)$y

  #Reconvert to mm
  Dx <- 2 ^ Dx_log

  return(Dx)
}

#' Calculates geometric standard deviation of grain size distribtion
#'
#' @param ps Grain size fractions
#' @param Ds Grain sizes (mm)
#'
#' @return Geometric standard deviation
#'
calc_sigmag <- function(ps, Ds){

  #convert Ds to log2
  phi <- log2(Ds)

  #find geometric mean
  phi_m <- sum(ps * phi)

  #find geometric sd
  phi_sd <- sum((phi - phi_m) ^ 2 * ps)

  #convert to mm
  sd <- 2 ^ sqrt(phi_sd)

  return(sd)
}
