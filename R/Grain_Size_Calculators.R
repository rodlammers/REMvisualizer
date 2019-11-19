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

#' Calculates critical stream power for the initial grain size distribution, by reach
#'
#' @param path Path to folder with model outputs
#' @param return_vals Logical. If TRUE, returns calculated critical stream powers by grain size and reach
#'
#' @return Critical stream power by grain size
#'
calc_omega_c <- function(path, return_vals = FALSE){

  #Get omega_c star from model inputs
  inputs <- read.table(file.path(path, "Model Inputs.txt"), header = TRUE, sep = "\t")
  omegac_star <- inputs$omegac_star
  b <- inputs$b

  #Get grain sizes and proportions
  ds <- unlist(read.table(file.path(path, "Input ds.txt"), header = FALSE, sep = " "))
  ps <- read.table(file.path(path, "Input ps.txt"), header = FALSE, sep = " ")

  #Calculate D50
  D50 <- apply(ps, 1, calc_Dx, ds, 0.5)

  #Get omegac_star_adj
  omegac_star_adj <- sapply(D50, function(x, omegac_star, ds, b){
    omegac_star * (ds / x) ^ b
  }, omegac_star, ds, b)

  #calculate omegac (W/m2)
  omegac <- apply(omegac_star_adj, 2, function(x, ds){
    x * 1000 * (9.81 * 1.65 * ds) ^ 1.5
  }, ds)

  #Plot critical stream power by reach and grain size
  ds_matrix <- matrix(rep(ds, ncol(omegac)), ncol = ncol(omegac))
  colors <- cRamp_legend(ncol(ds_matrix), palette = "viridis")
  par(mfrow = c(1,1), mar = c(4, 4, 1, 0.5))
  plot(as.numeric(ds_matrix) * 1000, as.numeric(omegac), pch = 21, bg = rep(colors, each = nrow(ds_matrix)), log = "x",
       xlab = "Grain size [mm]", ylab = "Critical Stream Power [W/m2]", las = 1)

  if (return_vals){
    return(omegac)
  }

}
