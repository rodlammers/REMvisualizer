#' Calculates a mass balance of modeled channel change
#'
#' Calculates changes in cross section area and compares that to sediment inputs
#' and outputs to determine if mass was conserved during the simulation. Note
#' the mass balance is not accurate if meandering was simulated.
#'
#' @param path Path to folder with model outputs
#'
#' @return Prints results of the mass balance and a plot of volume changes of
#'   channel cross sections.
#'
#'   \tabular{ll}{ Volume sum \tab Sum of total channel volume change ((-)
#'   indicates net erosion, (+) indicates net aggradation)\cr Bed vol out \tab
#'   Total volume of bed material load explorted from the watershed\cr Bed vol
#'   in \tab Total volume of bed material load imported to watershed\cr Bank
#'   tank \tab Volume of failed bank material in the bank "tank" (this is
#'   material that couldn't be deposited at the bank toe but is stored in a
#'   virtual "tank" to conserve mass.)\cr Bank washload \tab Volume of eroded
#'   bank washload\cr Bed washload (cohesive) \tab Volume of eroded washload
#'   from cohesive bed erosion\cr Knickpoint washload \tab Volume of eroded
#'   washload from knickpoint erosion\cr Knickpoint correction \tab A volume
#'   correction for when a knickpoint is initially located between two cross
#'   sections\cr Diff \tab Calculated volume difference between calculated
#'   change in sediment inflow and outflow and total channel change\cr Percent
#'   Diff \tab Calculated volume difference as a percentage of Volume sum\cr }
#'
#' @importFrom utils read.table
#' @importFrom dplyr %>%
#' @importFrom graphics par plot points abline
#' @export
#'
XS_areas <- function(path = ""){

  XS <- read.table(paste0(path, "/Output XS geometry all.txt"), header = FALSE,
                   sep= " ") %>%
    dplyr::filter_(~ V1 == 0 | V1 == max(V1))

  # if (is.na(XS[1, ncol(XS)])){
  #   XS <- XS[,1:(ncol(XS) - 1)]
  # }

  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE, sep = "\t")
  dx_final <- dx[dx[,1] == max(dx[,1]), 2:ncol(dx)]
  dx_initial <- dx[dx[,1] == min(dx[,1]), 2:ncol(dx)]
  n_xs_each <- apply(dx_final, 1, function(x){
    sum(x > 0, na.rm = TRUE)
  })

  dx_final <- as.vector(t(dx_final))
  dx_final <- dx_final[!is.na(dx_final)]
  dx_final <- dx_final[dx_final > 0]

  dx_initial <- as.vector(t(dx_initial))
  dx_initial <- dx_initial[!is.na(dx_initial)]
  dx_initial <- dx_initial[dx_initial > 0]

  n_xs <- sum(XS[,1] == 0)
  area_diff <- numeric(length = n_xs)

  for (i in 1:n_xs){
    initial <- as.numeric(XS[i, 2:ncol(XS)])
    final <- as.numeric(XS[n_xs + i, 2:ncol(XS)])

    area_initial <- 0
    area_final <- 0
    for (k in 2:10){
      area_initial <- area_initial + (initial[k] - initial[k - 1]) *
        (initial[10 + k - 1] + initial[10 + k]) / 2
      # print((initial[k] - initial[k - 1]) *
      #         (initial[10 + k - 1] + initial[10 + k]) / 2)
      area_final <- area_final + (final[k] - final[k - 1]) *
        (final[10 + k - 1] + final[10 + k]) / 2
      # print((final[k] - final[k - 1]) *
      #   (final[10 + k - 1] + final[10 + k]) / 2)
    }

    area_diff[i] <- area_final - area_initial
    #(-) means erosion; (+) means deposition

  }
  volume <- area_diff * dx_final

  # volume <- numeric()
  # for (i in 1:n_xs){
  #   initial <- as.numeric(XS[i, 2:ncol(XS)])
  #   final <- as.numeric(XS[n_xs + i, 2:ncol(XS)])
  #
  #   area_initial <- 0
  #   area_final <- 0
  #   max_z_initial <- max(initial[11:20])
  #   max_z_final <- max(final[11:20])
  #   for (k in 2:10){
  #     area_initial <- area_initial + (initial[k] - initial[k - 1]) *
  #       (2 * max_z_initial - initial[10 + k - 1] - initial[10 + k]) / 2
  #     # print((initial[k] - initial[k - 1]) *
  #     #         (initial[10 + k - 1] + initial[10 + k]) / 2)
  #     area_final <- area_final + (final[k] - final[k - 1]) *
  #       (2 * max_z_final - final[10 + k - 1] - final[10 + k]) / 2
  #     # print((final[k] - final[k - 1]) *
  #     #   (final[10 + k - 1] + final[10 + k]) / 2)
  #   }
  #
  #   volume[i] <- area_initial * dx_initial[i] - area_final * dx_final[i]
  #   #(-) means erosion; (+) means deposition
  #
  # }

  #What are the indices of the first XS of each reach?
  indices <- cumsum(c(1, n_xs_each))
  # geom <- read.table("U:/PhD Work/Toutle River Data/Initial Channel Geometry.txt",
  #                    header = TRUE, sep = "\t")
  #
  par(mfrow = c(1,1), mar = c(4,4,0.5,0.5))
  plot(volume/1e6, type = "l", las = 1, xaxt = "n", ylab = "Volume [million m^3]",
       lwd = 2, xlab = "")
  points(volume/1e6)
  # axis(1, at = indices[1:19], labels = rev(geom$XS))
  # abline(v = indices[1:19])
  abline(h = 0, lty = 2)
  vol_sum <- sum(volume)

  #Get all outputs
  output_sed <- read.table(paste0(path, "/Output sediment vol.txt"),
                           header = TRUE, sep = "\t")

  bed_vol_out <- max(output_sed$Bed_vol_out)
  bed_vol_in <- max(output_sed$Bed_vol_in)
  bank_tank_final <- max(output_sed$Bank_tank_vol, na.rm = TRUE)
  bank_washload <- max(output_sed$Bank_washload)
  bed_washload <- max(output_sed$Cohesive_bed_washload)
  knick_washload <- max(output_sed$Knickpoint_washload)
  knick_correction <-output_sed$Knick_correction[nrow(output_sed)]

  diff <- round(vol_sum + bed_vol_out / 0.6 - bed_vol_in / 0.6 +
                  bank_tank_final + bank_washload + bed_washload +
                  knick_washload + knick_correction, 0)

  writeLines(paste("Volume sum:", round(vol_sum, 0), "\n",
                   "Bed vol out:", round(bed_vol_out / 0.6, 0), "\n",
                   "Bed vol in:", round(bed_vol_in / 0.6, 0), "\n",
                   "Bank tank:", round(bank_tank_final, 0), "\n",
                   "Bank washload:", round(bank_washload, 0), "\n",
                   "Bed washload (cohesive):", round(bed_washload, 0), "\n",
                   "Knickpoint washload:", round(knick_washload, 0), "\n",
                   "Knickpoint correction:", round(knick_correction, 0), "\n",
                   "Diff:", diff, "\n",
                   "Percent Diff:", round(diff / (vol_sum - bed_vol_in / 0.6)
                                          * 100, 2), "%"))

  #Caclulate integral (assume all dx's are one?)
  int <- 0
  for (i in 2:length(volume)){
    int <- int + 0.5 * (volume[i] + volume[i - 1])
  }

  #What are the indices of the first XS of each reach?
  # indices <- cumsum(c(1, n_xs_each))
  #
  # area_nodes <- area_diff[indices[1:19]]
  #
  # XS_num <- 20
  # plot(as.numeric(XS[XS_num, 2:11]), as.numeric(XS[XS_num, 12:21]))
  # lines(as.numeric(XS[XS_num, 2:11]), as.numeric(XS[XS_num, 12:21]))
  # points(as.numeric(XS[n_xs + XS_num, 2:11]), as.numeric(XS[n_xs + XS_num, 12:21]), col = "red")
  # lines(as.numeric(XS[n_xs + XS_num, 2:11]), as.numeric(XS[n_xs + XS_num, 12:21]), col = "red")

}

#' Creates a series of points visualizing the channel network
#'
#' @param n_nodes The number of channel reaches or nodes
#' @param n_xs Number of cross sections per reach
#' @param link Matrix specifying reach layout
#' @param dx Cross section spacing
#' @param custom_sgn Specifies the direction each reach should be plotted (`-1`
#'   is left, `1` is right)
#'
#' @return List of x and y coordinates of channel network
#'
make_network <- function(n_nodes, n_xs, link, dx, custom_sgn){

  angle <- pi/8
  x <- matrix(0, nrow = n_nodes + 1, ncol = max(n_xs))
  y <- x
  branched <- rep(1, n_nodes)

  if (length(custom_sgn) == 0){
    sgn <- rep(1, n_nodes)
  }else{
    sgn <- custom_sgn
  }

  for (i in seq(n_nodes, 1, -1)){
    if (i < n_nodes){
      index <- which(i == link, arr.ind = TRUE)[2]
      ang <- angle / branched[index]
      for (j in n_xs[i]:1){
        if (j == n_xs[i]){
          x[i, j] <- x[index, 1] + sgn[index] * sin(ang) * dx[i, j + 1]
          y[i, j] <- y[index, 1] + cos(ang) * dx[i, j + 1]
        }else{
          x[i, j] <- x[i, j + 1] + sgn[index] * sin(ang) * dx[i, j + 1]
          y[i, j] <- y[i, j + 1] + cos(ang) * dx[i, j + 1]
        }
      }

      sgn[index] <- -1 * sgn[index]
      branched[i] <- branched[index] + 1
    }else{
      index <- n_nodes + 1
      for (j in n_xs[i]:1){
        if (j == n_xs[i]){
          x[i, j] <- x[i + 1, 1]
          y[i, j] <- y[i + 1, 1] + dx[i, j + 1]
        }else{
          x[i, j] <- x[i, j + 1]
          y[i, j] <- y[i, j + 1] + dx[i, j + 1]
        }
      }
    }
  }

  return(list(x = x, y = y))
}

#' Plots changes in bed elevation for each cross section in the network
#'
#' @param print Should the plot be printed (defaults to `FALSE`)
#' @param gif Should a gif be created (defaults to `FALSE`). Must have
#'   ImageMagick installed and a folder titled "Figs" in the `path` directory.
#' @param max_plots Maximum number of plots in gif
#' @param path Path to folder with model outputs
#' @param custom_sgn Specifies the direction each reach should be plotted (`-1`
#'   is left, `1` is right)
#' @param title Title to be printed on plot
#'
#' @importFrom utils read.table
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot rect points legend mtext grconvertX grconvertY text
#'
#' @export
#'
dz_plot <- function(print = FALSE, gif = FALSE, max_plots = 10,
                    path = "", custom_sgn = NULL,
                    title = NULL){
  bed_z <- read.table(paste0(path, "/Output z.txt"), header = FALSE)
  link <- read.table(paste0(path, "/Input link.txt"), header = FALSE, sep = " ")
  #L <- read.table("Input length.txt", header = FALSE, sep = " ")

  #L <- as.array(L[,1])
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE)


  times <- unique(bed_z[,1])
  times <- times[-length(times)]
  n_nodes <- ncol(link)
  n_xs <- apply(bed_z[1:n_nodes,2:ncol(bed_z)], 1, function(x){sum(x > 0)})

  coords <- make_network(n_nodes, n_xs, link, dx, custom_sgn)
  x <- coords$x
  y <- coords$y

  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)

  if (gif){
    if (substr(path, 1, 1) == "~"){stop("Please supply full path, no ~ allowed.")}
    #Get times where we are going to make plots
    # max_round <- max(times[times %% 365 == 0])
    # factors <- (length(times) + 1:1000) / (2:1001)
    # n_plots <- max(factors[factors/floor(factors) == 1 & factors <= max_plots])
    #if (length(n_plots) == 0) {stop("Try increasing the max number of plots.")}
    times_plot <- seq(0, max(times), length.out = max_plots)

    #Get total bed elevation change
    dz_tot <- bed_z[bed_z[,1] == max(times),2:ncol(bed_z)] -
      bed_z[bed_z[,1] == 0,2:ncol(bed_z)]
    max_dz <- max(abs(dz_tot))

    legend_cols <- cRamp_legend(7, "RdBu")

    png(paste0(path, "/Figs/dz Network %02d.png"), type = "cairo", units = "in", height = 4.5, width = 4, res = 500)

    for (k in 1:length(times_plot)){
      #calculated bed elevation changes
      time_round <- times[which(abs(times - times_plot[k]) == min(abs(times - times_plot[k])))]
      dz <- bed_z[bed_z[,1] == time_round,2:ncol(bed_z)] -
        bed_z[bed_z[,1] == min(times_plot),2:ncol(bed_z)]

      #Convert dz to vector to get colors - add maximum dz to get right scale
      dz_vector <- c(as.vector(t(dz)), max_dz, -max_dz)
      colors <- cRamp(dz_vector, "RdBu")

      #Convert colors back to matrix - removing max dz value
      colors <- matrix(colors[1:(length(colors) - 2)], nrow = n_nodes, byrow = TRUE)

      par(mar = c(0.5, 0.5, 4, 0.5), mfrow = c(1,1))
      plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

      for (i in seq(n_nodes, 1, -1)){
        if (i < n_nodes){
          index <- which(i == link, arr.ind = TRUE)[2]
        }
        else{
          index <- n_nodes + 1
        }
        for (j in n_xs[i]:1){
          # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
          #        lwd = 1)
          points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 1.2)
        }
        #Add numbers
        #text(xy[i, 1], xy[i, 2], i, pos = 4)
      }

      legend("bottomright", legend = paste(round(time_round / 365, 1), "Years"),
             pch = NA, bty = "n", text.col = "white", cex = 1.3)

      xval <- grconvertX(seq(0.1, 0.9, length.out = 7), from = "nfc")
      y_point <- rep(grconvertY(0.92, from = "nfc"), length(xval))
      y_lab <- rep(grconvertY(0.9, from = "nfc"), length(xval))

      points(xval, y_point, pch = 21, bg = legend_cols, cex = 1.7, xpd = NA)
      text(xval, y_lab, sprintf("%.1f", round(seq(-max_dz, max_dz, length.out = 7), 1)),
           pos = 1, xpd = NA, cex = 1.2)
      mtext("Elevation Change [m]", line = 3, side = 3, cex = 1.2)
    }
    dev.off()
    system("cmd.exe", input = c(paste("cd", path),
                      "magick convert -delay 100 -loop 0 Figs/dz*.png Figs/dz_out.gif"))
  }else{
    if (print){
      png(paste0(path, "/Network plot.png"), type = "cairo", units = "in", height = 4.5, width = 4, res = 500)
    }

    #Get total bed elevation change
    dz <- bed_z[bed_z[,1] == max(times),2:ncol(bed_z)] -
      bed_z[bed_z[,1] == 0,2:ncol(bed_z)]
    max_dz <- max(abs(dz))
    legend_cols <- cRamp_legend(7, "RdBu")

    #Convert dz to vector to get colors
    dz_vector <- c(as.vector(t(dz)), max_dz, -max_dz)
    colors <- cRamp(dz_vector, "RdBu")

    #Convert colors back to matrix
    colors <- matrix(colors[1:(length(colors) - 2)], nrow = n_nodes, byrow = TRUE)

    par(mar = c(0.5, 0.5, 4, 0.5), mfrow = c(1,1))
    if (length(title) > 0){par(mar = c(0.5, 0.5, 1.5, 0.5))}
    plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

    if(length(title) > 0){mtext(side = 3, title, line = 0.2, cex = 1.5)}

    for (i in seq(n_nodes, 1, -1)){
      if (i < n_nodes){
        index <- which(i == link, arr.ind = TRUE)[2]
      }
      else{
        index <- n_nodes + 1
      }
      for (j in n_xs[i]:1){
        # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
        #        lwd = 1)
        points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 1.2)
      }
      #Add numbers
      #text(xy[i, 1], xy[i, 2], i, pos = 4)
    }

    # legend("topleft", legend = round(seq(-max_dz, max_dz, length.out = 7), 2),
    #        fill = legend_cols, bty = "n", title = "Elevation Change [m]",
    #        text.col = "white")
    legend("bottomright", legend = paste(round(max(times) / 365, 1), "Years"),
           pch = NA, bty = "n", text.col = "white", cex = 1.3)

    xval <- grconvertX(seq(0.1, 0.9, length.out = 7), from = "nfc")
    y_point <- rep(grconvertY(0.92, from = "nfc"), length(xval))
    y_lab <- rep(grconvertY(0.9, from = "nfc"), length(xval))

    points(xval, y_point, pch = 21, bg = legend_cols, cex = 1.7, xpd = NA)
    text(xval, y_lab, sprintf("%.1f", round(seq(-max_dz, max_dz, length.out = 7), 1)),
         pos = 1, xpd = NA, cex = 1.2)
    mtext("Elevation Change [m]", line = 3, side = 3, cex = 1.2)

    if (print){
      dev.off()
    }
  }
}

#' Plots changes in channel width for all cross sections in the network
#'
#' @param print Should the plot be printed (defaults to `FALSE`)
#' @param gif Should a gif be created (defaults to `FALSE`). Must have
#'   ImageMagick installed and a folder titled "Figs" in the `path` directory.
#' @param max_plots Maximum number of plots in gif
#' @param path Path to folder with model outputs
#' @param custom_sgn Specifies the direction each reach should be plotted (`-1`
#'   is left, `1` is right)
#' @param title Title to be printed on plot
#'
#' @importFrom utils read.table
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot rect points legend mtext grconvertX grconvertY
#'   text
#'
#' @export
#'
width_plot <- function(print = FALSE, gif = FALSE, max_plots = 10,
                       path = "", custom_sgn = NULL,
                       title = NULL){
  width <- read.table(paste0(path, "/Output width.txt"), header = FALSE)
  link <- read.table(paste0(path, "/Input link.txt"), header = FALSE, sep = " ")

  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE, sep = "\t")

  times <- unique(width[,1])
  times <- times[-length(times)]
  n_nodes <- ncol(link)
  n_xs <- apply(width[1:n_nodes,2:ncol(width)], 1, function(x){sum(x > 0)})

  coords <- make_network(n_nodes, n_xs, link, dx, custom_sgn)
  x <- coords$x
  y <- coords$y

  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)

  if (gif){
    if (substr(path, 1, 1) == "~"){stop("Please supply full path, no ~ allowed.")}

    #Get times where we are going to make plots
    times_plot <- seq(0, max(times), length.out = max_plots)

    #Get total width change
    dwidth_tot <- width[width[,1] == max(times),2:ncol(width)] -
      width[width[,1] == 0,2:ncol(width)]
    max_dwidth <- max(abs(dwidth_tot))
    legend_cols <- rev(cRamp_legend(7, "RdBu"))

    png(paste0(path, "/Figs/Width Network %02d.png"), type = "cairo", units = "in", height = 4.5, width = 4, res = 500)
    for (k in 1:length(times_plot)){
      #calculated width changes
      time_round <- times[which(abs(times - times_plot[k]) == min(abs(times - times_plot[k])))]
      dwidth <- width[width[,1] == time_round, 2:ncol(width)] -
        width[width[,1] == min(times_plot),2:ncol(width)]

      #Convert dz to vector to get colors - add maximum dwidth to get right scale
      width_vector <- -c(as.vector(t(dwidth)), max_dwidth, -max_dwidth)
      colors <- cRamp(width_vector, "RdBu")

      #Convert colors back to matrix - removing max dwidth value
      colors <- matrix(colors[1:(length(colors) - 2)], nrow = n_nodes, byrow = TRUE)

      par(mar = c(0.5, 0.5, 4, 0.5), mfrow = c(1,1))
      plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

      for (i in seq(n_nodes, 1, -1)){
        if (i < n_nodes){
          index <- which(i == link, arr.ind = TRUE)[2]
        }
        else{
          index <- n_nodes + 1
        }
        for (j in n_xs[i]:1){
          # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
          #        lwd = 1)
          points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 1.2)
        }
        #Add numbers
        #text(xy[i, 1], xy[i, 2], i, pos = 4)
      }

      legend("bottomright", legend = paste(round(time_round / 365, 1), "Years"),
             pch = NA, bty = "n", text.col = "white", cex = 1.3)

      xval <- grconvertX(seq(0.1, 0.9, length.out = 7), from = "nfc")
      y_point <- rep(grconvertY(0.92, from = "nfc"), length(xval))
      y_lab <- rep(grconvertY(0.9, from = "nfc"), length(xval))

      points(xval, y_point, pch = 21, bg = legend_cols, cex = 1.7, xpd = NA)
      text(xval, y_lab, sprintf("%.1f", round(seq(-max_dz, max_dz, length.out = 7), 1)),
           pos = 1, xpd = NA, cex = 1.2)
      mtext("Elevation Change [m]", line = 3, side = 3, cex = 1.2)
    }
    dev.off()
    system("cmd.exe", input = c(paste("cd", path),
                                "magick convert -delay 100 -loop 0 Figs/Width*.png Figs/width_out.gif"))
  }else{
    if (print){
      png(paste0(path, "/Width plot.png"), type = "cairo", units = "in", height = 4.5, width = 4, res = 500)
    }

    #Get total bed elevation change
    dwidth <- width[width[,1] == max(times),2:ncol(width)] -
      width[width[,1] == 0,2:ncol(width)]
    max_dwidth <- max(abs(dwidth))
    legend_cols <- rev(cRamp_legend(7, "RdBu"))

    #Convert dz to vector to get colors
    width_vector <- -c(as.vector(t(dwidth)), max_dwidth, -max_dwidth)
    if (max_dwidth == 0){
      colors <- rep(legend_cols[4], length(width_vector))
    }else{
      colors <- cRamp(width_vector, "RdBu")
    }

    #Convert colors back to matrix - removing max dwidth value
    colors <- matrix(colors[1:(length(colors) - 2)], nrow = n_nodes, byrow = TRUE)


    par(mar = c(0.5, 0.5, 4, 0.5), mfrow = c(1,1))
    if(length(title) > 0){par(mar = c(0.5, 0.5, 1.5, 0.5))}
    plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

    if(length(title) > 0){mtext(side = 3, title, line = 0.2, cex = 1.5)}
    for (i in seq(n_nodes, 1, -1)){
      if (i < n_nodes){
        index <- which(i == link, arr.ind = TRUE)[2]
      }
      else{
        index <- n_nodes + 1
      }
      for (j in n_xs[i]:1){
        # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
        #        lwd = 1)
        points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 1.2)
      }
      #Add numbers
      #text(xy[i, 1], xy[i, 2], i, pos = 4)
    }

    # legend("topleft", legend = round(seq(-max_dwidth, max_dwidth, length.out = 7), 2),
    #        fill = legend_cols, bty = "n", title = "Width Change [m]",
    #        text.col = "white")
    legend("bottomright", legend = paste(round(max(times) / 365, 1), "Years"),
           pch = NA, bty = "n", text.col = "white", cex = 1.3)

    xval <- grconvertX(seq(0.1, 0.9, length.out = 7), from = "nfc")
    y_point <- rep(grconvertY(0.92, from = "nfc"), length(xval))
    y_lab <- rep(grconvertY(0.9, from = "nfc"), length(xval))

    points(xval, y_point, pch = 21, bg = legend_cols, cex = 1.7, xpd = NA)
    text(xval, y_lab, sprintf("%.1f", round(seq(-max_dwidth, max_dwidth, length.out = 7), 1)),
         pos = 1, xpd = NA, cex = 1.2)
    mtext("Width Change [m]", line = 3, side = 3, cex = 1.2)


    if (print){
      dev.off()
    }
  }
}

#' Plots changes over time in bed elevation for the most upstream cross section
#' in each reach
#'
#' @param path Path to folder with model outputs
#' @param type `type = 1` plots all lines on same plot, `type = 2` creates a
#'   separate plot for each reach
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines text grconvertX grconvertY
#' @importFrom grDevices adjustcolor
#'
#' @export
#'
dz_lines <- function(path = "", type = 1){
  bed_z <- read.table(paste0(path, "/Output z.txt"), header = FALSE, sep = "\t")

  #remove last bed_z column and change 0 values to NA
  bed_z <- bed_z[,1:(ncol(bed_z) - 1)]
  bed_z[bed_z == 0] <- NA
  bed_z$V1[is.na(bed_z$V1)] <- 0

  times <- unique(bed_z$V1)
  n_nodes <- sum(bed_z$V1 == 0)

  if (type == 1){
    #calculated bed elevation changes
    dz <- matrix(0, nrow = length(times), ncol = n_nodes)
    initial <- bed_z$V2[bed_z$V1 == 0]
    for (i in 1:length(times)){
      dz[i,] <- as.numeric((bed_z$V2[bed_z$V1 == times[i]] - initial))
    }
    colors <- cRamp_legend(n_nodes, "viridis")

    par(mfrow = c(1,1), mar = c(3, 3, 0.5, 0.5))
    min_dz <- min(apply(dz, 2, min, na.rm = TRUE))
    max_dz <- max(apply(dz, 2, max, na.rm = TRUE))
    plot(NA, xlim = range(times), ylim = c(min_dz, max_dz), las = 1,
         ylab = "Elevation Change [m]", xlab = "Year")
    for (i in 1:ncol(dz)){
      lines(times, dz[,i], lty = i, col = colors[i])
      text(max(times), dz[nrow(dz), i], labels = i, pos = 4, xpd = NA)
    }
  }else{
    par(mfrow = c(2,1), mar = c(3, 3, 0.5, 0.5), oma = c(0, 0, 0, 3))
    for (j in 1:n_nodes){
      #calculated bed elevation changes
      n_xs <- sum(!is.na(bed_z[j,2:ncol(bed_z)]))
      dz <- matrix(0, nrow = length(times), ncol = ncol(bed_z) - 1)
      initial <- bed_z[j, 2:ncol(bed_z)]
      for (i in 1:length(times)){
        dz[i,] <- as.numeric((bed_z[j + n_nodes * (i - 1), 2:ncol(bed_z)] - initial))
      }
      colors <- cRamp_legend(n_xs, "viridis")
      colors <- adjustcolor(colors, alpha.f = 0.7)

      min_dz <- min(apply(dz[,1:n_xs], 2, min, na.rm = TRUE))
      max_dz <- max(apply(dz[,1:n_xs], 2, max, na.rm = TRUE))
      plot(NA, xlim = c(0, nrow(dz) + 1), ylim = c(min_dz, max_dz), las = 1,
           ylab = "Elevation Change [m]", xlab = "Day", main = paste("Reach", j))
      for (i in 1:ncol(dz)){
        lines(1:nrow(dz), dz[,i], col = colors[i], lwd = 2)
        #text(nrow(dz), dz[nrow(dz), i], labels = i, pos = 4, xpd = NA)
      }

      if (i %% 2 != 0){
        par(xpd = NA)
        xvals <- grconvertX(c(1.0, 1.05), from = "nfc", to = "user")
        yvals <- grconvertY(c(0.4, 0.6), from = "ndc", to = "user")
        color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                        "","DS"),
                     align = "rb", gradient = "y",
                     rect.col = cRamp_legend(5, "viridis"), xpd = NA)
      }

    }

  }
}

#' Plots initial and final channel bed profile for each reach
#'
#' @param path Path to folder with model outputs
#' @param type `type = 1` plots all lines on same plot, `type = 2` creates a
#'   separate plot for each reach
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines points text
#' @importFrom grDevices rgb
#'
#' @export
#'
profiles <- function(path = "", type = 1){
  bed_z <- read.table(paste0(path, "/Output z.txt"), header = FALSE)

  #Change 0 values to NA
  bed_z[bed_z == 0] <- NA
  bed_z$V1[is.na(bed_z$V1)] <- 0

  times <- unique(bed_z$V1)
  n_nodes <- sum(bed_z$V1 == 0)

  #Get dx
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE)
  max_dx <- max(apply(dx[, 2:ncol(dx)], 2, max, na.rm = TRUE))
  xs <- 0:(ncol(bed_z) - 2)
  x_max <- max_dx * max(xs)

  min_z <- min(apply(bed_z[, 2:ncol(bed_z)], 2, min, na.rm = TRUE))
  max_z <- max(apply(bed_z[, 2:ncol(bed_z)], 2, max, na.rm = TRUE))

  if (type == 1){
    par(mar = c(4, 4, 0.5, 0.5), mfrow = c(1,1))
    plot(NA, xlim = c(0, x_max), ylim = c(min_z, max_z), las = 1,
         ylab = "Elevation [m]", xlab = "Station")
    for (i in 1:n_nodes){
      x <- cumsum(c(1, as.numeric(dx[i, 2:ncol(dx)])))
      x <- x[1:(length(x) - 1)]

      #initial profile
      lines(x, bed_z[i, 2:ncol(bed_z)])
      points(x, bed_z[i, 2:ncol(bed_z)], pch = 21, bg = rgb(0, 0, 0, 0.5))
      #final profile
      x <- cumsum(c(1, as.numeric(dx[nrow(dx) - n_nodes + i, 2:ncol(dx)])))
      x <- x[1:(length(x) - 1)]
      lines(x, bed_z[nrow(bed_z) - n_nodes + i, 2:ncol(bed_z)],
            col = "red")
      points(x, bed_z[nrow(bed_z) - n_nodes + i, 2:ncol(bed_z)],
             pch = 21, bg = rgb(1, 0, 0, 0.5))

      text(0.7, bed_z$V2[i], labels = i, pos = 2, xpd = NA)
    }
  }else{
    #n_plots <- ceiling(n_nodes / 6)
    par(mfrow = c(1, 2), mar = c(2, 2, 1, 0.5), oma = c(2, 2, 0, 0), mgp = c(2, 0.7, 0))
    for (i in 1:n_nodes){
      x <- cumsum(c(1, as.numeric(dx[i, 2:ncol(dx)])))
      x <- x[1:(length(x) - 1)]

      xlim <- c(0, max(x))
      ylim <- range(c(range(bed_z[i, 2:ncol(bed_z)], na.rm = TRUE),
                      range(bed_z[nrow(bed_z) - n_nodes + i, 2:ncol(bed_z)],
                            na.rm = TRUE)))
      plot(NA, xlim = xlim, ylim = ylim, las = 1,
           ylab = "", xlab = "", main = paste("Reach", i))

      #initial profile
      lines(x, bed_z[i, 2:ncol(bed_z)])
      points(x, bed_z[i, 2:ncol(bed_z)], pch = 21, bg = rgb(0, 0, 0, 0.5))
      #final profile
      x <- cumsum(c(1, as.numeric(dx[nrow(dx) - n_nodes + i, 2:ncol(dx)])))
      x <- x[1:(length(x) - 1)]
      lines(x, bed_z[nrow(bed_z) - n_nodes + i, 2:ncol(bed_z)],
            col = "red")
      points(x, bed_z[nrow(bed_z) - n_nodes + i, 2:ncol(bed_z)],
             pch = 21, bg = rgb(1, 0, 0, 0.5))

      if (i %% 2 == 0){
        mtext("Distance [m]", side = 1, outer = TRUE, line = 0.5)
        mtext("Elevation [m]", side = 2, outer = TRUE, line = 0.5)
      }
      #text(0.7, bed_z$V2[i], labels = i, pos = 2, xpd = NA)
    }
  }


}

#' Plots changes in channel width over time for all cross sections
#'
#' @param path Path to folder with model outputs
#' @param print Should the plot be printed to a file (defaults to `FALSE`)
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines mtext grconvertX grconvertY
#'
#' @export
#'
width_lines <- function(path = "", print = FALSE){
  width <- read.table(paste0(path, "/Output width.txt"), header = FALSE)

  times <- unique(width$V1)
  n_nodes <- sum(width$V1 == 0)

  #Change zeroes to NA
  width[width == 0] <- NA

  #Get dx
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE)
  xs <- 0:(ncol(width) - 1)
  x_max <- max(dx) * max(xs)

  min_width <- min(apply(width[, 2:ncol(width)], 2, min, na.rm = TRUE))
  max_width <- max(apply(width[, 2:ncol(width)], 2, max, na.rm = TRUE))

  if (n_nodes > 1) {
    n_col <- 2
    n_row <- 2
  } else{
    n_row <- 1
    n_col <- 1
  }

  par(mfrow = c(n_row, n_col), mar = c(2, 2, 1, 0.5), oma = c(2, 2, 3, 0))
  for (i in 1:n_nodes){
    row_seq <- seq(i, nrow(width), n_nodes)
    sub <- width[row_seq, ]

    n_colors <- sum(!is.na(sub[1, ])) - 1
    max_width <- max(apply(sub[, 2:(n_colors + 1)], 2, max, na.rm = TRUE))

    plot(NA, xlim = c(0, max(times)), ylim = c(0, max_width), las = 1,
         main = paste("Reach", i), xlab = "", ylab = "")
    colors <- cRamp_legend(n_colors, "viridis")
    for (j in 2:(n_colors + 1)){
      lines(times, sub[,j], lwd = 1, col = colors[j - 1])
      #points(times, sub[,j], pch = 16)
    }

    if (i %in% seq(from = 1, length.out = 10, by = 4)){
      mtext("Days", side = 1, line = 1, outer = TRUE)
      mtext("Width [m]", side = 2, line = 1, outer=  TRUE)

      xvals <- grconvertX(x = c(0.4, 0.6), from = "ndc", to = "user")
      yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
      color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                      "","DS"),
                   align = "rb", gradient = "x",
                   rect.col = cRamp_legend(5, "viridis"), xpd = NA)
    }
  }

}

#' Plots changes in channel slope over time for all cross sections
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines mtext grconvertX grconvertY
#'
#' @export
#'
slope_lines <- function(path = ""){
  slope <- read.table(paste0(path, "/Output slope.txt"), header = FALSE)
  times <- unique(slope$V1)
  slope_XS <- data_by_XS(slope)
  n_XS <- length(slope_XS)

  n_row <- min(c(n_XS, 4))
  par(mfrow = c(n_row, 1), mar = c(2, 3, 1.5, 0.5), oma = c(4, 3, 3, 0), mgp = c(2, 0.5, 0))
  for (i in 1:n_XS){
    n_reach <- sum(colMeans(slope_XS[[i]]) > 0) - 1
    colors <- cRamp_legend(n_reach, "viridis")
    ylim <- range(slope_XS[[i]][, 2:(n_reach + 1)])
    plot(NA, ylab = "", xlab = "", las = 1, ylim = ylim, xlim = range(times),
         main = paste("Reach", i))
    for (j in 1:n_reach){
      lines(times, slope_XS[[i]][,j + 1], col = colors[j])
    }

    if (i %in% seq(from = 1, length.out = 10, by = 4)){
      mtext("Days", side = 1, line = 2, outer = TRUE)
      mtext("Slope", side = 2, line = 1.5, outer=  TRUE)

      xvals <- grconvertX(x = c(0.4, 0.6), from = "ndc", to = "user")
      yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
      color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                      "","DS"),
                   align = "rb", gradient = "x",
                   rect.col = cRamp_legend(5, "viridis"), xpd = NA)
    }
  }

}

#' Plots initial and final cross section geometry for the most upstream cross
#' section in each reach
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines points mtext
#' @importFrom grDevices rgb
#'
#' @export
#'
XS_plots <- function(path = "") {

  output <- read.table(paste0(path, "/Output XS geometry.txt"))
  n_XS <- sum(output[,1] == 0)
  par(mfrow = c(1, 2), mar = c(2.5, 3, 1, 0.5), oma = c(1, 1.1, 0, 0))
  #par(mfrow = c(3, 2), mar = c(3,3,1,1))
  for (i in 1:n_XS){
    #file_nm <- paste0("U:/PhD Work/Toutle River Data/XS Plots/", xs_names[i], ".png")
    #png(file_nm, type = "cairo", units = "in", height = 4, width = 5, res = 500)
    #colors <- colorRampPalette(c("green", "purple"))(ncol)

    #model geometry
    modx1 <- as.numeric(output[i,3:10])
    modz1 <- as.numeric(output[i,13:20])
    modx2 <- as.numeric(output[nrow(output) - n_XS + i, 3:10])
    modz2 <- as.numeric(output[nrow(output) - n_XS + i, 13:20])

    xmin <- min(c(modx1, modx2), na.rm = TRUE)
    xmax <- max(c(modx1, modx2), na.rm = TRUE)
    ymin <- min(c(modz1, modz2), na.rm = TRUE)
    ymax <- max(c(modz1, modz2), na.rm = TRUE)

    plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
         xlab = "", ylab = "", las = 1, main = paste("Reach", i))
    # for (j in 1:ncol){
    #   lines(elev_m ~ dist_m, test[test$date == unique(test$date)[j], ], col = colors[j])
    # }

    #add simplified geometry plot
    #initial
    lines(modx1, modz1)
    points(modx1, modz1, pch = 16, col = rgb(0,0,0,0.7))
    #final
    lines(modx2, modz2, col = "red")
    points(modx2, modz2, pch = 16, col = rgb(1,0,0,0.7))

  }
  mtext("Station [m]", side = 1, line = 0, outer = TRUE)
  mtext("Elevation [m]", side = 2, line = 0, outer = TRUE)
}

#' Plots initial and final cross section geometry for any specified cross sections
#'
#' @param path Path to folder with model outputs
#' @param reach A numeric vector of reaches for which to plot cross sections
#' @param XS A numeric vector of the cross sections within each reach to be
#' plotted (one `XS` per `reach`)
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines points mtext
#' @importFrom grDevices rgb adjustcolor
#' @importFrom dplyr %>%
#'
#' @export
XS_plots2 <- function(path = "", reach = 1, XS = 1){

  output <- read.table(paste0(path, "/Output XS geometry all.txt")) %>%
    dplyr::filter_(~ V1 == 0 | V1 == max(V1))
  n_XS <- sum(output[,1] == 0)
  par(mfrow = c(2, 1), mar = c(2.5, 3, 1, 0.5), oma = c(1, 1.1, 0, 0))
  #par(mfrow = c(3, 2), mar = c(3,3,1,1))

  #Get all XS indices based on reach-XS pair
  #need to get #XS by reach
  dz_output <- read.table(file.path(path, "Output z.txt")) %>%
    dplyr::filter_(~ V1 == 0) %>%
    dplyr::select_(~ -V1)

  n_XS_reach <- apply(dz_output, 1, function(x){sum(x != 0)})

  XS_all <- purrr::map2_dbl(reach, XS, function(x, y, n_XS_reach){
    if (x == 1 & !(y > n_XS_reach[x])){
      XS_ID = y
    }else if(!(y > n_XS_reach[x])){
      XS_ID = sum(n_XS_reach[1:(x - 1)]) + y
    }else{
      XS_ID = NA
    }

    return(XS_ID)
  }, n_XS_reach)

  #Drop NAs from all vectors
  reach <- reach[!is.na(XS_all)]
  XS <- XS[!is.na(XS_all)]
  XS_all <- XS_all[!is.na(XS_all)]

  count <- 1
  for (i in XS_all){
    #file_nm <- paste0("U:/PhD Work/Toutle River Data/XS Plots/", xs_names[i], ".png")
    #png(file_nm, type = "cairo", units = "in", height = 4, width = 5, res = 500)
    #colors <- colorRampPalette(c("green", "purple"))(ncol)

    #model geometry
    modx1 <- as.numeric(output[i,3:10])
    modz1 <- as.numeric(output[i,13:20])
    modx2 <- as.numeric(output[n_XS + i, 3:10])
    modz2 <- as.numeric(output[n_XS + i, 13:20])

    xmin <- min(c(modx1, modx2), na.rm = TRUE)
    xmax <- max(c(modx1, modx2), na.rm = TRUE)
    ymin <- min(c(modz1, modz2), na.rm = TRUE)
    ymax <- max(c(modz1, modz2), na.rm = TRUE)

    plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
         xlab = "", ylab = "", las = 1, main = paste("Reach", reach[count], "XS", XS[count]))
    # for (j in 1:ncol){
    #   lines(elev_m ~ dist_m, test[test$date == unique(test$date)[j], ], col = colors[j])
    # }

    #add simplified geometry plot
    #initial
    lines(modx1, modz1, col = rgb(0,0,0,0.7), lwd = 2)
    points(modx1, modz1, col = rgb(0,0,0,0.7), pch = 16)
    #final
    lines(modx2, modz2, col = adjustcolor("red", 0.7), lwd = 2)
    points(modx2, modz2, col = adjustcolor("red", 0.7), pch = 16)

    count <- count + 1

  }
  mtext("Station [m]", side = 1, line = 0, outer = TRUE)
  mtext("Elevation [m]", side = 2, line = 0, outer = TRUE)
}

#' Plots cross section geometry over time for any specified cross section
#'
#' @param path Path to folder with model outputs
#' @param XS The number of the cross section to be plotted
#' @param n_plots The number of cross sections to plot
#' @param ts The time step (in seconds) between plottings
#' @param print Should the plot be printed to a file (defaults to `FALSE`)
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines points grconvertX grconvertY
#' @importFrom grDevices rgb png dev.off
#'
#' @export
XS_plots3 <- function(path = "", XS = 1, n_plots = 0, ts = 0.2, print = FALSE){
  output <- read.table(paste0(path, "/Output XS geometry.txt"))
  n_XS <- sum(output[,1] == 0)
  times <- unique(output[,1])
  par(mfrow = c(1, 1), mar = c(2.5, 3, 1, 0.5), oma = c(1, 1.1, 0, 0))
  #par(mfrow = c(3, 2), mar = c(3,3,1,1))

  if (n_plots == 0) {
    #model geometry
    modx1 <- as.numeric(output[XS,2:11])
    modz1 <- as.numeric(output[XS,12:21])
    modx2 <- as.numeric(output[nrow(output) - n_XS + XS, 2:11])
    modz2 <- as.numeric(output[nrow(output) - n_XS + XS, 12:21])

    xmin <- min(c(modx1, modx2), na.rm = TRUE)
    xmax <- max(c(modx1, modx2), na.rm = TRUE)
    ymin <- min(c(modz1, modz2), na.rm = TRUE)
    ymax <- max(c(modz1, modz2), na.rm = TRUE)

    plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
         xlab = "", ylab = "", las = 1)
    for (i in times) {
      row1 <- which(output$V1 == i)[XS]
      x <- as.numeric(output[row1, 2:11])
      z <- as.numeric(output[row1, 12:21])
      #initial
      lines(x, z, col = rgb(0, 0, 0, 0.3))
      points(x, z, col = rgb(0, 0, 0, 0.3))

      Sys.sleep(ts)

    }
  }else{
    #model geometry
    modx1 <- as.numeric(output[XS,3:10])
    modz1 <- as.numeric(output[XS,13:20])
    modx2 <- as.numeric(output[nrow(output) - n_XS + XS, 3:10])
    modz2 <- as.numeric(output[nrow(output) - n_XS + XS, 13:20])

    xmin <- min(c(modx1, modx2), na.rm = TRUE)
    xmax <- max(c(modx1, modx2), na.rm = TRUE)
    ymin <- min(c(modz1, modz2), na.rm = TRUE)
    ymax <- max(c(modz1, modz2), na.rm = TRUE)

    if (print){
      png(paste0(path, "/XS Evolution.png"), type = "cairo", units = "in",
          height = 3, width = 6, res = 500)
    }
    par(mfrow = c(1,1), mar = c(4, 4, 1, 0.5), oma = c(0, 0, 2, 0))

    plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
         xlab = "Station [m]", ylab = "Elevation[m]", las = 1, main = paste("XS", XS))

    times2 <- times[floor(seq(1, length(times), length.out = n_plots))]
    colors <- cRamp_legend(length(times2), "viridis")
    count <- 1
    for (i in times2) {
      row1 <- which(output$V1 == i)[XS]
      x <- as.numeric(output[row1, 3:10])
      z <- as.numeric(output[row1, 13:20])
      #initial
      lines(x, z, col = colors[count], lwd = 2)

      count <- count + 1

    }

    xvals <- grconvertX(x = c(0.4, 0.6), from = "npc", to = "user")
    yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
    color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c(expression(t[o]), "","",
                                                                    "",expression(t[f])),
                 align = "rb", gradient = "x",
                 rect.col = cRamp_legend(5, "viridis"), xpd = NA)

    if (print){
      dev.off()
    }
  }

  # rows <- seq(from = XS, by = n_XS, length.out = length(times))
  # y <- output$V20[rows]
  # plot(times, y, type = "l")
  # points(times, y)
  # mtext("Station [m]", side = 1, line = 0, outer = TRUE)
  # mtext("Elevation [m]", side = 2, line = 0, outer = TRUE)
}

#' Plots changes in channel sinuosity over time by reach
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom utils read.table
#' @importFrom grDevices adjustcolor
#' @importFrom graphics par plot lines grconvertX grconvertY
#'
#' @export
#'
sinuosity_plot <- function(path = ""){
  sinuosity <- read.table(paste0(path, "/Output sinuosity.txt"), header = FALSE)

  max_sinuosity <- max(apply(sinuosity,2,max,na.rm = TRUE))
  x <- 1:nrow(sinuosity)

  n_col <- ncol(sinuosity)
  colors <- cRamp_legend(max(n_col, 3), "viridis")
  colors <- adjustcolor(colors, alpha.f = 0.8)

  par(mar = c(4, 4, 3, 0.5))
  plot(NA, ylim = c(1, max_sinuosity), xlim = c(1, max(x)), las = 1, xlab = "Time",
       ylab = "Sinuosity")
  for (i in 1:ncol(sinuosity)){
    lines(x, sinuosity[ ,i], col = colors[i], lwd = 2)
  }

  xvals <- grconvertX(x = c(0.4, 0.6), from = "ndc", to = "user")
  yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
  color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                  "","DS"),
               align = "rb", gradient = "x",
               rect.col = cRamp_legend(5, "viridis"), xpd = NA)

  sinuosity_diff <- sinuosity[nrow(sinuosity), ] - sinuosity[1, ]

  return(sinuosity_diff)
}

#' Plots changes in bend radius of curvature over time for all cross sections
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom graphics par plot lines mtext grconvertX grconvertY
#' @importFrom utils read.table
#'
#' @export
#'
Rc_lines <- function(path = ""){
  Rc <- read.table(paste0(path, "/Output Rc.txt"), header = FALSE)
  times <- unique(Rc$V1)
  Rc[Rc == 0] <- NA
  Rc$V1[is.na(Rc$V1)] <- 0

  Rc <- data_by_XS(Rc)
  n_reaches <- length(Rc)

  par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.5), oma = c(2, 2, 0, 3))
  for (i in 1:n_reaches){
    ylim <- range(Rc[[i]][,2:ncol(Rc[[i]])], na.rm = TRUE)
    n_col <- sum(!is.na(apply(Rc[[i]][,2:ncol(Rc[[i]])], 2, sum)))
    colors <- cRamp_legend(n_col, "viridis")
    plot(NA, xlim = range(times) / 365, ylim = ylim, las = 1, ylab = "",
         xlab = "", main = paste("Reach", i))
    for(j in 1:n_col){
      lines(times / 365, Rc[[i]][,1+j], lwd = 2, col = colors[j])
    }

    if (i %in% seq(from = 1, by = 4, length.out = 10)){
      mtext("Rc [m]", side = 2, line = 0, outer = TRUE)
      mtext("Years", side = 1, line = 0, outer = TRUE)

      xvals <- grconvertX(x = c(0.95, 0.98), from = "ndc", to = "user")
      yvals <- grconvertY(y = c(0.4, 0.6), from = "ndc", to = "user")
      color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                      "","DS"),
                   align = "lt", gradient = "y",
                   rect.col = cRamp_legend(5, "viridis"), xpd = NA)

    }
  }


}

#' Plots changes in bed median grain size over time for all cross sections
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics par plot lines
#' @importFrom utils read.table
#'
#' @export
#'
D50_plot <- function(path = ""){
  D50 <- read.table(paste0(path, "/Output D50.txt"), header = FALSE, sep = "\t")

  #remove last column and change 0 values to NA
  D50 <- D50[,1:(ncol(D50) - 1)]
  D50[D50 == 0] <- NA
  D50$V1[is.na(D50$V1)] <- 0

  times <- unique(D50$V1) / 365
  n_nodes <- sum(D50$V1 == 0)

  #Get dx
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE, sep = "\t")
  xs <- 0:(ncol(D50) - 2)
  x_max <- max(dx, na.rm = TRUE) * max(xs)

  min_D50 <- min(apply(D50[, 2:ncol(D50)], 2, min, na.rm = TRUE))
  max_D50 <- max(apply(D50[, 2:ncol(D50)], 2, max, na.rm = TRUE))

  # if (n_nodes > 1) {
  #   n_col <- 2
  #   n_row <- min(ceiling(n_nodes / 2), 3)
  # } else{
  #   n_row <- 1
  #   n_col <- 1
  # }
  n_row <- 2
  n_col <- 1

  par(mfrow = c(n_row, n_col), mar = c(3, 3.5, 1, 0.5), mgp = c(2, 1, 0))
  for (i in 1:n_nodes){
    row_seq <- seq(i, nrow(D50), n_nodes)
    sub <- D50[row_seq, ]
    n_colors <- sum(!is.na(sub[1, ])) - 1

    max_D50<- max(apply(sub[, 2:(n_colors + 1)], 2, max, na.rm = TRUE))

    colors <- cRamp_legend(n_colors, "viridis")
    colors <- adjustcolor(colors, alpha.f = 0.8)

    plot(NA, xlim = c(0, max(times)), ylim = c(0, max_D50), las = 1,
         main = paste("Reach", i), ylab = expression(D[50]~"[mm]"),
         xlab = "Years")
    for (j in 2:(n_colors + 1)){
      lines(times, sub[,j], lwd = 2, col = colors[j - 1])
      #points(times, sub[,j], pch = 16, col = colors[j])
    }

    if (i %% 2 != 0){
      xvals <- max(times) * c(0.1, 0.4)
      yvals <- max_D50 * c(0.15, 0.2)
      plotrix::color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                               "","DS"),
                            align = "rb", gradient = "x",
                            rect.col = cRamp_legend(5, "viridis"))
    }
  }
}


# bank_loading <- function(path = "", start_year = 0,
#                          n_MC = 0,
#                          MC_path = ""){
#   bank_loads <- read.table(paste0(path, "/Output bank loading.txt"), header = TRUE)
#
#   #find number of days
#   days <- 1:nrow(bank_loads)
#   n_days <- length(days)
#
#   p_loads <- bank_loads[,which(substr(colnames(bank_loads), 1, 1) == "P")]
#   sed_loads <- bank_loads[,which(substr(colnames(bank_loads), 1, 3) == "Sed")]
#
#   years <- ceiling(days / 365) + start_year
#
#   p_loads <- apply(p_loads, 1, sum)
#   sed_loads <- apply(sed_loads, 1, sum)
#
#   combined <- data.frame(p_loads, sed_loads, years)
#
#   combined <- combined %>%
#     group_by(years) %>%
#     summarize(annual_P = sum(p_loads),
#               annual_sed = sum(sed_loads),
#               n = n()) %>%
#     filter(n > 300)
#
#   # stats <- combined %>%
#   #   summarize(mean_P = mean(annual_P),
#   #             mean_sed = mean(annual_sed))
#
#   if (n_MC > 0){
#     MC_loading <- list()
#     for (i in 0:(n_MC - 1)){
#       MC_loading[[i + 1]] <- read.table(paste0(MC_path,
#                                                "/MC Outputs/Output bank loading", i, ".txt"),
#                                         sep = "\t", header = TRUE)
#     }
#     MC_P <- lapply(MC_loading, function(x, years){
#       sub <- x[,which(substr(colnames(x), 1, 1) == "P")]
#       p_day <- apply(sub, 1, sum)
#       comb <- data.frame(p_day, years) %>%
#         group_by(years) %>%
#         summarise(p_load = sum(p_day),
#                   n_days = n()) %>%
#         filter(n_days > 300)
#       return(comb$p_load)
#     }, years)
#     MC_P <- do.call("rbind", MC_P)
#     MC_sed <- lapply(MC_loading, function(x, years){
#       sub <- x[,which(substr(colnames(x), 1, 3) == "Sed")]
#       sed_day <- apply(sub, 1, sum)
#       comb <- data.frame(sed_day, years) %>%
#         group_by(years) %>%
#         summarise(sed_load = sum(sed_day),
#                   n_days = n()) %>%
#         filter(n_days > 300)
#       return(comb$sed_load)
#     }, years)
#     MC_sed <- do.call("rbind", MC_sed)
#
#     P_med <- apply(MC_P, 2, median)
#     P_05 <- apply(MC_P, 2, quantile, probs = 0.05)
#     P_95 <- apply(MC_P, 2, quantile, probs = 0.95)
#
#     sed_med <- apply(MC_sed, 2, median)
#     sed_05 <- apply(MC_sed, 2, quantile, probs = 0.05)
#     sed_95 <- apply(MC_sed, 2, quantile, probs = 0.95)
#
#     # MC_stats <- data.frame(mean_pmed = mean(P_med),
#     #                        mean_p05 = mean(P_05),
#     #                        mean_p95 = mean(P_95),
#     #                        mean_sed = mean(sed_med),
#     #                        mean_sed05 = mean(sed_05),
#     #                        mean_sed95 = mean(sed_95))
#     #
#     # stats <- cbind(stats, MC_stats)
#   }
#
#   colors <- RColorBrewer::brewer.pal(5, "PuOr")
#   png(paste0(path, "/Loading plot.png"), type = "cairo", units = "in",
#       height = 4.5, width = 6.5, res = 500)
#   par(mfrow = c(2, 1), mar = c(2, 4.5, 0, 0.5), oma = c(1.5, 0, 1.5, 0),
#       mgp = c(2.5, 0.7, 0))
#   ylim <- c(0, max(combined$annual_P))
#   if (n_MC > 0){
#     ylim <- c(0, max(P_95))
#   }
#   plot(annual_P ~ years, combined, type = "n", ylab = "Annual P\nLoading [kg]",
#        las = 1, xlab = "", xpd = NA,
#        ylim = ylim)
#   if (n_MC > 0){
#     polygon(c(combined$years, rev(combined$years)), c(P_05, rev(P_95)),
#             col = colors[2], border = NA)
#     lines(combined$years, P_med, col = colors[1], lty = 2, lwd = 2)
#   }
#   lines(annual_P ~ years, combined, col = colors[1], lwd = 2)
#   if (n_MC > 0){
#     legend("top", legend = c("Single Run", "MC Median", "90% CI"),
#            pch = c(NA, NA, 15), pt.cex = 1.5,
#            lwd = c(2, 2, NA), lty = c(1, 2, NA),
#            col = c(rep(colors[1], 2), colors[2]),
#            bty = "n", horiz = TRUE, inset = c(0, -0.23), xpd = NA)
#   }
#
#   ylim <- c(0, max(combined$annual_sed / 1000))
#   if (n_MC > 0){
#     ylim <- c(0, max(sed_95 / 1000))
#   }
#   plot(I(annual_sed / 1000) ~ years, combined, type = "n",
#        ylab = "Annual Sed\nLoading [tons]", xlab = "Year", xpd = NA,
#        las = 1, ylim = ylim)
#   if (n_MC > 0){
#     polygon(c(combined$years, rev(combined$years)), c(sed_05 / 1000,
#                                                       rev(sed_95) / 1000),
#             col = colors[4], border = NA)
#     lines(combined$years, sed_med / 1000, col = colors[5], lty = 2, lwd = 2)
#   }
#   lines(I(annual_sed / 1000) ~ years, combined, lwd = 2, col = colors[5])
#   dev.off()
#
#   if(n_MC == 0){
#     return(combined)
#   }else{
#     return(list(modeled = combined, MC_P = MC_P, MC_sed = MC_sed))
#   }
# }
#


#' Plots knickpoint location over time
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot mtext legend points
#' @importFrom dplyr %>%
#'
#' @export
#'
knickpoint_plot <- function(path = ""){

  knick_x <- read.table(paste0(path, "/Output knick locations.txt"),
                        header = FALSE, sep = "\t")
  knick_x <- knick_x[,1:(ncol(knick_x) - 1)]

  n_xs <- sum(knick_x[,1] == 0)

  #Change zeros to NA
  knick_x[knick_x == 0] <- NA
  #Change initial times back to zero
  knick_x[1:n_xs, 1] <- 0

  reorganize <- function(start, x, n_row){
    rows <- seq(start, nrow(x) - n_row + start, n_row)

    x_out <- x[rows,]
    return(x_out)
  }

  plot_points <- function(y, color, x){
    points(x, -y, pch = 21, bg = color)
  }

  by_reach <- lapply(1:n_xs, reorganize, knick_x, n_xs)
  names(by_reach) <- 1:length(by_reach)

  par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.5), oma = c(1, 2, 0, 0))
  for (i in 1:length(by_reach)){
    x <- by_reach[[i]]
    #find columns with data
    columns <- colSums(x[,2:ncol(x)], na.rm = TRUE) > 0
    plot_num <- 0
    if (sum(columns) > 0){
      plot_num <- plot_num + 1
      x_val <- x[,1] / 365
      ylim <- -rev(range(x[,2:ncol(x)], na.rm = TRUE))
      n_col <- sum(columns)
      colors <- cRamp_legend(n_col, "viridis")
      colors2 <- rep(NA, ncol(x) - 1)
      colors2[columns] <- colors
      plot(NA, xlim = range(x_val), ylim = ylim, las = 1, xlab = "",
           ylab = "", main = paste("Reach", i))
      if (plot_num %in% seq(from = 1, by = 4, length.out = 100)){
        mtext("Time [years]", side = 1, line = 0, outer = TRUE)
        mtext("Knick Location [m]", side = 2, line = 0, outer = TRUE)
      }
      out <- mapply(plot_points, as.data.frame(x[,2:ncol(x)]), colors2,
                    MoreArgs = list(x = x_val))

      x2 <- reshape2::melt(x, id.vars = "V1", value.name = "Knick_x") %>%
        dplyr::filter_(~ !is.na(Knick_x)) %>%
        dplyr::arrange_(~ V1) %>%
        dplyr::rename_(.dots = list(Days = "V1"))

      migration_rate <- -round(stats::lm(Knick_x ~ I(Days/365), x2)$coefficients[2], 1)
      legend("bottomright", legend = bquote(.(migration_rate) ~ "m/yr"),
             pch = NA, bty = "n")
    }

  }

}


#' Calculates average of top and bottom channel width
#'
#' @param path Path to folder with model outputs
#' @param plot Should channel width be plotted (defaults to `FALSE`)
#'
#' @return A list of average channel widths by reach
#'
#' @importFrom utils read.table
#' @importFrom graphics lines
#'
avg_widths <- function(path = "", plot = FALSE){

  XS_output <- read.table(paste0(path, "/Output XS geometry.txt"))
  n_XS <- sum(XS_output[,1] == 0)
  times <- unique(XS_output[,1])

  XS <- lapply(1:n_XS, function(x, data, n_XS){
    data[seq(x, nrow(data) - n_XS + x, n_XS), ]
  }, XS_output, n_XS)

  #average of top and bottom width
  TW <- lapply(XS, function(data){
    apply(data, 1, function(x){
      mean(c(x[9] - x[4], x[7] - x[6]))
    })
  })


  if (plot){
    ylim <- range(lapply(TW, range))

    plot(NA, xlim = range(times/365), ylim = ylim, las = 1, xlab = "Years",
         ylab = "Average Width [m]")
    for(i in 1:length(TW)){
      lines(times / 365, TW[[i]])
      fit <- stats::lm(log(TW[[i]][2:length(times)]) ~ log(times[2:length(times)] / 365))
      lines(times[2:length(times)] / 365, exp(fit$fitted.values), lty = 2)
    }
  }else{
    return(TW)
  }
}

#' Transforms model outputs into data by each cross section
#'
#' @param data A matrix of model output data
#'
#' @return A list of data by cross section
data_by_XS <- function(data){
  n_XS <- sum(data[,1] == min(data[,1]))
  times <- unique(data[,1])

  output <- lapply(1:n_XS, function(x, data, n_XS){
    data[seq(x, nrow(data) - n_XS + x, n_XS), ]
  }, data, n_XS)

  return(output)
}


#' Plots specific stream power
#'
#' @param path Path to folder with model outputs
#' @param type `type = 1` plots stream power over time for each reach
#'   separately; `type = 2` plots stream power longitudinally by reach
#'
#' @importFrom dplyr %>%
#' @importFrom utils read.table
#' @importFrom graphics par plot grconvertY grconvertX lines points mtext
#'
#' @export
#'
plot_omega <- function(path = "", type = 1){

  omega <- read.table(paste0(path, "/Output stream power.txt"), header = FALSE) %>%
    dplyr::filter_(~ V1 > 0)

  times <- unique(omega$V1)

  XS <- 1:sum(omega$V1 == min(omega$V1))

  y_max <- max(omega[,2])
  y_min <- min(omega[,2][which(omega[,2] > 0)])

  colors <- cRamp_legend(length(times), "viridis")

  if (type == 1){
    par(mfrow = c(1,1), mar = c(4.5,4,2.5,0.5), oma = rep(0, 4))
    plot(NA, xlim = c(1, length(XS)), ylim = c(y_min, y_max),
         ylab = expression("Omega [W/"*m^2~"]"), xlab = "Reach", las = 1, log = "y", mgp = c(3, 1, 0))

    yvals <- grconvertY(c(1.0, 1.05), from = "npc", to = "user")
    xvals <- grconvertX(c(0.4, 0.6), from = "ndc", to = "user")
    color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c(expression(t[0]), "","",
                                                                    "",expression(t[f])),
                 align = "lt", gradient = "x",
                 rect.col = cRamp_legend(5, "viridis"), xpd = NA)

    for (i in 1:length(times)){
      subset <- omega[omega$V1 == times[i], 2]
      lines(subset, col = colors[i])
      points(subset, pch = 21, bg = colors[i])
    }
  }else if (type == 2){
    omega_XS <- data_by_XS(omega)
    par(mfrow = c(length(XS), 1), mar = c(0.5, 2, 0.5, 0.5), oma = c(3, 2.5, 3, 0))
    for (i in 1:length(XS)){
      n_reach <- sum(colMeans(omega_XS[[i]]) > 0) - 1
      colors <- cRamp_legend(n_reach, "viridis")
      ylim <- range(omega_XS[[i]][, 2:(n_reach + 1)])
      plot(NA, ylab = "", xlab = "", las = 1, ylim = ylim, xlim = range(times))
      for (j in 1:n_reach){
        lines(times, omega_XS[[i]][,j + 1], col = colors[j])
      }
    }
    mtext("Days", side = 1, line = 2, outer = TRUE)
    mtext(expression("Omega [W/"*m^2*"]"), side = 2, line = 1, outer = TRUE)

    yvals <- grconvertY(c(0.95, 0.97), from = "ndc", to = "user")
    xvals <- grconvertX(c(0.4, 0.6), from = "ndc", to = "user")
    color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                    "", "DS"),
                 align = "rb", gradient = "x",
                 rect.col = cRamp_legend(5, "viridis"), xpd = NA)
  }

}

#' Plots the channel network showing changes in bed elevation, width, and width-depth ratio over time
#'
#' @param path Path to folder with model outputs
#' @param XS Number of cross sections to label
#' @param pos Position of labels for cross sections (`right` or `left`)
#' @param years A numeric vectors of years of the simulation to plot
#' @param print Should the plot be printed to a file (defaults to `FALSE`)
#'
#' @importFrom utils read.table
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot rect points arrows text legend
#'
#' @export
#'
network_XS_plot <- function(path = "", XS = NULL,
                            pos = c("right", "right"),
                            years = c(1, 3, 5, 10, 20),
                            print = FALSE){

  bed_z <- read.table(paste0(path, "/Output z.txt"), header = FALSE, sep = "\t")
  link <- read.table(paste0(path, "/Input link.txt"), header = FALSE, sep = " ")
  #L <- read.table("Input length.txt", header = FALSE, sep = " ")

  #L <- as.array(L[,1])
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE, sep = "\t")

  #remove last 2 bed_z column
  bed_z <- bed_z[,1:(ncol(bed_z) - 1)]

  times <- unique(bed_z[,1])
  n_nodes <- ncol(link)

  angle <- pi/8

  n_xs <- apply(bed_z[1:n_nodes,2:ncol(bed_z)], 1, function(x){sum(x > 0)})
  x <- matrix(0, nrow = n_nodes + 1, ncol = max(n_xs))
  y <- x
  sgn <- rep(1, n_nodes)
  branched <- rep(1, n_nodes)

  for (i in seq(n_nodes, 1, -1)){
    if (i < n_nodes){
      index <- which(i == link, arr.ind = TRUE)[2]
      ang <- angle / branched[index]
      for (j in n_xs[i]:1){
        if (j == n_xs[i]){
          x[i, j] <- x[index, 1] + sgn[index] * sin(ang) * dx[i, j + 1]
          y[i, j] <- y[index, 1] + cos(ang) * dx[i, j + 1]
        }else{
          x[i, j] <- x[i, j + 1] + sgn[index] * sin(ang) * dx[i, j + 1]
          y[i, j] <- y[i, j + 1] + cos(ang) * dx[i, j + 1]
        }
      }

      sgn[index] <- -1 * sgn[index]
      branched[i] <- branched[index] + 1
    }else{
      index <- n_nodes + 1
      for (j in n_xs[i]:1){
        if (j == n_xs[i]){
          x[i, j] <- x[i + 1, 1]
          y[i, j] <- y[i + 1, 1] + dx[i, j + 1]
        }else{
          x[i, j] <- x[i, j + 1]
          y[i, j] <- y[i, j + 1] + dx[i, j + 1]
        }
      }
    }
  }

  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)

  #Get total bed elevation change
  dz_tot <- bed_z[bed_z[,1] == max(times),2:ncol(bed_z)] -
    bed_z[bed_z[,1] == 0,2:ncol(bed_z)]
  max_dz <- max(abs(dz_tot))

  legend_cols <- cRamp_legend(7, "RdYlBu")

  find_nearest <- function(find, data, shorten = FALSE){
    index <- which(abs(data - find) == min(abs(data - find)))
    if (shorten){
      index <- index[1]
    }
  }
  n_plot <- 5
  rows <- sapply(years, find_nearest, times / 365, TRUE)
  times_plot <- times[rows]
  line_cols <- cRamp_legend(n_plot, "viridis")

  if (print){
    png(paste0(path, "/Network XS Plot.png"), type = "cairo", units = "in", width = 6.5,
        height = 4.5, res = 500)
  }
  #par(mar = c(1, 1, 3, 1))
  par(mar = c(0.5, 0.5, 3, 0.5), mfrow = c(3, 5))
  for (k in 1:length(times_plot)){
    #calculated bed elevation changes
    dz <- bed_z[bed_z[,1] == times_plot[k],2:ncol(bed_z)] -
      bed_z[bed_z[,1] == 0,2:ncol(bed_z)]

    #Convert dz to vector to get colors - add maximum dz to get right scale
    dz_vector <- c(as.vector(t(dz)), max_dz, -max_dz)
    if (max_dz == 0){
      colors <- rep(legend_cols[4], length(dz_vector))
    }else{
      colors <- cRamp(dz_vector, "RdYlBu")
    }

    #Convert colors back to matrix - removing max dz value
    colors <- matrix(colors[1:(length(colors) - 2)], nrow = n_nodes, byrow = TRUE)

    plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

    for (i in seq(n_nodes, 1, -1)){
      if (i < n_nodes){
        index <- which(i == link, arr.ind = TRUE)[2]
      }
      else{
        index <- n_nodes + 1
      }
      for (j in n_xs[i]:1){
        # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
        #        lwd = 1)
        points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 1.2)
      }
      #Add arrow if XS specified
      if (k == 1){
        if (i %in% XS){
          index <- which(XS == i)
          x_diff1 <- (xmax - xmin) * 0.1
          x_diff2 <- (xmax - xmin) * 0.05
          x0 <- x[i, 1]
          y0 <- y[i, 1]
          if (pos[index] == "right"){
            x0 <- x0 + x_diff2
            x1 <- x0 + x_diff1
            txt_pos <- 4
          }else if(pos[index] == "left"){
            x0 <- x0 - x_diff2
            x1 <- x0 - x_diff1
            txt_pos <- 2
          }
          arrows(x0 = x0, x1 = x1, y0 = y0, y1 = y0, code = 1,
                 length = 0.05, col = "white")
          text(x = x1, y = y0, labels = LETTERS[index], cex = 0.7,
               col = "white", pos = txt_pos)
        }
      }
    }
    #add_label(-0.05, 0.07, paste0("(", letters[k], ")"), col = "white")
    if (k == 3){
      legend("top", legend = round(seq(-max_dz, max_dz, length.out = 7), 1),
             fill = legend_cols, bty = "n", title = "Elevation Change [m]",
             horiz = TRUE, inset = c(0, -0.35), xpd = NA)
    }
    legend("bottomright", legend = paste(round(times_plot[k] / 365, 1), "Years"),
           pch = NA, bty = "n", text.col = "white", cex = 1)
  }


  XS_output <- read.table(paste0(path, "/Output XS geometry all.txt"), header = FALSE,
                          sep= " ")

  #Get average channel width
  n_XS <- sum(XS_output[,1] == 0)
  times <- unique(XS_output[,1])

  XS <- lapply(1:n_XS, function(x, data, n_XS){
    data[seq(x, nrow(data) - n_XS + x, n_XS), ]
  }, XS_output, n_XS)

  #average of top and bottom width
  width <- lapply(XS, function(data){
    apply(data, 1, function(x){
      mean(c(x[9] - x[4], x[7] - x[6]))
    })
  })

  width <- do.call("rbind", width)

  dwidth <- apply(width, 2, function(x, init){
    (x - init) / init
  }, init = width[,1])
  dwidth[dwidth > 1] <- 1

  #Width_depth ratio
  w_h <- width_depth(path = path, plot = FALSE)

  dwh <- apply(w_h, 2, function(x, init){
    (x - init) / init
  }, init = w_h[,1])
  dwh[dwh > 1] <- 1

  plot_network <- function(data, times_plot, xmin, xmax, ymax, link,
                           x, y, type, times){
    for (k in 1:length(times_plot)){
      #Convert dz to vector to get colors - add maximum dz to get right scale
      max_dwh <- max(abs(data[,which(times %in% times_plot)]))
      dwh_vector <- c(data[,which(times == times_plot[k])], max_dwh, -max_dwh)
      if (max_dwh == 0){
        colors <- rep(legend_cols[4], length(dwh_vector))
      }else{
        colors <- cRamp(dwh_vector, "RdYlBu")
      }
      colors <- colors[1:(length(colors) - 2)]

      plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

      count <- length(colors)
      for (i in seq(n_nodes, 1, -1)){
        if (i < n_nodes){
          index <- which(i == link, arr.ind = TRUE)[2]
        }
        else{
          index <- n_nodes + 1
        }
        for (j in n_xs[i]:1){
          # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
          #        lwd = 1)
          points(x[i,j], y[i,j], pch = 16, col = colors[count], cex = 1.2)
          count <- count - 1
        }
        #Add numbers
        #text(xy[i, 1], xy[i, 2], i, pos = 4)
      }

      if (type == "width"){
        #add_label(-0.05, 0.07, paste0("(", letters[k + length(times_plot)], ")"),
        #          col = "white")
      }else if (type == "width_depth"){
        #add_label(-0.05, 0.07, paste0("(", letters[k + length(times_plot) * 2],
        #                              ")"), col = "white")
      }
      if (k == 3){
        if (type == "width"){
          title = "Percent Width Change"
          vals <- round(seq(-max_dwh * 100, max_dwh * 100, length.out = 7))
        }else if (type == "width_depth"){
          title = "Percent Width-Depth Ratio Change"
          vals <- round(seq(-max_dwh * 100, max_dwh * 100, length.out = 7))
        }
        legend("top", legend = vals,
               fill = legend_cols, bty = "n", title = title,
               horiz = TRUE, inset = c(0, -0.35), xpd = NA)
      }
      legend("bottomright", legend = paste(round(times_plot[k] / 365, 1), "Years"),
             pch = NA, bty = "n", text.col = "white", cex = 1)
    }
  }

  plot_network(data = dwidth, times_plot = times_plot, xmin = xmin, xmax = xmax,
               ymax = ymax, link = link, x = x, y = y, type = "width", times = times)
  plot_network(data = dwh, times_plot = times_plot, xmin = xmin, xmax = xmax,
               ymax = ymax, link = link, x = x, y = y, type = "width_depth",
               times = times)

  if (print){
    dev.off()
  }



}

#' Calculates channel width-depth ratio
#'
#' @param path Path to folder with model outputs
#' @param plot Should the data be plotted (defaults to `FALSE`)
#'
#' @return A dataframe with width-depth ratio by reach
#'
#' @importFrom utils read.table
#' @importFrom dplyr %>%
#' @importFrom graphics par lines mtext grconvertX grconvertY
width_depth <- function(path = "", plot = FALSE){

  XS_output <- read.table(paste0(path, "/Output XS geometry all.txt"), header = FALSE,
                          sep= " ")
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE) %>%
    dplyr::filter_(~ V1 == 0) %>%
    dplyr::select(-dplyr::one_of("V1"))
  n_XS_reach <- apply(dx, 1, function(x){sum(x>0)})

  #Get average channel width
  n_XS <- sum(XS_output[,1] == 0)
  times <- unique(XS_output[,1])

  XS <- lapply(1:n_XS, function(x, data, n_XS){
    data[seq(x, nrow(data) - n_XS + x, n_XS), ]
  }, XS_output, n_XS)

  #average of top and bottom width
  width <- lapply(XS, function(data){
    apply(data, 1, function(x){
      mean(c(x[9] - x[4], x[7] - x[6]))
    })
  })

  depth <- lapply(XS, function(data){
    apply(data, 1, function(x){
      mean(c(x[14] - x[16], x[19] - x[17]))
    })
  })

  w_h <- purrr::map2(width, depth, function(x, y){
    x/y
  })

  w_h <- do.call("rbind", w_h)
  rownames(w_h) <- unlist(purrr::map2(1:length(n_XS_reach), n_XS_reach, rep))

  if (plot){
    par(mfrow = c(2,2), mar = c(2, 2, 1, 0.5), oma = c(3, 3, 3, 0))
    xlim <- range(times)
    start_row <- 1
    for (i in 1:length(n_XS_reach)){
      reach_data <- w_h[start_row:(start_row + n_XS_reach[i] - 1), ]
      ylim <- range(reach_data)
      plot(NA, xlim = xlim, ylim = ylim, las = 1, ylab = "n", xlab = "n",
           main = paste("Reach", i))
      colors <- cRamp_legend(n_XS_reach[i], "viridis")

      for (j in 1:nrow(reach_data)){
        lines(times, reach_data[j, ], col = colors[j])
      }
      start_row <- start_row + n_XS_reach[i]

      if (i %in% seq(from = 1, length.out = 10, by = 4)){
        mtext("Days", side = 1, line = 2, outer = TRUE)
        mtext("width/Depth", side = 2, line = 1.5, outer=  TRUE)

        xvals <- grconvertX(x = c(0.4, 0.6), from = "ndc", to = "user")
        yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
        color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                        "","DS"),
                     align = "rb", gradient = "x",
                     rect.col = cRamp_legend(5, "viridis"), xpd = NA)
      }
    }
  }else{
    return(w_h)
  }
}

#' Plots sediment inflow and outflow over time
#'
#' @param path Path to folder with model outputs
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines legend
#'
#' @export
#'
sed_lines <- function(path = ""){
  sed <- read.table(paste0(path, "/Output sediment vol.txt"), header = TRUE)

  sed_in <- diff(sed$Bed_vol_in)
  sed_out <- diff(sed$Bed_vol_out)

  ymax <- max(c(sed_in, sed_out))

  colors <- plot_colors()
  par(mfrow = c(1,1), mar = c(4, 4, 0.5, 0.5), oma = rep(0, 4))
  plot(sed_in, type = "l", ylim = c(-ymax, ymax), col = colors[1], lwd = 2, las = 1,
       ylab = expression("Sediment Load ["*m^3*"/day]"), xlab = "Day", mgp = c(2.8, 0.8, 0))
  lines(-sed_out, col = colors[2], lwd = 2)
  legend("topright", legend = c("Sed in", "Sed out"), lty = 1, lwd = 2, col = colors,
         bty = "n")
}

#' Plots network showing changes in channel width, with uncertainty
#'
#' @param print Should the plot be printed to a file (defaults to `FALSE`)
#' @param n_MC Number of Monte Carlo simulations
#' @param path Path to folder with model outputs
#' @param MC_path Path to "MC Outputs" folder (only if different than `path`)
#' @param custom_sgn Specifies the direction each reach should be plotted (`-1`
#'   is left, `1` is right)
#' @param prob Numeric vector of percentiles of Monte Carlo results to plot in
#'   addition to the median (defaults to 0.05 and 0.95)
#' @param use_files Logical. Should results files that have been loaded be used
#'   (defaults to `TRUE`)
#'
#' @importFrom utils read.table
#' @importFrom stats quantile median
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot rect points legend grconvertX grconvertY text mtext
#'
#' @export
#'
width_MC_plot <- function(print = FALSE, n_MC, path = "",
                          MC_path = NULL, custom_sgn = NULL, prob = c(0.05, 0.95),
                          use_files = TRUE){

  link <- read.table(paste0(path, "/Input link.txt"), header = FALSE, sep = " ")
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE, sep = "\t")

  if (is.null(MC_path)){MC_path <- path}

  #Get MC results
  #check if it already exists
  if (exists("MC_width") & use_files){
    print("Using data already read from file.")
  }else{
    print("Reading data from file...")
    MC_width <- list()
    for (i in 0:(n_MC - 1)){
      MC_width[[i + 1]] <- read.table(paste0(MC_path, "/MC Outputs/Output width", i, ".txt"), header = FALSE,
                                      sep = "\t")
      #drop last column
      MC_width[[i + 1]] <- MC_width[[i + 1]][,1:(ncol(MC_width[[i + 1]]) - 1)]
    }
  }

  times <- unique(MC_width[[1]][,1])
  n_nodes <- ncol(link)

  dwidth <- list()
  for (i in 1:n_MC){
    dwidth[[i]] <- MC_width[[i]][MC_width[[i]][,1] == max(times),] -
      MC_width[[i]][MC_width[[i]][,1] == 0,]
  }

  n_xs <- apply(MC_width[[1]][1:n_nodes,2:ncol(MC_width[[1]])], 1, function(x){sum(x > 0)})

  coords <- make_network(n_nodes, n_xs, link, dx, custom_sgn)
  x <- coords$x
  y <- coords$y

  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)

  #Bed elevation change - 5th, 50th, and 95th percentiles
  dwidth_05 <- matrix(0,nrow = n_nodes, ncol = max(n_xs))
  dwidth_med <- dwidth_05
  dwidth_95 <- dwidth_05
  for (j in 1:n_nodes){
    for (k in 1:n_xs[j]){
      subset <- sapply(dwidth, '[', j, 1 + k)
      dwidth_05[j, k] <- quantile(subset, prob[1])
      dwidth_med[j, k] <- median(subset)
      dwidth_95[j, k] <- quantile(subset, prob[2])
    }
  }

  #Get max bed elevation change
  max_dwidth <- max(abs(c(dwidth_05, dwidth_med, dwidth_95)))
  legend_cols <- rev(cRamp_legend(7, "RdBu"))

  #Convert dz to vector to get colors - tack on max dz
  vector_05 <- -c(as.vector(t(dwidth_05)), max_dwidth, -max_dwidth)
  colors_05 <- cRamp(vector_05, "RdBu")[1:(length(vector_05) - 2)]
  vector_med <- -c(as.vector(t(dwidth_med)), max_dwidth, -max_dwidth)
  colors_med <- cRamp(vector_med, "RdBu")[1:(length(vector_med) - 2)]
  vector_95 <- -c(as.vector(t(dwidth_95)), max_dwidth, -max_dwidth)
  colors_95 <- cRamp(vector_95, "RdBu")[1:(length(vector_95) - 2)]

  #Convert colors back to matrix
  colors_05 <- matrix(colors_05, nrow = n_nodes, byrow = TRUE)
  colors_med <- matrix(colors_med, nrow = n_nodes, byrow = TRUE)
  colors_95 <- matrix(colors_95, nrow = n_nodes, byrow = TRUE)

  if (print){
    png(paste0(path, "/Width MC Network.png"), type = "cairo", units = "in",
        height = 4.5, width = 10, res = 800)
  }

  label <- c(paste0(sprintf("%.1f", prob[1]*100), "%"), "Median",
             paste0(sprintf("%.1f", prob[2]*100), "%"))

  par(mar = c(1, 1, 5, 1), mfrow = c(1,3))
  for (k in 1:3){
    if (k == 1){
      colors <- colors_05
    }else if (k == 2){
      colors <- colors_med
    }else{
      colors <- colors_95
    }

    plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

    for (i in seq(n_nodes, 1, -1)){
      if (i < n_nodes){
        index <- which(i == link, arr.ind = TRUE)[2]
      }
      else{
        index <- n_nodes + 1
      }
      for (j in n_xs[i]:1){
        # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
        #        lwd = 1)
        points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 2)
      }
      #Add numbers
      #text(xy[i, 1], xy[i, 2], i, pos = 4)
    }

    # legend("topleft", legend = round(seq(-max_dwidth, max_dwidth, length.out = 7), 2),
    #        fill = legend_cols, bty = "n", title = "Width Change [m]",
    #        cex = 1.3, text.col = "white")
    legend("bottomright", legend = label[k], pch = NA, bty = "n", cex = 2,
           text.col = "white")
    # add_label(0, 0.05, paste0("(", letters[k], ")"), col = "white",
    #           cex = 2)
    if (k == 2){
      xval <- grconvertX(seq(-0.1, 1.1, length.out = 7), from = "nfc")
      y_point <- rep(grconvertY(0.935, from = "nfc"), length(xval))
      y_lab <- rep(grconvertY(0.92, from = "nfc"), length(xval))

      points(xval, y_point, pch = 21, bg = legend_cols, cex = 2.5, xpd = NA)
      text(xval, y_lab, sprintf("%.1f", round(seq(-max_dwidth, max_dwidth, length.out = 7), 1)),
           pos = 1, xpd = NA, cex = 1.7)
      mtext("Width Change [m]", line = 3.5, side = 3, cex = 1.2)
    }
  }

  if (print){
    dev.off()
  }
}

#' Plots network showing changes in channel bed elevation, with uncertainty
#'
#' @param print Should the plot be printed to a file (defaults to `FALSE`)
#' @param n_MC Number of Monte Carlo simulations
#' @param path Path to folder with model outputs
#' @param MC_path Path to "MC Outputs" folder (only if different than `path`)
#' @param custom_sgn Specifies the direction each reach should be plotted (`-1`
#'   is left, `1` is right)
#' @param prob Numeric vector of percentiles of Monte Carlo results to plot in
#'   addition to the median (defaults to 0.05 and 0.95)
#' @param use_files Logical. Should results files that have been loaded be used
#'   (defaults to `TRUE`)
#'
#' @importFrom utils read.table
#' @importFrom stats quantile median
#' @importFrom graphics par plot rect points legend grconvertX grconvertY text mtext
#' @importFrom grDevices png dev.off
#'
#' @export
#'
dz_MC_plot <- function(print = FALSE, n_MC, path = "",
                       MC_path = NULL, custom_sgn = NULL, prob = c(0.05, 0.95),
                       use_files = TRUE){

  link <- read.table(paste0(path, "/Input link.txt"), header = FALSE, sep = " ")
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE, sep = "\t")

  if (is.null(MC_path)){MC_path <- path}

  #Get MC results
  #check if it already exists
  if (exists("MC_bed_z") & use_files){
    print("Using data already read from file.")
  }else{
    print("Reading data from file...")
    MC_bed_z <- list()
    for (i in 0:(n_MC - 1)){
      MC_bed_z[[i + 1]] <- read.table(paste0(MC_path, "/MC Outputs/Output z", i, ".txt"), header = FALSE,
                                      sep = "\t")
      #drop last column
      MC_bed_z[[i + 1]] <- MC_bed_z[[i + 1]][,1:(ncol(MC_bed_z[[i + 1]]) - 1)]
    }
  }

  #Assign MC_bed_z to environment
  assign("MC_bed_z", MC_bed_z, .GlobalEnv)

  times <- unique(MC_bed_z[[1]][,1])
  n_nodes <- ncol(link)

  dz <- list()
  for (i in 1:n_MC){
    dz[[i]] <- MC_bed_z[[i]][MC_bed_z[[i]][,1] == max(times),] -
      MC_bed_z[[i]][MC_bed_z[[i]][,1] == 0,]
  }

  n_xs <- apply(MC_bed_z[[1]][1:n_nodes,2:ncol(MC_bed_z[[1]])], 1, function(x){sum(x > 0)})

  coords <- make_network(n_nodes, n_xs, link, dx, custom_sgn)
  x <- coords$x
  y <- coords$y

  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)

  #Bed elevation change - 5th, 50th, and 95th percentiles
  dz_05 <- matrix(0,nrow = n_nodes, ncol = max(n_xs))
  dz_med <- dz_05
  dz_95 <- dz_05
  for (j in 1:n_nodes){
    for (k in 1:n_xs[j]){
      subset <- sapply(dz, '[', j, 1 + k)
      dz_05[j, k] <- quantile(subset, prob[1])
      dz_med[j, k] <- median(subset)
      dz_95[j, k] <- quantile(subset, prob[2])
    }
  }

  #Get max bed elevation change
  max_dz <- max(abs(c(dz_05, dz_med, dz_95)))
  legend_cols <- cRamp_legend(7, "RdBu")

  #Convert dz to vector to get colors - tack on max dz
  vector_05 <- c(as.vector(t(dz_05)), max_dz, -max_dz)
  colors_05 <- cRamp(vector_05, "RdBu")[1:(length(vector_05) - 2)]
  vector_med <- c(as.vector(t(dz_med)), max_dz, -max_dz)
  colors_med <- cRamp(vector_med, "RdBu")[1:(length(vector_med) - 2)]
  vector_95 <- c(as.vector(t(dz_95)), max_dz, -max_dz)
  colors_95 <- cRamp(vector_95, "RdBu")[1:(length(vector_95) - 2)]

  #Convert colors back to matrix
  colors_05 <- matrix(colors_05, nrow = n_nodes, byrow = TRUE)
  colors_med <- matrix(colors_med, nrow = n_nodes, byrow = TRUE)
  colors_95 <- matrix(colors_95, nrow = n_nodes, byrow = TRUE)

  if (print){
    png(paste0(path, "/dz MC Network.png"), type = "cairo", units = "in",
        height = 4.5, width = 10, res = 800)
  }

  label <- c(paste0(sprintf("%.1f", prob[1]*100), "%"), "Median",
             paste0(sprintf("%.1f", prob[2]*100), "%"))

  par(mar = c(1, 1, 5, 1), mfrow = c(1,3))
  for (k in 1:3){
    if (k == 1){
      colors <- colors_05
    }else if (k == 2){
      colors <- colors_med
    }else{
      colors <- colors_95
    }

    plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

    for (i in seq(n_nodes, 1, -1)){
      if (i < n_nodes){
        index <- which(i == link, arr.ind = TRUE)[2]
      }
      else{
        index <- n_nodes + 1
      }
      for (j in n_xs[i]:1){
        # lines(c(x[i, j], x[index, 1]), c(y[i, j], y[index, 1]), col = "black",
        #        lwd = 1)
        points(x[i,j], y[i,j], pch = 16, col = colors[i,j], cex = 2)
      }
      #Add numbers
      #text(xy[i, 1], xy[i, 2], i, pos = 4)
    }

    # legend("topleft", legend = round(seq(-max_dz, max_dz, length.out = 7), 2),
    #        fill = legend_cols, bty = "n", title = "Elevation Change [m]",
    #        cex = 1.6, text.col = "white")
    legend("bottomright", legend = label[k], pch = NA, bty = "n", cex = 2,
           text.col = "white")
    # add_label(0, 0.05, paste0("(", letters[k], ")"), col = "white",
    #           cex = 2)
    if (k == 2){
      xval <- grconvertX(seq(-0.1, 1.1, length.out = 7), from = "nfc")
      y_point <- rep(grconvertY(0.935, from = "nfc"), length(xval))
      y_lab <- rep(grconvertY(0.92, from = "nfc"), length(xval))

      points(xval, y_point, pch = 21, bg = legend_cols, cex = 2.5, xpd = NA)
      text(xval, y_lab, sprintf("%.1f", round(seq(-max_dz, max_dz, length.out = 7), 1)),
           pos = 1, xpd = NA, cex = 1.7)
      mtext("Elevation Change [m]", line = 3.5, side = 3, cex = 1.2)
    }
  }

  if (print){
    dev.off()
  }
}

#' Plots mass sediment loading rates by reach, with uncertainty
#'
#'
#' @param path Path to folder with model outputs
#' @param custom_sgn Specifies the direction each reach should be plotted (`-1`
#'   is left, `1` is right)
#' @param MC_path Path to "MC Outputs" folder (only if different than `path`)
#' @param n_MC Number of Monte Carlo simulations
#' @param units Character specifying units to be used in plot ("kg", "ton", or
#'   "1000 ton")
#' @param prob Numeric vector of percentiles of Monte Carlo results to plot in
#'   addition to the median (defaults to 0.05 and 0.95)
#' @param print Should the plot be printed to a file (defaults to `FALSE`)
#' @param type Whether sediment loading (`type = "sed"`, default) or pollutant
#'   loading (`type =  "p"`) should be plotted
#' @param use_files Logical. Should results files that have been loaded be used
#'   (defaults to `TRUE`)
#'
#' @importFrom utils read.table
#' @importFrom stats quantile median
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot rect lines legend grconvertX grconvertY points
#'   text mtext
#'
#' @export
#'
reach_loads <- function(path = "", custom_sgn = NULL,
                        MC_path = NULL, n_MC = 0, units = "ton", prob = c(0.05, 0.95),
                        print = FALSE, type = "sed", use_files = TRUE){

  if (is.null(MC_path)){MC_path <- path}

  bank_loads <- read.table(paste0(path, "/Output bank loading.txt"), header = TRUE)
  if (type == "sed"){
    sed_loads <- bank_loads[,which(substr(colnames(bank_loads), 1, 3) == "Sed")]
  }else{
    sed_loads <- bank_loads[,which(substr(colnames(bank_loads), 1, 1) == "P")]
  }

  link <- read.table(paste0(path, "/Input link.txt"), header = FALSE, sep = " ")
  dx <- read.table(paste0(path, "/Output dx.txt"), header = FALSE)
  dx_init <- dx[which(dx$V1 == 0),]
  lengths <- apply(dx_init, 1, sum) / 1000 #Reach length, km

  n_nodes <- ncol(bank_loads) / 2
  n_xs <- apply(dx[1:n_nodes,2:ncol(dx)], 1, function(x){sum(x > 0)})

  coords <- make_network(n_nodes, n_xs, link, dx, custom_sgn)
  x <- coords$x
  y <- coords$y
  xmax <- max(x, na.rm = TRUE)
  xmin <- min(x, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)

  scale <- dplyr::case_when(units == "kg" ~ 1,
                     units == "ton" ~ 1000,
                     units == "1000 ton" ~ 1e6)

  #Get MC results
  if (n_MC > 0){
    #check if it already exists
    if (exists("MC_loading") & use_files){
      print("Using data already read from file.")
    }else{
      print("Reading data from file...")
      MC_loading <- list()
      for (i in 0:(n_MC - 1)){
        MC_loading[[i + 1]] <- read.table(paste0(MC_path, "/MC Outputs/Output bank loading", i, ".txt"), header = TRUE)
        #Keep sediment or P
        if (type == "sed"){
          MC_loading[[i + 1]] <- MC_loading[[i + 1]][,1:n_nodes]
        }else{
          MC_loading[[i + 1]] <- MC_loading[[i + 1]][,(n_nodes + 1):ncol(MC_loading[[i + 1]])]
        }
      }
    }

    #Assign MC_bed_z to environment
    assign("MC_loading", MC_loading, .GlobalEnv)

    loads_by_reach <- lapply(MC_loading, function(x, lengths, scale){
      apply(x, 2, sum) / scale / lengths / (nrow(x) / 365) #ton (or kg) / km / yr
    }, lengths, scale)
    loads_by_reach <- do.call("rbind", loads_by_reach)
    loads_05 <- apply(loads_by_reach, 2, quantile, prob[1])
    loads_med <- apply(loads_by_reach, 2, median)
    loads_95 <- apply(loads_by_reach, 2, quantile, prob[2])
    max_load <- max(c(loads_05, loads_med, loads_95))

    colors_05 <- cRamp(c(loads_05, max_load, 0), "Reds")
    colors_med <- cRamp(c(loads_med, max_load, 0), "Reds")
    colors_95 <- cRamp(c(loads_95, max_load, 0), "Reds")

    colors_legend <- cRamp_legend(7, "Reds")
    legend_lab <- dplyr::case_when(units == "kg" ~ "kg/km/yr",
                            units == "ton" ~ "ton/km/yr",
                            units == "1000 ton" ~ "1000 ton/km/yr")

    if (print){
      png(paste0(path, "/Loading MC Network.png"), type = "cairo", units = "in",
          height = 4.5, width = 10, res = 800)
    }

    label <- c(paste0(sprintf("%.1f", prob[1]*100), "%"), "Median",
               paste0(sprintf("%.1f", prob[2]*100), "%"))

    par(mar = c(1, 1, 5, 1), mfrow = c(1,3))
    for (k in 1:3){
      if (k == 1){
        colors <- colors_05
      }else if (k == 2){
        colors <- colors_med
      }else{
        colors <- colors_95
      }

      plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n",
           xaxt = "n", xlab = "", ylab = "")
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

      for (i in 1:n_nodes){
        lines(x[i,1:n_xs[i]], y[i,1:n_xs[i]], lwd = 4, col = colors[i])
      }

      legend("bottomright", legend = label[k], pch = NA, bty = "n", cex = 2,
             text.col = "white")
      # add_label(0, 0.05, paste0("(", letters[k], ")"), col = "white",
      #           cex = 2)
      if (k == 2){
        xval <- grconvertX(seq(-0.1, 1.1, length.out = 7), from = "nfc")
        y_point <- rep(grconvertY(0.935, from = "nfc"), length(xval))
        y_lab <- rep(grconvertY(0.92, from = "nfc"), length(xval))

        points(xval, y_point, pch = 22, bg = colors_legend, cex = 2.5, xpd = NA)
        text(xval, y_lab, sprintf("%.0f", round(seq(0, max_load, length.out = 7), 0)),
             pos = 1, xpd = NA, cex = 1.7)
        if (type == "sed"){
          mtext(paste0("Sediment Loading [", legend_lab, "]"), line = 3.5,
                side = 3, cex = 1.2)
        }else{
          mtext(paste0("Pollutant Loading [", legend_lab, "]"), line = 3.5,
                side = 3, cex = 1.2)
        }
      }
    }
    if (print){
      dev.off()
    }
  }else{

    loads_by_reach <- apply(sed_loads, 2, sum) / scale / lengths / (nrow(sed_loads) / 365) #ton (or kg) / km / yr
    max_load <- max(loads_by_reach)

    colors_legend <- cRamp_legend(7, "Reds")
    legend_lab <- dplyr::case_when(units == "kg" ~ "kg/km/yr",
                            units == "ton" ~ "ton/km/yr",
                            units == "1000 ton" ~ "1000 ton/km/yr")

    if (max_load == 0){
      colors <- rep(colors_legend[1], length(loads_by_reach))
    }else{
      colors <- cRamp(c(loads_by_reach, max_load, 0), "Reds")
    }

    if (print){
      png(paste0(path, "/Loading Network.png"), type = "cairo", units = "in",
          height = 4.5, width = 4.5, res = 800)
    }

    par(mar = c(1, 1, 5, 1), mfrow = c(1,1))

    plot(NA, xlim = c(xmin, xmax * 1.05), ylim = c(0, ymax), yaxt = "n",
         xaxt = "n", xlab = "", ylab = "")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray40")

    for (i in 1:n_nodes){
      lines(x[i,1:n_xs[i]], y[i,1:n_xs[i]], lwd = 4, col = colors[i])
    }

    xval <- grconvertX(seq(0.1, 0.9, length.out = 7), from = "nfc")
    y_point <- rep(grconvertY(0.89, from = "nfc"), length(xval))
    y_lab <- rep(grconvertY(0.87, from = "nfc"), length(xval))

    points(xval, y_point, pch = 22, bg = colors_legend, cex = 1.8, xpd = NA)
    text(xval, y_lab, sprintf("%.0f", round(seq(0, max_load, length.out = 7), 0)),
         pos = 1, xpd = NA, cex = 1.2)
    if (type == "sed"){
      mtext(paste0("Sediment Loading [", legend_lab, "]"), line = 3.5,
            side = 3, cex = 1.2)
    }else{
      mtext(paste0("Pollutant Loading [", legend_lab, "]"), line = 3.5,
            side = 3, cex = 1.2)
    }

    if (print){
      dev.off()
    }

  }

}

#' Plots cumulative or daily pollutant loading
#'
#' @param path Path to folder with model outputs
#' @param type `type = 1` plots cumulative loads, `type = 2` plots daily loads
#' @param returnvals Should the cumulative load values be returned (logical)
#'
#' @importFrom utils read.table
#' @importFrom graphics par plot lines grconvertX grconvertY
#'
#' @export
#'
pollutant_loading <- function(path = "", type = 1, returnvals = FALSE){

  loads <- read.table(file.path(path, "Output bank loading.txt"), header = TRUE)
  sed_loads <- loads[,which(substr(colnames(loads), 1, 3) == "Sed")]
  P_loads <- loads[,which(substr(colnames(loads), 1, 1) == "P")]

  total_sed <- apply(sed_loads, 1, sum)
  max_sed <- sum(total_sed)

  total_P <- apply(P_loads, 1, sum)
  max_P <- sum(total_P)

  get_scales <- function(max){
    scale <- dplyr::case_when(log10(max) < 3 ~ 1,
                       log10(max) < 6 ~ 1000,
                       log10(max) < 9 ~ 1e6,
                       log10(max) < 12 ~ 1e9)
    units <- dplyr::case_when(log10(max) < 3 ~ "kg",
                       log10(max) < 6 ~ "Mg",
                       log10(max) < 9 ~ "Gg",
                       log10(max) < 12 ~ "Tg")

    return(list(scale = scale, units = units))
  }

  sed_scale <- get_scales(max_sed)

  if (type == 1){
    par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1.5, 0.5), oma = c(0, 0, 1.5, 0),
        mgp = c(2, 0.8, 0))
    colors <- cRamp_legend(ncol(sed_loads), "viridis")
    plot(1:length(total_sed) / 365, cumsum(total_sed) / sed_scale$scale, type = "l",
         ylab = bquote("Sediment Load ["*.(a)*"]", list(a = sed_scale$units)), xlab = "Years",
         las = 1, main = "Sediment Loads", lwd = 2)
    for (i in 1:ncol(sed_loads)){
      lines(1:nrow(sed_loads) / 365, cumsum(sed_loads[,i]) / sed_scale$scale, col = colors[i], lwd = 2)
    }

    P_scale <- get_scales(max_P)
    plot(1:length(total_P) / 365, cumsum(total_P) / P_scale$scale, type = "l",
         ylab = bquote("Pollutant Load ["*.(a)*"]", list(a = P_scale$units)), xlab = "Years",
         las = 1, main = "Pollutant Loads", lwd = 2)
    for (i in 1:ncol(P_loads)){
      lines(1:nrow(P_loads) / 365, cumsum(P_loads[,i]) / P_scale$scale, col = colors[i], lwd = 2)
    }

    #Legend
    xvals <- grconvertX(x = c(0.4, 0.6), from = "ndc", to = "user")
    yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
    color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                    "","DS"),
                 align = "rb", gradient = "x",
                 rect.col = cRamp_legend(5, "viridis"), xpd = NA)
  }
  else if (type == 2){
    sed_loads[sed_loads == 0] <- NA
    P_loads[P_loads == 0] <- NA
    par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1.5, 0.5), oma = c(0, 0, 1.5, 0),
        mgp = c(2, 0.8, 0))
    colors <- cRamp_legend(ncol(sed_loads), "viridis")
    sed_scale <- get_scales(max(total_sed, na.rm = TRUE))
    plot(1:length(total_sed) / 365, total_sed / sed_scale$scale, type = "l",
         ylab = bquote("Sediment Load ["*.(a)*"]", list(a = sed_scale$units)), xlab = "Years",
         las = 1, main = "Sediment Loads", lwd = 2)
    for (i in 1:ncol(sed_loads)){
      lines(1:nrow(sed_loads) / 365, sed_loads[,i] / sed_scale$scale, col = colors[i], lwd = 2)
    }

    P_scale <- get_scales(max(P_loads, na.rm = TRUE))
    plot(1:length(total_P) / 365, total_P / P_scale$scale, type = "l",
         ylab = bquote("Pollutant Load ["*.(a)*"]", list(a = P_scale$units)), xlab = "Years",
         las = 1, main = "Pollutant Loads", lwd = 2)
    for (i in 1:ncol(P_loads)){
      lines(1:nrow(P_loads) / 365, P_loads[,i] / P_scale$scale, col = colors[i], lwd = 2)
    }

    #Legend
    xvals <- grconvertX(x = c(0.4, 0.6), from = "ndc", to = "user")
    yvals <- grconvertY(y = c(0.96, 0.98), from = "ndc", to = "user")
    color.legend(xvals[1], yvals[1], xvals[2], yvals[2], legend = c("US", "","",
                                                                    "","DS"),
                 align = "rb", gradient = "x",
                 rect.col = cRamp_legend(5, "viridis"), xpd = NA)
  }

  if (returnvals){
    return(data.frame(Sed = total_sed,
                      P = total_P))
  }
}
