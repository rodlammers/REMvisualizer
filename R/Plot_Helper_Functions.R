#' @importFrom utils read.table
#' @importFrom graphics par text
#'
add_label <- function(xfrac, yfrac, label, pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, xpd = NA, ...)
}

range02 <- function(x)(x + max(abs(x)))/(2 * max(abs(x)))

range01 <- function(x)(x-min(x))/diff(range(c(x, max(abs(x)))))

range03 <- function(x)(x-min(x))/diff(range(c(x, max(x))))

#' Creates a color ramp from a set of values
#'
#' @param x Series of values to create color ramp for
#' @param palette Palette of colors to use (either `viridis` or a palette from
#'   `RColorBrewer`)
#' @param alpha Transparency factor (defaults to `1`)
#'
#' @return Set of colors corresponding to each supplied value
#'
#' @importFrom grDevices colorRamp rgb adjustcolor
#'
cRamp <- function(x, palette, alpha = 1){
  #cols <- colorRamp(c("red", "gray40", "blue"))(range02(x))
  #if (sum(x < 0) == length(x)){
  range <- range03(x)
  # }else{
  #   range <- range01(x)
  # }
  if (palette == "viridis"){
    cols <- colorRamp(viridis::viridis_pal()(10))(range)
    cols <- apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
    cols <- adjustcolor(cols, alpha.f = alpha)
  }else if(palette == "custom"){
    cols <- colorRamp(colorspace::diverge_hcl(n = length(range), c = c(100, 0),
                                              l = c(50, 90)))(range)
    cols <- apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
    cols <- adjustcolor(cols, alpha.f = alpha)
  }else{
    cols <- suppressWarnings(colorRamp(RColorBrewer::brewer.pal(11, palette))(range))
    cols <- apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
    cols <- adjustcolor(cols, alpha.f = alpha)
  }
  return(cols)
}

#' Creates a color ramp of specified length
#'
#' @param x Length of color palette
#' @param palette Palette of colors to use (either `viridis` or a palette from
#'   `RColorBrewer`)
#' @param alpha Transparency factor (defaults to `1`)
#'
#' @return Set of colors of length `x`
#'
cRamp_legend <- function(x, palette, alpha = 1){
  range <- seq(0, 1, length.out = x)

  if (palette != "viridis"){
    cols <- suppressWarnings(colorRamp(RColorBrewer::brewer.pal(11, palette))(range))
    cols <- apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
    cols <- adjustcolor(cols, alpha.f = alpha)
  }else{
    cols <- colorRamp(viridis::viridis_pal()(10))(range)
    cols <- apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
    cols <- adjustcolor(cols, alpha.f = alpha)
  }

  return(cols)
}

#' A series of nice, color-blind friendly, colors for plotting
#'
#' @param alpha Transparency factor (defaults to `1`)
#' @param plot Should the colors be plotted (defaults to `FALSE`)
#'
#' @return A set of eight named colors
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics par lines points text
#'
plot_colors <- function(alpha = 1, plot = FALSE){
  colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
              "#D55E00", "#CC79A7", "#999999")

  colors <- adjustcolor(colors, alpha.f = alpha)

  names(colors) <- c("orange", "sky blue", "bluish green", "yellow", "blue",
                     "vermillion", "reddish purple", "gray")

  if (plot){
    par(mfrow = c(1,1), mar = c(0,0,0,0))
    plot(NA, xlim = c(0,1), ylim = c(0, length(colors) + 1), xaxt = "n",
         yaxt = "n", xlab = "", ylab = "", axes = FALSE)
    for (i in 1:length(colors)){
      lines(0:1, rep(i, 2), lwd = 3, col = colors[i])
      points(0.5, i, pch = 21, bg = colors[i], cex = 2)
      text(0.2, i, names(colors)[i], pos = 3)
    }
  }

  return(colors)
}

#' @importFrom graphics par rect
gradient.rect<-function(xleft,ybottom,xright,ytop,reds,greens,blues,
                        col=NULL,nslices=50,gradient="x",border=par("fg")) {

  if(is.null(col)) col <- plotrix::color.gradient(reds, greens, blues, nslices)
  else nslices<-length(col)
  nrect<-max(unlist(lapply(list(xleft,ybottom,xright,ytop),length)))
  if(nrect > 1) {
    if(length(xleft) < nrect) xleft<-rep(xleft,length.out=nrect)
    if(length(ybottom) < nrect) ybottom<-rep(ybottom,length.out=nrect)
    if(length(xright) < nrect) xright<-rep(xright,length.out=nrect)
    if(length(ytop) < nrect) ytop<-rep(ytop,length.out=nrect)
    for(i in 1:nrect)
      gradient.rect(xleft[i],ybottom[i],xright[i],ytop[i],
                    reds,greens,blues,col,nslices,gradient,border=border)
  }
  else {
    if (gradient == "x") {
      xinc <- (xright - xleft)/nslices
      xlefts <- seq(xleft, xright - xinc, length = nslices)
      xrights <- xlefts + xinc
      rect(xlefts,ybottom,xrights,ytop,col=col,lty=0, xpd = NA)
      rect(xlefts[1],ybottom,xrights[nslices],ytop,border=border, xpd = NA)
    }
    else {
      yinc <- (ytop - ybottom)/nslices
      ybottoms <- seq(ybottom, ytop - yinc, length = nslices)
      ytops <- ybottoms + yinc
      rect(xleft,ybottoms,xright,ytops,col=col,lty=0, xpd = NA)
      rect(xleft,ybottoms[1],xright,ytops[nslices],border=border, xpd = NA)
    }
  }
  invisible(col)
}

#' @importFrom graphics par strheight strwidth text
#'
color.legend<-function (xl,yb,xr,yt,legend,rect.col,cex=1,align="lt",
                        gradient="x",...) {

  oldcex<-par("cex")
  par(xpd=TRUE,cex=cex)
  gradient.rect(xl,yb,xr,yt,col=rect.col,nslices=length(rect.col),
                gradient=gradient)
  if(gradient == "x") {
    xsqueeze<-(xr-xl)/(2*length(rect.col))
    textx<-seq(xl+xsqueeze,xr-xsqueeze,length.out=length(legend))
    if(match(align,"rb",0)) {
      texty<-yb-0.2*strheight("O")
      textadj<-c(0.5,1)
    }
    else {
      # assume it's the default
      texty<-yt+0.2*strheight("O")
      textadj<-c(0.5,0)
    }
  }
  else {
    ysqueeze<-(yt-yb)/(2*length(rect.col))
    texty<-seq(yb+ysqueeze,yt-ysqueeze,length.out=length(legend))
    if(match(align,"rb",0)) {
      textx<-xr+0.2*strwidth("O")
      textadj<-c(0,0.5)
    }
    else {
      # assume it's the default
      textx<-xl-0.2*strwidth("O")
      textadj<-c(1,0.5)
    }
  }
  text(textx,texty,labels=legend,adj=textadj,...)
  par(xpd=FALSE,cex=oldcex)
}
