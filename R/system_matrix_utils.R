# Trickery to handle sourcing from an Rmd
if(isTRUE(file.info("R")$isdir)){
  # working directory contains R
  source("R/utilities.R")
} else if(isTRUE(file.info("../R")$isdir)){
  # parent of working directory contains R
  source("../R/utilities.R")
}

# Specify a test image by number of rows, number of columns, spacing,
# and a default value for slowness.
testImage <- function(nrows, ncols, spacing, slowness=1/1500){
  img <- matrix(slowness, nrows, ncols)
  attr(img, "spacing") <- spacing
  attr(img, "slowness") <- slowness
  attr(img, "class") <- c("test image", attr(img, "class"))
  img
}

# Overlay grid lines on a plot of a test image
gridOverlay <- function(test_image, col="lightblue", lwd=2, ...){
  if(!is(test_image, "test image"))stop("Argument is not a test image.")
  spacing <- attr(test_image, "spacing")
  oset <- -.5*spacing
  xmax <- spacing*(nrow(test_image)-.5)
  xgrids <- seq(oset, xmax, by=spacing)
  ymax <- spacing*(ncol(test_image)-.5)
  ygrids <- seq(oset, ymax, by=spacing)
  segments(xgrids, rep(oset,nrow(test_image)), xgrids, rep(ymax, nrow(test_image)),
           col=col, lwd=lwd, ...)
  segments(rep(oset,ncol(test_image)), ygrids, rep(xmax, ncol(test_image)), ygrids,
           col=col, lwd=lwd, ...)
}

# Return the y coordinates of n transducers aligned at the left and right
# of a test image.
#
# The transducer array will span the image vertically from the middle
# of the bottom pixel to the middle of the top pixel. In the real world,
# of course, transducer spacing will be fixed and pixels chosen for
# analytic convenience. The convention used here is an artifice suited
# to abstract analysis of system matrices.
probeYs <- function(n, test_image){
  if(!is(test_image, "test image"))stop("Argument is not a test image.")
  if(!is.numeric(n) | length(n) != 1)stop("Argument, n, is not a numerical scalar.")
  spacing <- attr(test_image,"spacing") # pixel spacing
  height <- ncol(test_image)*spacing    # columns index the y dimension
  ys <- seq(spacing/2, height-spacing/2, length.out = n)
  attr(ys, "class") <- c("probe ys", attr(ys, "class"))
  ys
}

# Overlay probe positions on a plot of a test image
probeOverlay <- function(probe_ys, test_image, col = c("green3", "red"), pch=19, ...){
  if(!is(probe_ys, "probe ys"))stop("First argument is not an object of class 'probe ys'")
  if(!is(test_image, "test image"))stop("Second argument is not a test image.")
  n <- length(probe_ys)
  spacing <- attr(test_image, "spacing")
  width <- nrow(test_image)*attr(test_image, "spacing")
  points(rep(-spacing/2,n), probe_ys-spacing/2, col=head(col,1), pch=pch, ...)
  points(rep(width-spacing/2, n), probe_ys-spacing/2, col=tail(col,1), pch=pch, ...)
}

# Overlay a transmitter-to-receiver path on a plot of a test image
pathOverlay <- function(probe_ys, test_image, transmitter, receiver, col="black", lwd=2, ...){
  if(!is(probe_ys, "probe ys"))stop("First argument is not an object of class 'probe ys'")
  if(!is(test_image, "test image"))stop("Second argument is not a test image.")
  n <- length(probe_ys)
  if(!is.numeric(transmitter) | !(length(transmitter)==1) |
     !(transmitter > 0) | ! transmitter <= n)stop(
    paste("Argument 'transmitter' must be an integer between 1 and ", n))
  if(!is.numeric(receiver) | !(length(receiver)==1) |
     !(receiver > 0) | !(receiver <= n))stop(
    paste("Argument 'receiver' must be an integer between 1 and ", n))
  delta <- attr(test_image, "spacing")/2 # for plot
  y1 <- probe_ys[round(transmitter)]
  y2 <- probe_ys[round(receiver)]
  width <- nrow(test_image)*attr(test_image, "spacing")
  segments(-delta, y1+delta, width-delta, y2+delta, col=col, lwd=lwd, ...)
}

# Representative values for speeds of sound in tissues of interest
# Reference: 
# [Nebeker and Nelson](http://www.jultrasoundmed.org/content/31/9/1389.full) 
# "Mean values from published sound speed reports are as follows: fat, 1478 m/s; 
# glandular breast, 1510 m/s; benign breast tumors, 1513 m/s; 
# and malignant breast tumors, 1548 m/s"
speed <- c(fat=1478, gland=1510, benign=1513, malignant=1548, 
           silicone=1000, alcohol=1180, water=1480)

# Returns a system matrix, or Matrix, for the given probe positions and test image.
# If parameter 'check' is TRUE, entries will be checked to ensure they sum to
# appropriate path lengths, i.e., that no pixels were missed due to precision issues.
# If paramter 'MatrixPackage' is TRUE, the returned object will be a Matrix
# from the package of the same name, which is advantageous in the case of large,
# sparse matrices. Otherwise, an ordinary R matrix will be returned.
sysMat <- function(probe_ys, test_image, check=FALSE, MatrixPackage=FALSE){
  if(!is(probe_ys, "probe ys"))stop("First argument is not an object of class 'probe ys'")
  if(!is(test_image, "test image"))stop("Second argument is not a test image.")
  if(MatrixPackage & is(try(find.package("Matrix"), silent=TRUE),"try-error")){
    if(is(err, "try-error")){
      stop("Argument MatrixPackage is TRUE, but the Matrix package is not installed.")
    } else {
      library(Matrix)
    }
  }
  spacing <- attr(test_image, "spacing")
  nprobe <- length(probe_ys)
  nx <- nrow(test_image)
  ny <- ncol(test_image)
  npix <- nx*ny
  x0 <- 0
  x1 <- nx*spacing
  S <- if(MatrixPackage){
    Matrix(0, nrow=nprobe^2, ncol=npix)
  } else {
    matrix(0, nrow=nprobe^2, ncol=npix)
  }
  for (i in 1:nprobe) 
    for (j in 1:nprobe){
      u <- c(x0, probe_ys[i], 0)
      v <- c(x1, probe_ys[j], 0)
      segs <- segmentLengths(nx, ny, u, v, spacing, zero_origin=TRUE)
      if(check & 
         !isTRUE(all.equal(sum(segs[,"segment_length"]), 
                           sqrt(x1^2 + (v[2]-u[2])^2)))){
        stop(paste("Line segment(s) missing in path",i, j))
      }
      for (k in 1:length(segs$segment_length)) {
        idx <- (segs$y_index[k]-1) * nx + segs$x_index[k]
        S[(j-1)*nprobe + i, idx] <- segs$segment_length[k]
      }
    }
  S
}
