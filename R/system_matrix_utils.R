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

# Modifies createSij to work independently of handy scanner.
formSij <- function(setup){
  npix <- setup$npix
  numrows <- setup$n*(setup$n-1)
  Sij <- matrix(0,nrow=numrows,ncol=npix^2)
  for (i in 1:(setup$n-1)) 
    for (j in 1:setup$n){
      row <- pixLengths(i,j,setup)
      for (k in 1:length(row$segment_length)) {
        idx <- (row$y_index[k]-1) * npix + row$x_index[k]
        # Sij[(i-1)*setup$n + j,idx] <- row$segment_length[k]
        Sij[(j-1)*(setup$n-1) + i, idx] <- row$segment_length[k]
      }
    }
  #now do some checks on Sij
  rsums <- rowSums(Sij)
  
  for (i in 1:numrows) {
    nz <- sum(Sij[i,]>0)
    if (nz<npix)stop("too few nonzero elements in row ",i)
    cl <- floor((i-1)/(setup$n-1))+1
    rw <- (i-1)%%(setup$n-1) + 1
    dist <- sqrt((setup$align$transmitters[1,rw]-setup$align$receivers[1,cl])^2 +
                   (setup$align$transmitters[2,rw]-setup$align$receivers[2,cl])^2 +
                   (setup$align$transmitters[3,rw]-setup$align$receivers[3,cl])^2)
    if (!isTRUE(all.equal(dist,rsums[i])))stop("pathlengths not equal in row ",i," ",dist," ",rsums[i])
  }
  Sij
}
