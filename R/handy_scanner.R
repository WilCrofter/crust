# Code to create a phantom and opposing probes, create associated scan data, and return elements of
# a system matrix. The phantom will be a cube, 30 mm on a side, with embedded water and alcohol spheres.
# The probes will be aligned, initially, in a horizontal plane at a given z coordinate. Their
# alignments may be changed if desired. The number of transducers can be specified. They will 
# be spaced to cover a 30 mm side.

# Trickery to handle sourcing from an Rmd
if(isTRUE(file.info("R")$isdir)){
  # working directory contains R
  source("R/box_scanner.R")
} else if(isTRUE(file.info("../R")$isdir)){
  # parent of working directory contains R
  source("../R/box_scanner.R")
}

# Create a phantom, two opposing probes, and a grid.
#' @param n the number of transducers per probe
#' @param z the z coordinate of the probe plane
#' @param npix the number of cells in a grid which covers the horizontal section
#' @return a list of the phantom, the probes (xmitr and rcvr,) an alignment utility (align,) npix, and z. 
handySetup <- function(n, z, npix=n-1){
  phantom <- newBoxPhantom(dimensions=c(30, 30, 30), 
                           ctr_alcohol=c(8, 15-1.5, 18), r_alcohol=6, 
                           ctr_water=c(24, 16+1.5, 12), r_water=5,
                           alignment_in_world = diag(1,4,4))
  xmitr <- newProbe(n=n, spacing = 30/n)
  rcvr <- newProbe(n=n, spacing = 30/n)
  setup <- list(phantom=phantom, xmitr=xmitr, rcvr=rcvr, align=NULL, n=n, npix=npix, z=z)
  setup <- alignX(setup)
  setup
}

alignX <- function(setup){
  pt <- c(setup$phantom$dimensions[1:2]/2,setup$z)
  setup$align <- alignProbes(setup$phantom, pt, c(1,0,0), setup$xmitr, setup$rcvr)
  setup
}

alignY <- function(setup){
  pt <- c(setup$phantom$dimensions[1:2]/2,setup$z)
  setup$align <-alignProbes(setup$phantom, pt, c(0,1,0), setup$xmitr, setup$rcvr)
  setup
}

# Create a figure showing the horizontal section, the transmitter and receiver
# positions, and the grid.
showSetup <- function(setup,legends=TRUE){
  plotSectionAndArrays(setup$phantom, setup$align$transmitters, setup$align$receivers,
                       setup$z, by=ceiling(setup$npix/8), legends=legends)
  plotGrid(setup$npix, setup$npix, spacing=30/setup$npix, add=TRUE, col="green")
}

# Add to an existing figure the path from transmitter i to receiver j.
addPath <- function(i, j, setup, col="magenta", lwd=2){
  u <- setup$align$transmitters[,i]
  v <- setup$align$receivers[,j]
  segments(u[1], u[2], v[1], v[2], col=col, lwd=lwd)
}

# Return a matrix of times of flight. Row indices correspond to
# transmitters, column indices to receivers.
# NOTE: for a color plot, use plotScan on the result.
scanSetup <- function(setup){
  doScan(setup$phantom, setup$align$transmitters, setup$align$receivers)
}

# Return the intersected pixels and lengths of the intersections for the
# path between transmitter i and receiver j.
pixLengths <- function(i, j, setup){
  segmentLengths(setup$npix, setup$npix, 
                 setup$align$transmitters[,i], 
                 setup$align$receivers[,j], 
                 spacing = 30/setup$npix)
}

# Create the Sij matrix by calling pixLengths in double loop
createSij <- function(setup){
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
