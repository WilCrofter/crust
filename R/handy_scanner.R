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
  source("R/box_scanner.R")
}

# Create a phantom, two opposing probes, and a grid.
#' @param n the number of transducers per probe
#' @param z the z coordinate of the probe plane
#' @param npix the number of cells in a grid which covers the horizontal section
#' @return a list of the phantom, the probes (xmitr and rcvr,) an alignment utility (align,) npix, and z. 
handySetup <- function(n, z, npix=n){
  phantom <- newBoxPhantom(dimensions=c(30, 30, 30), 
                           ctr_alcohol=c(8, 15-1.5, 18), r_alcohol=6, 
                           ctr_water=c(24, 16+1.5, 12), r_water=5,
                           alignment_in_world = diag(1,4,4))
  xmitr <- newProbe(n=n, spacing = 30/n)
  rcvr <- newProbe(n=n, spacing = 30/n)
  align <- alignProbes(phantom, c(15,15,z), c(1,0,0), xmitr, rcvr)
  list(phantom=phantom, xmitr=xmitr, rcvr=rcvr, align=align, npix=npix, z=z)
}

# Create a figure showing the horizontal section, the transmitter and receiver
# positions, and the grid.
showSetup <- function(setup){
  plotSectionAndArrays(setup$phantom, setup$align$transmitters, setup$align$receivers,
                       setup$z, by=ceiling(setup$npix/8), legends=TRUE)
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
