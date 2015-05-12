# Automates scan simulations as outlined in www/sim_scan_1.Rmd

# Trickery to handle sourcing from an Rmd
if(isTRUE(file.info("R")$isdir)){
  # working directory contains R
  source("R/utilities.R")
  source("R/box_phantom.R")
  source("R/probe.R")
} else if(isTRUE(file.info("../R")$isdir)){
  # parent of working directory contains R
  source("../R/utilities.R")
  source("../R/box_phantom.R")
  source("../R/probe.R")
}

#' Given a box phantom, an interior point, and a direction, align two probes in
#' contralateral horizontal positions on the phantom boundaries centered on the
#' line defined by the interior point and direction.
#' @param phantom a box phantom object
#' @param point a point in the phantom, specified in the phantom's internal reference frame.
#' @param direction a vector specified in the phantom's internal frame
#' @param xmitr the transmitting probe object
#' @param rcvr the receiving probe object
#' @return positions of transducer array in the phantom's internal reference frame and z for axial section.
alignProbes <- function(phantom, point, direction, xmitr, rcvr){
  if(!is(phantom, "box phantom"))stop("Argument phantom is not a box phantom object.")
  if(!is(xmitr, "probe"))stop("xmitr is not a probe object.")
  if(!is(rcvr, "probe"))stop("rcvr is not a probe object.")
  
  # Find center points on phantom boundary
  centers <- phantom$opposingPoints(point, direction)
  # Convert to world coordinates
  centers[[1]] <- phantom$local2world(centers[[1]])
  centers[[2]] <- phantom$local2world(centers[[2]])
  
  # Position the transmitting probe
  v <- centers[[1]]-centers[[2]]
  t1 <- diag(1,4,4)
  # rotate probe's z axis until it aligns with v
  t1[1:3, 1:3] <- rotateU2V(c(0,0,1), v)
  # rotate about v until transformed x axis is perpendicular to world z
  #   find the angle transformed x axis makes with world z
  theta <- acos(t1[3,1])
  #   rotate about v by theta-pi/2 so that transformed x is parallel world x/y plane.
  t1[1:3, 1:3] <- rotationMatrix(v, theta-pi/2) %*% t1[1:3, 1:3]
  # translate to centers[[2]]
  t1 <- affineTransform(c(0,0,0), 0, centers[[2]]) %*% t1 
  xmitr$moveTo(t1)
  # Position the receiving probe
  t2 <- diag(1,4,4)
  t2[1:3, 1:3] <- rotateU2V(c(0,0,1), -v) # rotate z axis in the direction -v
  #   find the angle which the transformed x axis makes with world z
  theta <- acos(t2[3,1])
  t2[1:3, 1:3] <- rotationMatrix(-v, theta-pi/2) %*% t2[1:3, 1:3]
  t2 <- affineTransform(c(0,0,0), 0, centers[[1]]) %*% t2 # translate
  rcvr$moveTo(t2)
  
  # Transform transducer array positions from internal probe coordinates to
  # internal phantom coordinates.
  # First convert from probe to world coordinates
  transmitters <- sapply(xmitr$transmitters, function(x)xmitr$local2world(c(x, 0, 0)))
  receivers <- sapply(rcvr$receivers, function(x)rcvr$local2world(c(x, 0, 0)))
  # Next convert from world to phantom
  transmitters <- sapply(1:(ncol(transmitters)), function(k){phantom$world2local(as.vector(transmitters[,k]))})
  receivers <- sapply(1:(ncol(receivers)), function(k){phantom$world2local(as.vector(receivers[,k]))})
  
  # Return a list of transmitter and receiver coordinates to support a scan,
  # and a z coordinate to support plotting an axial section.
  list(transmitters=transmitters, receivers=receivers, z=point[3])
}

plotSectionAndArrays <- function(phantom, transmitters, receivers, z, by=4){
  phantom$plotAxialSection(z)
  idx <- seq(1, ncol(transmitters), by=by)
  points(transmitters[1,idx], transmitters[2,idx], pch=19, col="blue")
  pos <- c(1, 4)[which.max(apply(transmitters[1:2,], 1, sd))]
  text(transmitters[1,idx], transmitters[2,idx], labels=as.character(idx), pos=pos, col="blue" )
  idx <- seq(1, ncol(receivers), by=by)
  points(receivers[1,idx], receivers[2,idx], pch=19, col="red")
  pos <- c(3, 2)[which.max(apply(transmitters[1:2,], 1, sd))]
  text(receivers[1,idx], receivers[2,idx], labels=as.character(idx), pos=pos, col="red" )
  legend("topright", c("transmitter", "receiver"), pch=19, col=c("blue", "red"))
}

doScan <- function(phantom, transmitters, receivers){
  ans <- matrix(0, ncol(transmitters), ncol(receivers))
  for(i in 1:ncol(transmitters)){
    for(j in 1:ncol(receivers)){
      ans[i,j] <- phantom$timeOfFlight(as.vector(transmitters[,i]), as.vector(receivers[,j]), warn=FALSE)
    }
  }
  ans
}

plotScan <- function(scan){
  image(1:nrow(scan), 1:ncol(scan), scan, asp=1, col=rainbow(100), xlab="Transmitter i", ylab="Receiver j", main="Time of flight vs path i,j")
  temp <- round(seq(min(scan), max(scan), length.out = 10), 4)
  legend('topright', as.character(temp), fill=rainbow(10), title="ToF (msec)")
}