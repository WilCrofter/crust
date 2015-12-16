####################################################################################
# Functions to position a rectangle over a pixelated image and to calculate
# transducer coordinates based on the rectangle. The pixels of the image in question
# must be labeled by tissue type, (as opposed to slowness or speed,) with pixels 
# external to tissue labeled 0.
#
# The main functions implemented below are displayImage, minRect (which computes
# a rectangle,) and transducerPos (which computes transducer positions from
# the output of minRect.) A rectangle can be displayed over an image using the
# standard R function, polygon.
#
# These functions depend on R/utilities.R. Whole thigh images are stored as csv
# files but must be converted to matrices when loaded. For instance,
#               img <- as.matrix(read.csv("data/whole_thigh_32.csv"))
####################################################################################

# Display an image, given its grid size in mm. In the case of the whole thigh images
# for which these functions were written, grid size can be estimated from the
# file names. For file data/whole_thigh_K.csv, a grid size of 
# 30/(K-.5) + 30/(K-1.5) would be reasonable. (Don't ask. Anyway, grid size is 
# not critical. It merely affects the scale of the phantom. A typical thigh might 
# be 20 cm in diameter.)
displayImage <- function(img, gridsize, xlab="x (mm)", ylab="y (mm)", asp=1, ...){
  image(gridsize*(1:nrow(img)-.5), gridsize*(1:ncol(img)-.5), img, 
        xlab=xlab, ylab=ylab, asp=asp, ...)
}

# Given a grid size, convert from x,y coordinates, c(x,y), to row-column coordinates
# assuming the pixel, img[1,1], has corners c(0,0) and c(gridsize, gridsize)
xy2pixel <- function(xycoords, gridsize){
  # The floor and ceiling of an integer are equal, hence the use of floor + 1
  1+floor(xycoords/gridsize)
}

# Given grid size and the indices, c(i,j), of a pixel, return the x,y
# coordinates of its midpoint.
pixel2xy <- function(ijpix, gridsize){
  (ijpix-0.5)*gridsize
}

# Find the points of a grid crossed by a line segment whose two endpoints are u and v
# both of which are 2-tuples. The pixel grid is defined by a gridsize and
# assumed parallel to the x and y axes. Crossings will be in order from u to v.
gridCrossings <- function(u, v, gridsize, img, eps = (.Machine$double.eps)^.75){
  kx <- seq(0, 1+nrow(img))*gridsize
  ky <- seq(0, 1+ncol(img))*gridsize
  # wCrossings requires vectors with 3 components:
  lambdas <- sort(unique(c(wCrossings(c(u,0), c(v,0), c(1,0,0), kx), 
                           wCrossings(c(u,0), c(v,0), c(0,1,0), ky))))
  # Since lambdas are sorted in ascending order, crossings will be
  # in order from u to v. Eliminate any which are too close together.
  idx <- which(diff(lambdas) < eps)
  if(length(idx)>0)lambdas <- lambdas[-(1+idx)]
  crossings <- sapply(lambdas, function(lambda)u+lambda*(v-u))
}

# Given coordinates of grid crossings, find the pixels crossed.
pixelsCrossed <- function(crossings, gridsize){
  # The average of two consecutive crossings should be in the interior
  # of the pixel which they bound
  midpoints <- (crossings[,-ncol(crossings)] + crossings[,-1])/2
  apply(midpoints, 2, function(u)xy2pixel(u,gridsize))
}

# Given a 2xn array of pixel indices, return only the indices
# which satisfy the given function
pixelsSatisfying <- function(indices, img, fct){
  indices[,which(apply(indices, 2, function(idx)fct(img[idx[1], idx[2]])))]
}

# Given a 2xn array of pixel indices, cull those out of bounds with respect
# to the dimensions of img
cullIndices <- function(indices, img){
  idx <- indices[1,] > 0 & indices[1,] <= nrow(img) & 
    indices[2,] > 0 & indices[2, ] <= ncol(img)
  indices[,idx]
}

# Given two endpoints, u and v, of a line segment, return the indices of the
# nonzero pixels which it crosses.
nonzeroPix <- function(u, v, img, gridsize){
  # Grid crossings will be returned in order from u to v
  crossings <- gridCrossings(u, v, gridsize, img)
  # In-bound pixels crossed by the line in order from u to v
  pix <- cullIndices(pixelsCrossed(crossings, gridsize), img)
  # Of these, the pixels with nonzero entries
  pixelsSatisfying(pix, img, function(x)x!=0)
}

# Given a center and direction, determine a point of the form
# center + lambda*direction on the boundary of the image,
# where lambda >=0
boundaryPt <- function(center, direction, img, gridsize){
  ur <- gridsize*c(nrow(img),ncol(img))
  ll <- c(0,0)
  lambda <- c( (ur-center)/direction, (ll-center)/direction )
  center + min(lambda[lambda >= 0])*direction
}

# Find a minimal rectangle of specified position, orientation, and height which encloses
# a portion of the non-zero pixels in an image. The gridsize, center coordinates, and 
# height must be in the same units, usually mm.
minRect <- function(img, gridsize, center, angle_in_radians, height){
  # Form unit vector at the given angle.
  direction <- c(cos(angle_in_radians), sin(angle_in_radians))
  # Find u and v on the line determined by center and direction, which
  # include the entire image.
  u <- boundaryPt(center, -direction, img, gridsize)
  v <- boundaryPt(center, +direction, img, gridsize)
  # Indices of nonzero pixels crossed by the segment
  nzpix <- nonzeroPix(u, v, img, gridsize)
  # That the line segment should intersect tissue is a programming convenience.
  # It's not strictly necessary, but would require more code than it's worth.
  if(length(nzpix) == 0)stop("Given line segment does not intersect tissue.")
  # The endpoints of nzpix approximately bound the tissue.
  aprxbounds <- nzpix[,c(1,ncol(nzpix))]
  # Find their midpoints in x,y coordinates
  midpts <- apply(aprxbounds, 2, function(w)pixel2xy(w,gridsize))
  # Find the closest points to them on the line determined by u and v.
  lambdaA <- t(v-u)%*%(midpts[,1]-u)/sum((v-u)^2)
  lambdaB <- t(v-u)%*%(midpts[,2]-u)/sum((v-u)^2)
  ptA <- (1-lambdaA)*u + lambdaA*v
  ptB <- (1-lambdaB)*u + lambdaB*v
  # Find their average and a unit vector pointing from their average toward
  # ptA. This will be useful for moving both points outward.
  avg <- (ptA+ptB)/2
  outA <- (ptA-avg)/sqrt(sum((ptA-avg)^2))
  # Find a unit vector at right angles (in a counter-clockwise sense) to the 
  # line determined by u and v (equivalently, by ptA and ptB.)
  perp <- c(cos(angle_in_radians + pi/2), sin(angle_in_radians + pi/2)) 
  # Form the endpoints of a line segment of the given height, parallel to perp
  # and centered at the origin.
  ep1 <- +height*perp/2
  ep2 <- -height*perp/2
  # Translate the endpoint centers to ptA then move ptA outward 
  # until the associated line segment fails to intersect tissue.
  cnr1 <- ep1+ptA
  cnr2 <- ep2+ptA
  while(length(nonzeroPix(cnr1+gridsize*outA, cnr2+gridsize*outA, img, gridsize)) >0){
    cnr1 <- cnr1+gridsize*outA
    cnr2 <- cnr2+gridsize*outA
  }
  # Similarly for ptB
  cnr4 <- ep1+ptB
  cnr3 <- ep2+ptB
  while(length(nonzeroPix(cnr3-gridsize*outA, cnr4-gridsize*outA, img, gridsize)) >0){
    cnr3 <- cnr3-gridsize*outA
    cnr4 <- cnr4-gridsize*outA
  }
  # Return the derived rectangle's 4 corners as a 4x2 array.
  rbind(cnr1, cnr2, cnr3, cnr4)
}

# Given a number, n, of equally spaced transducers and two points, u and v, representing
# the two endpoints of a probe, return transducer positions as an n by 2 matrix. 
# Transducers are placed in the centers of the n intervals of equal length which 
# partition the line segment between u and v.
#        u--1--|--2--|...|--n--v
transducerPos <- function(n, u, v){
  # First transducer is half a space from u.
  # Subsequent transducers are full spaces from one another.
  t(sapply(1:n, function(k)u+(k-.5)*(v-u)/n))
}

# Given a rectangle, r, as calculated by minRect, above, and a number, n, of transducers,
# return a list of 2 nx2 matrices of probe positions, one for each of two probes,
# one positioned between corners 1 and 2 of the rectangle, the other between corners
# 3 and 4. The probes will be numbered in a counter-clockwise direction, the first
# probe numbering from corner 1 to corner 2, the second from corner 3 to corner 4.
# probe numbering is easily reversed, e.g., probe <- probe[nrow(probe):1,], if
# desired.
probePos <- function(n, r){
  list(probe1=transducerPos(n, r[1,], r[2, ]),
       probe2=transducerPos(n, r[3,], r[4, ]))
}