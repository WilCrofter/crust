# Given horizontal grid lines crossing the x (vertical) axis at points x=hg[i], 
# vertical grid lines crossing the z (horizontal) axis at points z=vg[j], 
# and the coefficients of a line, b, and m  where x = b + m*z, 
# derive the sequence of points at which the line crosses a grid line.
findPixelCrossings <- function(hg, vg, b, m){
  # The line will hit every vertical grid line; z is the horizontal coordinate,
  # so the associated grid lines are vertical. We just solve for the corresponding
  # x (vertical) coordinate.
  crossings <- data.frame(z=vg, x=b+m*vg)
  # If m == 0, the line will hit every vertical grid line but will not enter or exit
  # pixels in the vertical direction. If m is not zero, we solve for the points
  # at which the line hits horizontal grid lines.
  if(m != 0){
    temp <- data.frame(z =(hg-b)/m , x=hg)
    # Not every z coordinate will be within range.
    idx <- temp[,"z"] >= min(hg) & temp[,"z"] <= max(hg)
    # Append those which are to the points of intersection
    crossings <- rbind(crossings, temp[idx,])
  }
  # Order the crossings by the z (horizontal) coordinates
  crossings <- crossings[order(crossings[,"z"]), ]
}

# Find the transmitter-receiver pairs whose US signal will cross a given pixel.
# cornerA, cornerB -- (x, z) coordinates of opposing pixel corners, e.g., lower left and upper right.
# xmit -- the x (vertical) coordinates of the transmitter array, xmit[1], xmit[2], etc.
# zxmit -- the z coordinate of the transmitter array.
# rcv -- the x coordinates of the receiver array.
# zrcv -- the z coordinate of the receiver array.
findLinesCrossingPixel <- function(cornerA, cornerB, xmit, zxmit, rcv, zrcv){
  # Find coordinates of grid lines which bound the pixel
  zgrids <- range(cornerA[1], cornerB[1])
  xgrids <- range(cornerA[2], cornerB[2])
  # A line segment drawn between transmitter and receiver is a convex
  # combination of the two endpoints: (1-lambda)*v1 + lambda*v2,
  # where lambda is between 0 and 1, and v1, v2 are the endpoint coordinate
  # vectors. Such a lambda will exist if and only if the the endpoint coordinates
  # satisfy certain inequalities. A range of possible lambdas can be derived using
  # inequality constraints on the z coordinate.
  lrange <- range((zgrids - zxmit)/(zrcv-zxmit))
  # Given lrange and the x coordinates of a transmitter and receiver, find the range
  # of x coordinates which can be obtained with lambda in lrange.
  x_mins <- apply(cbind((xgrids[1]-(1-lrange[1])*xmit)/lrange[1], 
                           (xgrids[1]-(1-lrange[2])*xmit)/lrange[2]),
                     1, min)
  x_maxs <- apply(cbind((xgrids[2]-(1-lrange[1])*xmit)/lrange[1], 
                           (xgrids[2]-(1-lrange[2])*xmit)/lrange[2]),
                     1, max)
  
  # Return a list of receivers for each transmitter. Strict inequalities insure
  # that a line actually enters the pixel rather than just grazing a boundary.
  lapply(1:length(xmit), function(n)rcv[rcv > x_mins[n] & rcv < x_maxs[n]])
}


# compute margins which make plot window square
adj_margins <- function(){
  fin <- par("fin") # (window in inches)
  mai <- par("mai") # (margins in inches)
  # We want fin[1]-(mai[2]+mai[4]) == fin[2] - (mai[1]+mai[3])
  w <- fin[1]-(mai[2]+mai[4]) # plot window width
  h <- fin[2]-(mai[1]+mai[3]) # plot window height
  # We'll preserve the smallest of these two dimensions
  if(w > h){
    delta <- (w-h)/2
    mai[c(2,4)] <- mai[c(2,4)] + delta
  } else {
    delta <- (h-w)/2
    mai[c(1,3)] <- mai[c(1,3)] + delta
  }
  mai
}

# Returns a 2x51 matrix representing the vertices of a
# polygon with approximates the ellipse formed by 
# applying the 2x2 matrix, M, to a circle of radius
# r which is centered at the origin.
# If result <- ellipse(r, M), use
# polygon(x+result[1,], y+result[2,]) to plot the
# the ellipse centered at x, y.
ellipse <- function(r, M=matrix(c(1,0,0,1), 2, 2)){
  theta <- seq(0, 2*pi, length.out = 51)
  temp <- r*cbind(cos(theta), sin(theta)) %*% M
  colnames(temp) <- c("x", "y")
  return(temp)
}

# Given an axis of rotation, u, and an angle, theta, return a 3x3 matrix
# representing rotation about u by theta, assuming the right-hand rule.
rotationMatrix <- function(u, theta){
  if(!is.numeric(u) | !is.vector(u) | 
       length(u) != 3)stop("axis of rotation must be a 3-vector")
  if(!is.numeric(theta) | length(theta) != 1)stop("theta must be a number")
  normu <- sqrt(sum(u^2))
  if(!isTRUE(all.equal(theta,0)) & 
       isTRUE(all.equal(normu, 0)))stop("axis of rotation cannot be zero unless theta is zero")
  if(normu > 0) u <- u/normu
  N <- matrix(c(0, u[3], -u[2], -u[3], 0, u[1], u[2], -u[1], 0), 3, 3)
  diag(1, 3, 3) + sin(theta)*N + (1-cos(theta))*(N %*% N)
}

# Return an affine transform (member of SE(3)) defined by an axis, u,
# an angle of rotation, theta, and a translation.
affineTransform <- function(u, theta, translation){
  ans <- diag(1, 4, 4)
  ans[1:3, 1:3] <- rotationMatrix(u, theta)
  ans[1:3, 4] <- translation
  ans
}

# Return TRUE if and only if the argument is an affine transform
isAffine <- function(transform){
  if(!is.numeric(transform))return(FALSE)
  if(!is.matrix(transform))return(FALSE)
  if(!isTRUE(all.equal(dim(transform), c(4,4))))return(FALSE)
  if(!isTRUE(all.equal(transform[1:3,1:3] %*% t(transform[1:3, 1:3]), diag(1,3,3))))return(FALSE)
  if(!isTRUE(all.equal(transform[4,1:4], c(0,0,0,1))))return(FALSE)
  TRUE
}

# Return the inverse an affine transform
invAffine <- function(transform){
  if(!isAffine(transform))stop("Attempt to take affine inverse of non-affine transform.")
  # Inverse transform of base2relative
  inv <- diag(1, 4, 4)
  inv[1:3,1:3] <- t(transform[1:3, 1:3])
  inv[1:3,4] <- -inv[1:3,1:3] %*% transform[1:3, 4]
  inv
}


# Cross product uXv
crossProduct <- function(u, v){
  if(!is.numeric(u) | !is.numeric(v) | length(u) != 3 | length(v) != 3)
    stop("u and v must be numeric 3-vectors")
  u[c(2,3,1)]*v[c(3,1,2)] - u[c(3,1,2)]*v[c(2,3,1)]
}

# Return the rotation matrix which rotates u into the direction of v
# about the axis uXv.
rotateU2V <- function(u, v){
  if(!is.numeric(u) | !is.numeric(v) | length(u) != 3 | length(v) != 3)
    stop("u and v must be numeric 3-vectors")
  normu <- sqrt(sum(u^2))
  normv <- sqrt(sum(v^2))
  if(isTRUE(all.equal(normu, 0)))stop("u must be nonzero")
  if(isTRUE(all.equal(normv, 0)))stop("v must be nonzero")
  u <- u/normu
  v <- v/normv
  theta <- acos(sum(u*v))
  rotationMatrix(crossProduct(u,v), theta)
}

# Generates n points uniformly distributed on the unit sphere
# or half sphere (half=TRUE).
# Reference: http://mathworld.wolfram.com/SpherePointPicking.html
rusphere <- function(n, half=FALSE){
  theta <- runif(n, 0, 2*pi)
  cosu <- runif(n, half-1, 1) # 'half' is automatically cast from logical to numerical
  sinu <- sqrt(1-cosu^2)
  cbind(x=sinu*cos(theta), y=sinu*sin(theta), z=cosu)
}

# Returns an affine transform representing random alignment error. Translation errors are
# normally distributed. Rotation errors are described by a slight rotation about a 
# random axis. The angle of rotation is generated by a symmetric beta distribution.
# Default standard deviations follow the reference, "we can expect a translational error 
# less than 4 mm; assuming the rotational error is less than 0.5 degrees, the expected 
# error in each axis will be around 3 mm"
# Reference: Aalamifar et. al., Co-robotic ultrasound tomography: dual arm setup and 
# error analysis
alignmentError <- function(sd_translation=3/2, sd_angle=pi/360){
  alpha <- (pi^2/sd_angle^2 - 1)/2
  affineTransform(as.vector(rusphere(1)), 
                  pi*(rbeta(1, alpha, alpha)-.5), 
                  rnorm(3,sd = sd_translation))
}

