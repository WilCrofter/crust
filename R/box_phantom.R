#' Returns a box (rectilinear) phantom composed of silicon with two embedded spheres containing
#' alcohol and water, respectively. The phantom's left, front, bottom corner is the
#' the origin of the local reference frame.
#' @param dimensions --a 3-vector of width, depth, and height, the x, y, z dimensions respectively
#' @param ctr_alcohol --a 3-vector specifying the center of the alcohol component
#' @param r_alcohol --the radius of the alcohol component
#' @param ctr_water --a 3-vector specifying the center of the water component
#' @param r_water --the radius of the water component
newBoxPhantom <- function(dimensions, ctr_alcohol, r_alcohol, ctr_water, r_water,
                          alignment_in_world=diag(1, 4, 4)){
  ### Speed of sound
  c_silicone <- 1000.0 # mm/msec
  c_alcohol <- 1180.0 # mm/msec (ethyl alcohol)
  c_water <- 1480.0 # mm/msec
    
  ### Argument checks
  # Ensure all arguments are numeric
  for(a in list(dimensions, ctr_alcohol, r_alcohol, ctr_water, r_water)){
    if(!is.numeric(a))stop(paste("Argument ", a, "is not numeric."))
  }
  rm(a)
  # Ensure that dimensions and centers are 3 vectors
  if(!(length(dimensions)==3))stop("dimensions must be a 3 vector")
  if(!(length(ctr_alcohol)==3))stop("ctr_alcohol must be a 3 vector")
  if(!(length(ctr_water)==3))stop("ctr_water must be a 3 vector")
  # Ensure that dimensions are positive
  if(!(sum(dimensions > 0)))stop("Phantom must have positive width, depth, and height.")
  # Ensure that radii are positive numbers
  if(!isTRUE(r_alcohol > 0) | length(r_alcohol) != 1)stop("r_alcohol is not a positive number")
  if(!isTRUE(r_water > 0) | length(r_water) != 1)stop("r_water is not a positive number")
  # Ensure that the spheres are contained within the phantom body.
  if(min(ctr_alcohol-r_alcohol) < 0 |
       max(ctr_alcohol-dimensions+r_alcohol) > 0)stop("Alcohol sphere is out of bounds.")
  if(min(ctr_water-r_water) < 0 |
       max(ctr_water-dimensions+r_water) > 0)stop("Water sphere is out of bounds.")
  # Ensure that the spheres don't overlap
  if(r_alcohol+r_water > sqrt(sum((ctr_alcohol-ctr_water)^2)))stop("Water and alcohol spheres overlap.")

  # Convert a 3-vector OR an affine transform from world coordinates
  # to local coordinates.
  world2local <- function(u){
    if(isAffine(u)){
      return(invAffine(alignment_in_world) %*% u)
    } else if(is.numeric(u) & length(u) == 3){
      return(as.vector(invAffine(alignment_in_world) %*% c(u,1))[1:3])
    } else {
      stop("Argument is neither a 3-vector nor an affine transform.")
    }
  }
  
  # Convert a 3-vector OR an affine transform from local coordinates
  # to world coordinates.
  local2world <- function(u){
    if(isAffine(u)){
      return(alignment_in_world %*% u)
    } else if(is.numeric(u) & length(u) == 3){
      return(as.vector(alignment_in_world %*% c(u,1))[1:3])
    } else {
      stop("Argument is neither a 3-vector nor an affine transform.")
    }
  }
  
  # Within tolerance, is a given point inside or on the boundary of the phantom?
  isContained <- function(u){
    sum(sapply(1:3, function(k){
      (isTRUE(all.equal(u[k],0)) | u[k] > 0) & 
        (isTRUE(all.equal(u[k], dimensions[k])) | u[k] < dimensions[k])
    })) == 3
  }
  
  # Given a point and a vector indicating direction, return a list of two
  # diametrically opposite points on the phantom's surface which fall on
  # the line L(m) = point + m*direction.
  opposingPoints <- function(point, direction){
    if(!is.numeric(point) | !(length(point)==3))stop("point must be numeric of length 3")
    if(!is.numeric(direction) | !(length(direction)==3))stop("direction must be numeric of length 3")
    if(isTRUE(all.equal(direction,c(0,0,0))))stop("direction cannot be zero")
    # Solve for multiples of direction which cross the faces
    r0 <- -point/direction
    r1 <- (dimensions - point)/direction
    # Order them according to absolute value
    r0 <- r0[order(abs(r0))]
    r1 <- r1[order(abs(r1))]
    # Check the associated points for containment
    i0 <- sapply(r0, function(r)isContained(point + r*direction))
    i1 <- sapply(r1, function(r)isContained(point + r*direction))
    # Pick the first contained point
    pt1 <- point + r0[i0][1]*direction
    pt2 <- point + r1[i1][1]*direction
    return(list(pt1, pt2))
  }
  
  ## Find the length of the intersection of a line segment with a sphere.
  lengthThruSphere <- function(u, v, ctr, radius){
    # Translate to origin
    u <- u - ctr
    v <- v - ctr
    # Solve for x in [0,1] such that (1-x)u + xv
    # is on the surface of the sphere
    pnom <- c(sum(u^2) - radius^2,
              2*sum(u*(v-u)),
              sum((v-u)^2))
    roots <- polyroot(pnom)
    # polyroot returns numbers which are nominally complex
    # even if their imaginary parts are zero.
    x <- Re(roots)
    if(!isTRUE(all.equal(Im(roots), c(0,0)))){
      # Roots have non-zero imaginary parts (within tolerance)
      # hence line does not intersect sphere.
      return(0)
    } else {
      # If either root is outside of the interval [0,1], then
      # either u or v is inside the sphere. Truncate roots
      # to the unit interval.
      x <- pmin(pmax(x, 0), 1)
      # Calculate the length of the line segment between
      # (1-x[1])*u + x[1]*v and (1-x[2])*u + x[2]*v
      delta <- x[1]-x[2]
      segment <- delta*(v-u)
      return(sqrt(sum(segment^2)))
    }
  }
  
  #' Plot an axial (horizontal) section of the phantom at z=k
  #' TODO: generalize to a planar section of any sort
  plotAxialSection <- function(k){
    # Argument checks
    if(!is.numeric(k) | length(k) != 1)stop("k must be a number")
    if(!(k>=0) | !(k <= dimensions[3]))stop("section is out of bounds")
    plot(c(0,dimensions[1]), c(0,dimensions[2]), type='n', xlab="x", ylab="y", asp=1,
         main=paste("Axial section of box phantom at z =", k))
    polygon(c(0, dimensions[1], dimensions[1], 0, 0), 
            c(0, 0, dimensions[2], dimensions[2], 0), col="lightyellow", lwd=3)
    if(abs(k-ctr_alcohol[3]) < r_alcohol){
      r <- sqrt(r_alcohol^2 - (ctr_alcohol[3]-k)^2)
      theta <- (pi/25)*(0:50)
      polygon(ctr_alcohol[1]+r*cos(theta), ctr_alcohol[2]+r*sin(theta), col="pink", lwd=3)
    }
    if(abs(k-ctr_water[3]) < r_water){
      r <- sqrt(r_water^2 - (ctr_water[3]-k)^2)
      theta <- (pi/25)*(0:50)
      polygon(ctr_water[1]+r*cos(theta), ctr_water[2]+r*sin(theta), col="lightblue", lwd=3)
    }
    legend('topleft', c("silicone", "alcohol", "water"), bg="white", fill=c("lightyellow", "pink", "lightblue"))
  }
  
  #' Given a line segment defined by two endpoints, u and v, in the phantom or
  #' on its boundary, determine the time of flight between them.
  timeOfFlight <- function(u, v){
    # Argument checks
    if(!is.numeric(u) | length(u) != 3)stop("Argument u must be a numeric 3-vector.")
    if(!is.numeric(v) | length(v) != 3)stop("Argument v must be a numeric 3-vector.")
    if(!isContained(u))stop("Argument u is outside of the phantom body.")
    if(!isContained(v))stop("Argument v is outside of the phantom body.")
    l_total <- sqrt(sum((u-v)^2))
    l_alcohol <- lengthThruSphere(u, v, ctr_alcohol, r_alcohol)
    l_water <- lengthThruSphere(u, v, ctr_water, r_water)
    (l_total-l_alcohol-l_water)/c_silicone + l_alcohol/c_alcohol + l_water/c_water 
  }

  # Return a reference to the current environment after adding
  # "rigid body" and "box phantom" class attributes
  e <- environment(timeOfFlight)
  class(e) <- c("box phantom", "rigid body", class(e))
  e
}