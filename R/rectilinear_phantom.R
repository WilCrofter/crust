#' Returns a rectilinear phantom composed of silicon with two embedded spheres containing
#' alcohol and water, respectively. The phantom's left, front, bottom corner is the
#' the origin of the local reference frame.
#' @param dimensions --a 3-vector of width, depth, and height, the x, y, z dimensions respectively
#' @param ctr_alcohol --a 3-vector specifying the center of the alcohol component
#' @param r_alcohol --the radius of the alcohol component
#' @param ctr_water --a 3-vector specifying the center of the water component
#' @param r_water --the radius of the water component
makeRectilinearPhantom <- function(dimensions, 
                                   ctr_alcohol, r_alcohol, ctr_water, r_water){
  ### Argument checks
  # Ensure all arguments are numeric
  for(a in list(dimensions, ctr_alcohol, r_alcohol, ctr_water, r_water)){
    if(!is.numeric(a))stop(paste("Argument ", a, "is not numeric."))
  }
  # Ensure that dimensions and centers are 3 vectors
  if(!(length(dimensions)==3))stop("dimensions must be a 3 vector")
  if(!(length(ctr_alcohol)==3))stop("ctr_alcohol must be a 3 vector")
  if(!(length(ctr_water)==3))stop("ctr_water must be a 3 vector")
  # Ensure that dimensions are positive
  if(!(sum(dimensions > 0)))stop("Phantom must have positive width, depth, and height.")
  # Ensure that radii are positive numbers
  if(!isTRUE(r_alcohol > 0) | length(r_alcohol) != 1)stop("r_alcohol is not a positive number")
  if(!isTRUE(r_water > 0) | length(r_water) != 1)stop("r_water is not a positive number")
  ## Ensure that the spheres are contained within the phantom body.
  if(min(ctr_alcohol-r_alcohol) < 0 |
       max(ctr_alcohol-dimensions+r_alcohol) > 0)stop("Alcohol sphere is out of bounds.")
  if(min(ctr_water-r_water) < 0 |
       max(ctr_water-dimensions+r_water) > 0)stop("Water sphere is out of bounds.")
  
  # Is a given point within the phantom?
  isContained <- function(point){
    sum(point >= 0 & point <= dimensions) == length(point)
  }
  
  # Given a point and a vector indicating direction, return a list of two
  # diametrically opposite points on the phantom's surface
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
  
  
  
  
}