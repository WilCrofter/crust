#' Returns a "calibration link" between the reference frames of two objects.
#' A calibration link represents a putative, possibly erroneous, relationship
#' between the reference frames, not their true relationship in the physical world.
#' @param base2relative, an affine transform which converts coordinates in base frame to those in the relative.
newLink <- function(base_object, relative_object, base2relative){
  # Argument checks
  if(!is.numeric(base2relative) | !isTRUE(all.equal(dim(base2relative), c(4,4))) | 
       !isTRUE(all.equal(base2relative[1:3, 1:3] %*% t(base2relative[1:3, 1:3]), diag(1,3,3) )) |
       !isTRUE(all.equal(base2relative[4, 1:3], c(0,0,0)))){
    stop("Argument 'base2relative' is not an affine transform")
  }
  if(!is(base_object, "rigid body"))stop("Base object is not a rigid body")
  if(!is(relative_object, "rigid body"))stop("Relative object is not a rigid body")
  if(identical(relative_object, base_object))stop("Base and relative objects cannot be identical.")
    
  # Inverse transform of base2relative
  relative2base <- invAffine(base2relative)
  
  # Given a point, u, in base coordinates, convert to relative coordinates
  base2Relative <- function(u){
    as.vector(base2relative %*% c(u,1))[1:3]
  }
  
  # Given a point, u, in relative coordinates, convert to base coordinates
  relative2Base <- function(u){
    as.vector(relative2base %*% c(u,1))[1:3]
  }
  
  # Replace the current relationship between base and relative coordinate systems
  # with that given by the argument
  moveRelative <- function(transform){
    this <- environment(moveRelative)
    this$relative2base <- invAffine(transform)
    this$base2relative <- transform
  }
  
  # Adjust the current relationship between base and coordinate systems by
  # left-multiplying base2relative by the given transform and right-multiplysing
  # its inverse to relative2base.
  adjustRelative <- function(transform){
    this <- environment(adjustRelative)
    this$relative2base <- this$relative2base %*% invAffine(transform)
    this$base2relative <- transform %*% this$base2relative
  }
  
  # Return a reference to the current environment, after adding
  # "link" to its class attribute.
  e <- environment(relative2Base)
  class(e) <- c("link", class(e))
  e
}