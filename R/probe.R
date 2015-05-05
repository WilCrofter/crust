#' Creates a probe of n transducers equally spaced along the x axis
#' in the local reference frame, and at a physical position and
#' orientation given by frame_in_world. 
newProbe <- function(n, spacing, alignment_in_world){
  n <- as.integer(round(n))
  if(length(n) != 1)stop("n must be an integer")
  if(!is.numeric(spacing) | length(spacing) != 1)stop("spacing must be a real number")
  if(!isAffine(alignment_in_world))stop("position and orientation in world is not an affine transform")
  
  receivers <- spacing*(seq(1, n, by=1)-n/2)
  transmitters <- spacing*(seq(1.5, n-.5, by=1)-(n-1)/2)
  
  # Change the real position and orientation of this object in the physical world.
  # CAUTION: relative positions as they appear in links will not be affected.
  moveTo <- function(transform){
    if(!isAffine(transform))stop("Argument is not an affine transform")
    this <- environment(moveTo)
    this$alignment_in_world <- transform
    invisible()
  }
  
  # Adjust the real position and orientation of this object in the physical world.
  # I.e., left multiply its current position and orientation by the give transform.
  # CAUTION: relative positions as they appear in links will not be affected.
  adjustBy <- function(transform){
    if(!isAffine(transform))stop("Argument is not an affine transform")
    this <- environment(adjustBy)
    this$alignment_in_world <- transform %*% this$alignment_in_world
    invisible()
  }
  
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
  
  # Return a reference to the current environment after adding
  # "rigid body" and "box phantom" class attributes
  e <- environment(world2local)
  class(e) <- c("probe", "rigid body", class(e))
  e
  
}