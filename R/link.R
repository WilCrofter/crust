#' Returns a "link" between the reference frames of two objects,
#' given an affine transform, base2relative, which converts coordinates in the
#' base object to those in the relative object.
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
  relative2base <- diag(1, 4, 4)
  relative2base[1:3,1:3] <- t(base2relative[1:3, 1:3])
  relative2base[1:3,4] <- -relative2base[1:3,1:3] %*% base2relative[1:3, 4]
  
  # Given a point, u, in base coordinates, convert to relative coordinates
  base2Relative <- function(u){
    as.vector(base2relative %*% c(u,1))[1:3]
  }
  
  # Given a point, u, in relative coordinates, convert to base coordinates
  relative2Base <- function(u){
    as.vector(relative2base %*% c(u,1))[1:3]
  }
  
  # Return a reference to the current environment, after adding
  # "link" to its class attribute.
  e <- environment(relative2Base)
  class(e) <- c("link", class(e))
  e
}