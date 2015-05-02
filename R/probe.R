#' Creates a probe of n transducers equally spaced along the x axis
#' in the local reference frame. 
newProbe <- function(n, spacing){
  n <- as.integer(round(n))
  if(length(n) != 1)stop("n must be an integer")
  if(!is.numeric(spacing) | length(spacing) != 1)stop("spacing must be a real number")
  
  recievers <- spacing*seq(1, n, by=1)
  transmitters <- spacing*seq(1.5, n-.5, by=1)
  links <- list()
  
  addLink <- function(link){
    this <- environment(addLink)
    if(!is(link, "link") | 
         (!identical(link$base_object, this) & 
            !identical(link$relative_object, this))){
      warning("Argument to addLink is not a link with this object. It will not be added.")
      return(FALSE)
    } else {
      links <- c(links, link)
      return(TRUE)
    }
  }
  
  removeLink <- function(link){
    if(!is(link, "link")){
      warning("Argument to removeLink is not a link. Attempted removal failed.")
      return(FALSE)
    } else {
      idx <- sapply(links, function(x)identical(x,link))
      links <- links[!idx]
      return(TRUE)
    }
  }
  
  # Apply the given affine transform to the local reference frame.
  # In practice this means applying it to every link which connects
  # this body to another, thus changing its relative position with
  # respect to the other bodies.
  transformOrigin <- function(transform){
    this <- environment(transformOrigin)
    # Inverse of the given transform
    inv <- diag(1, 4, 4)
    inv[1:3,1:3] <- t(transform[1:3, 1:3])
    inv[1:3,4] <- -inv[1:3,1:3] %*% transform[1:3, 4]
    
    for(link in links){
      if(identical(link$base_object), this){
        # This body is the base object hence
        link$base2relative <- transform %*% b2r
        link$relative2base <- link$relative2base %*% inv
      } else if(identical(link$relative_object), this){
        # This body is the relative object, hence
        link$base2relative <- inv %*% link$base2relative
        link$relative2base <- link$relative2base %*% transform
      }
    }
  }
  
}