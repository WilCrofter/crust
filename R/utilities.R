# Given horizontal grid lines crossing the x (vertical) axis at, hg[i] = i*dx, 
# vertical grid lines crossing the z (horizontal) axis at vg[j] = j*dz, 
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