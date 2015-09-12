

# Use ICM to remove the noise from the given image.
# * covar is the known covariance of the Gaussian noise.
# * max_diff is the maximum contribution to the potential
#   of the difference between two neighbouring pixel values.
# * weight_diff is the weighting attached to the component of the potential
#   due to the difference between two neighbouring pixel values.
# * iterations is the number of iterations to perform.

restore_image <- function(src, covar, max_diff, weight_diff, iterations, vals){
  
  # Maintain two buffer images.
  # In alternate iterations, one will be the
  # source image, the other the destination.
  buffer <- array(0, c(dim(src),2))
  buffer[,,1] <- src
  s <- 2
  d <- 1
  
  # This value is guaranteed to be larger than the
  # potential of any configuration of pixel values.
  V_max <- dim(src)[1] * dim(src)[2]  * ((4)^2 / (2*covar) + 4 * weight_diff * max_diff)
  
  for (i in  1 : iterations){
    
    # Switch source and destination buffers.
    if (s == 1){
      s <- 2
      d <- 1
    }
    else{
      s <- 1
      d <- 2
    }
    
    # Vary each pixel individually to find the
    # values that minimise the local potentials.
    for (r in 1 : dim(src)[1]){
      for (c in 1 : dim(src)[2]){
        
        V_local <- V_max
        min_val <- -1
        for (val in vals){
          
          # The component of the potential due to the known data.
          V_data <- (val - src[r,c])^2 / (2 * covar)
          
          # The component of the potential due to the
          # difference between neighbouring pixel values.
          V_diff <- 0;
          if (r > 1){
            V_diff <- V_diff + min( (val - buffer[r-1,c,s])^2, max_diff )
          }
          if (r < dim(src)[1]){
            V_diff <- V_diff + min( (val - buffer[r+1,c,s])^2, max_diff )
          }
          if (c > 1){
            V_diff <- V_diff + min( (val - buffer[r,c-1,s])^2, max_diff )
          }
          if (c < dim(src)[2]){
            V_diff <- V_diff + min( (val - buffer[r,c+1,s])^2, max_diff )
          }
          
          V_current <- V_data + weight_diff * V_diff
          
          if (V_current < V_local){
            min_val <- val
            V_local <- V_current
          }
          
        }
        
        buffer[r,c,d] <- min_val
      }
    }
    
    #i
  }
  
  buffer[,,d]
}
