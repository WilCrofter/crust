library(Matrix)

# Determine S matrix column associated with pixel i,j
#
# pixrow is the pixel ROW in the IMAGE as it appears in a call to image,
# numbered from bottom up. It corresponds to the COLUMN of the
# image regarded as a MATRIX. In other words, ipix varies in the
# vertical dimension of the image from bottom up.
#
# pixcol is the pixel COLUMN (matrix row), and varies from right
# to left in the image.
# 
# height is the HEIGHT of the IMAGE, the number of pixel rows.
# width is , the WIDTH of the image, the number of pixel columns.
ij2col <- function(pixrow, pixcol, height, width){
  # (j-1)*nr + i
  (pixrow-1)*width + pixcol
}

# Determine the S matrix row associated with transducer pair i, j
ij2row <- function(i, j, ntransducers){
  (j-1)*ntransducers + i
}

formH <- function(height, width){
  d <- height*width
  H <- Matrix(0, d, d)
  for(pixrow in 1:height){
    for(pixcol in 1:(width-1)){
      r <- ij2col(pixrow,pixcol,height,width)
      H[r,r] <- 1
      s <- ij2col(pixrow,pixcol+1,height,width)
      H[r,s] <- -1
    }
  }
  H
}

formV <- function(height, width){
  d <- height*width
  V <- Matrix(0, d, d)
  for(pixrow in 1:(height-1)){
    for(pixcol in 1:width){
      r <- ij2col(pixrow,pixcol,height,width)
      V[r,r] <- 1
      s <- ij2col(pixrow+1,pixcol,height,width)
      V[r,s] <- -1
    }
  }
  V
}

# NOTE: height is the number of rows of pixels in the IMAGE, width the number of columns.
genS <- function(height, width, gridsize){
  S <- Matrix(0, height^2, height*width)
  for(i in 1:height){
    for(j in 1:height){
      u <- c(-1e-6, (i-0.5)*gridsize, 0)
      v <- c(width*gridsize+1e-6, (j-0.5)*gridsize, 0)
      segs <- segmentLengths(width, height, u, v, gridsize)
      ridx <- ij2row(i,j,height)
      for(k in 1:nrow(segs)){
        # segs[k,2] is the IMAGE Y dimension (height) which corresponds to pixel ROW
        # segs[k,1] is the IMAGE X dimension (width) which corresponds to pixel COL
        S[ridx, ij2col(segs[k,2], segs[k,1], height, width)] <- segs[k,3]
      }
    }
  }
  S
}

pcmean <- function(img, beta){
  mean(abs( (as.vector(beta)-as.vector(img))/as.vector(img) ))
}

# Placing the patch at pixrow, slide it horizontally, apply S at each offset, and
# save the correlation of the result with tau.
slider <- function(img, patch, pixrow, tau){
  width <- dim(img)[1]
  height <- dim(img)[2]
  pwidth <- dim(patch)[1]
  pheight <- dim(patch)[2]
  ans <- numeric()
  for(pixcol in 1:(width-pwidth+1)){
    pimg <- matrix(1/1450, width, height)
    pimg[pixcol:(pixcol+pwidth-1), pixrow:(pixrow+pheight-1)] <- patch
    ans <- c(ans, cor(as.vector(S %*% as.vector(pimg)), as.vector(tau)))
  }
  ans
}
