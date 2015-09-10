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

formH <- function(nr, nc){
  d <- nr*nc
  H <- Matrix(0, d, d)
  for(i in 1:nr){
    for(j in 1:(nc-1)){
      r <- ij2col(i,j,nr,nc)
      H[r,r] <- 1
      s <- ij2col(i,j+1,nr,nc)
      H[r,s] <- -1
    }
  }
  t(H) %*% H
}

formV <- function(nr, nc){
  d <- nr*nc
  V <- Matrix(0, d, d)
  for(i in 1:(nr-1)){
    for(j in 1:nc){
      r <- ij2col(i,j,nr,nc)
      V[r,r] <- 1
      s <- ij2col(i+1,j,nr,nc)
      V[r,s] <- -1
    }
  }
  t(V) %*% V
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
