library(Matrix)

# Determine S matrix column associated with pixel i,j
ij2col <- function(i, j, nr, nc){
  (j-1)*nr + i
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

# NOTE: n is the number of rows of pixels in the IMAGE, m the number of columns.
# BUT, n is the number of COLUMNS in the corresponding matrix.
genS <- function(n, m, gridsize){
  S <- Matrix(0, n^2, n*m)
  for(i in 1:n){
    for(j in 1:n){
      u <- c(-1e-6, (i-0.5)*gridsize, 0)
      v <- c(m*gridsize+1e-6, (j-0.5)*gridsize, 0)
      segs <- segmentLengths(m, n, u, v, gridsize)
      ridx <- ij2row(i,j,n)
      for(k in 1:nrow(segs)){
        S[ridx, ij2col(segs[k,2], segs[k,1], n, m)] <- segs[k,3]
      }
    }
  }
  S
}

pcsd <- function(img, beta){
  sd( (as.vector(beta)-as.vector(img))/as.vector(img) )
}
