library(Matrix)

ij2col <- function(i, j, nr, nc){
  (j-1)*nr + i
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
