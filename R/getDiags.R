wdiag <- function(mat,w){
  if (nrow(mat)!=ncol(mat)) stop("matrix is not square")
  n <- nrow(mat)
  if (w>n) stop("specified diagonal is too big")
  
  ans <- numeric()
  
  for (i in 1:(n-w+1)){
    ans <- c(ans,mat[i,i+w-1])
  }
  if (w>1){
    for (i in (n-w+2):n){
      ans <- c(ans,mat[i,1+i-(n-w+2)])
    }
  }
  ans
  
  
}