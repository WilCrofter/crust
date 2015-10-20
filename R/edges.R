edges <- list(m1=matrix(1,4,4), m2=matrix(c(1,1,-1,-1), 4, 4), m3=t(matrix(c(1,1,-1,-1),4,4)))
# M <- sapply(edges, as.vector)
# ihat <- matrix(0,365,61)
# for(i in seq(1, 361, by=4)){
#   for(j in seq(1, 57, by=4)){
#     e <- img[i:(i+3), j:(j+3)]
#     ihat[i:(i+3), j:(j+3)] <- matrix(M %*% solve(t(M) %*% M, t(M) %*% as.vector(e)), 4, 4)
#   }
# }


tileProjection <- function(edges, height, width){
    by <- dim(edges[[1]])[1]
    M <- numeric()
    for(j in seq(1, height, by=by)){
      for(i in seq(1, width, by=by)){
        if((by + i - 1 > width) | (by + j-1 > height))next
        for(edge in edges){
          temp <- matrix(0, width, height)
          temp[i:(by+i-1), j:(by+j-1)] <- edge
          M <- cbind(M, as.vector(temp))
        }
      }
    }
    M %*% solve(t(M) %*% M) %*% t(M)
}

stripeProjection <- function(height, width){
  M <- numeric()
  for(i in 2:width){
    temp <- matrix(0, width, height)
    temp[1,] <- 1
    temp[i,] <- -1
    M <- cbind(M,as.vector(temp))
  }
  M %*% solve(t(M) %*% M) %*% t(M)
}

tryall <- function(tau, simg, S, ITP, SP, tWt, sWt){
  STS <- t(S) %*% S
  bhat <- numeric(dim(S)[2])
  bhat <- solve(STS + tWt*ITP + sWt*SP,t(S) %*% tau)
  cor(as.vector(simg),as.vector(bhat))
}

