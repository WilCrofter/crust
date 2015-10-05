edges <- list(m1=matrix(1,4,4), m2=matrix(c(1,1,-1,-1), 4, 4), m3=t(matrix(c(1,1,-1,-1),4,4)))
M <- sapply(edges, as.vector)
ihat <- matrix(0,365,61)
for(i in seq(1, 361, by=4)){
  for(j in seq(1, 57, by=4)){
    e <- img[i:(i+3), j:(j+3)]
    ihat[i:(i+3), j:(j+3)] <- matrix(M %*% solve(t(M) %*% M, t(M) %*% as.vector(e)), 4, 4)
  }
}
# image(x=(1:364)/2, y=(1:60)/2, z=ihat[1:364,1:60], asp =1)

changeOfVariable <- function(edges, height, width){
  M <- sapply(edges, as.vector)
  p <- solve(t(M) %*% M) %*% t(M)
  by <- dim(edges[[1]])[1]
  P <- matrix(0, length(edges)*ceiling(height/by)*ceiling(width/by), height*width)
  rowP <- 0
  tbl <- cbind(rep(1:by,by),rep(1:by,each=by))
  for(j in seq(1, height, by=by)){
    for(i in seq(1, width, by=by)){
      if((tbl[dim(p)[2],1] + i - 1 > width) | (tbl[dim(p)[2],2]+j-1 > height))continue
      for(k in 1:(dim(p)[1])){
        rowP <- rowP + 1
        for(m in 1:(dim(p)[2])){
          iCol <- tbl[m,1] + i - 1
          jCol <- tbl[m,2] + j - 1
          colP <- (jCol-1)*width + iCol
          P[rowP, colP] <- p[k,m]
        }
      }
    }
  }
  P
}

recoverImg <- function(coefs, edges, height, width){
  M <- sapply(edges, as.vector)
  by <- dim(edges[[1]])[1]
  img <- matrix(0, width, height)
  patches <- M %*% matrix(coefs, ncol(M), length(coefs)/ncol(M))
  patchIdx <- 0
  for(j in seq(1, height, by=by)){
    for(i in seq(1,width, by=by)){
      if((i+by-1) > width | (j+by-1) > height ){
        print("hi")
      }
      patchIdx <- patchIdx + 1
      img[i:(i+by-1), j:(j+by-1)] <- matrix(patches[,patchIdx], by, by)
    }
  }
  img
}

# Deprecated stuff related to projecting an image on the orthogonal complement
# of the stripe space. Actually, it's easy. If img is the image **matrix**, i.e.,
# its rows are x, and its columns y in a picture of the image, then
#  1) subtract the mean of each row from the row
#  2) add the projection of the image on the normalized constant image

formPerp <- function(m){
  sapply(1:(m-1), function(k)c(rep(0,k-1),1,-1,rep(0,m-(k+1))))
}

formPerpProj <- function(m){
  M <- formPerp(m)
  M %*% solve(t(M) %*% M, t(M))
}

img1 <- img
for(i in 1:nrow(img))img1[i,] <- img[i,]-mean(img[i,])
