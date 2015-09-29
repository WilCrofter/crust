# random horizontal difference
dij <- function(n, img, delta=1){
  i <- sample(dim(img)[1]-delta, n, replace=TRUE)
  j <- sample(dim(img)[2], n, replace=TRUE)
  sapply(1:n, function(n){(img[i[n],j[n]] - img[i[n]+delta, j[n]])^2})
}

sumdij <- function(img, n=100, delta=1){
  sum(dij(n, img, delta))
}

wH <- 1e-12 # 1e-12 and wH/2 were best
wV <- wH/3
bhat2 <- solve(STS + wH*HTH+wV*VTV, t(S) %*% tau)
ihat2 <- as.testImage(matrix(bhat2,  48, 16), 1/2)
plotTestImage(ihat2, main=paste("wH = ", wH, "wV = ", wV))
cor(as.vector(bhat2), as.vector(img))
