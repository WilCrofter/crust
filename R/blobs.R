# source('~/dev/r/crust/C/ctest.R')
# source('~/dev/r/crust/R/smooth.R')
# S <- genStest(32, 200, 4*.46)
# Sp <- S[rangeRows(0, 32, 32), ]

blobs <- function(Splus, x, y, size=8, gridsize=.92){
  img <- matrix(0, 200, 32)
  temp <- img[x:(x+size-1), y:(y+size-1)] <- 1
  par(mfrow=c(2,1))
  image(x=(1:200)*.46, y=(1:32)*.46, z=img, asp = 1, xlab="x (mm)", ylab="y (mm)", main="image")
  z=matrix(P%*%as.vector(img),200,32)
  image(x=(1:200)*.46, y=(1:32)*.46, z=z, asp = 1, xlab="x (mm)", ylab="y (mm)", main="projection")
  par(mfrow=c(1,1))
  print(cor(as.vector(img), as.vector(z)))
  mean(z[x:(x+size-1), y:(y+size-1)])
}