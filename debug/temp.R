source("R/smooth.R")
source("C/ctest.R")
source("R/thigh_script.R")
source("R/utilities.R")
res <- 16
img <- tissue2Slowness(readWideScan(3,res))
gridsize <- .46*128/res
height <- res
width <- nrow(img)
x <- seq(1,width)*gridsize
y <- seq(1,height)*gridsize
image(x,y,img,asp=1)
S <- genS(height, width, gridsize)
Sc <- genStest(height, width, gridsize)
Sp <- S[rangeRows(0,height-1,height),]
P <- t(Sp) %*% solve(Sp %*% t(Sp), Sp)
proj <- matrix(P %*% as.vector(img), width, height)
image(x,y,proj,asp=1)
temp <- t(apply(proj,1,diff))
image(x,y[-1],temp,asp=1)
