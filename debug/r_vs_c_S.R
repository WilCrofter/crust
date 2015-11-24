source("R/smooth.R")
source("C/ctest.R")
source("R/thigh_script.R")
source("R/utilities.R")
## Example in which R and C generate different system matrices. 
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
max(abs(S-Sc))/max(abs(S))

