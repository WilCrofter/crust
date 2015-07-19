# Trickery to handle sourcing from an Rmd
if(isTRUE(file.info("R")$isdir)){
  # working directory contains R
  source("R/utilities.R")
} else if(isTRUE(file.info("../R")$isdir)){
  # parent of working directory contains R
  source("../R/utilities.R")
}

# Specify a grid by number of rows, number of columns, and spacing
# The horizontal direction is x, 
newGrid <- function(nrows, ncols, spacing){
  
}


# Modifies createSij to work independently of handy scanner.
formSij <- function(setup){
  npix <- setup$npix
  numrows <- setup$n*(setup$n-1)
  Sij <- matrix(0,nrow=numrows,ncol=npix^2)
  for (i in 1:(setup$n-1)) 
    for (j in 1:setup$n){
      row <- pixLengths(i,j,setup)
      for (k in 1:length(row$segment_length)) {
        idx <- (row$y_index[k]-1) * npix + row$x_index[k]
        # Sij[(i-1)*setup$n + j,idx] <- row$segment_length[k]
        Sij[(j-1)*(setup$n-1) + i, idx] <- row$segment_length[k]
      }
    }
  #now do some checks on Sij
  rsums <- rowSums(Sij)
  
  for (i in 1:numrows) {
    nz <- sum(Sij[i,]>0)
    if (nz<npix)stop("too few nonzero elements in row ",i)
    cl <- floor((i-1)/(setup$n-1))+1
    rw <- (i-1)%%(setup$n-1) + 1
    dist <- sqrt((setup$align$transmitters[1,rw]-setup$align$receivers[1,cl])^2 +
                   (setup$align$transmitters[2,rw]-setup$align$receivers[2,cl])^2 +
                   (setup$align$transmitters[3,rw]-setup$align$receivers[3,cl])^2)
    if (!isTRUE(all.equal(dist,rsums[i])))stop("pathlengths not equal in row ",i," ",dist," ",rsums[i])
  }
  Sij
}
