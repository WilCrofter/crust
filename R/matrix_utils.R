#find identical rows in nr x nc matrix mat
findIdentiRows <- function(mat,nr,nc){
  for (i in 1:(nr-1))
    for (j in (i+1):nr){
  #    ctr <- sum(mat[i,]==mat[j,])
   #   if (ctr==nc) print(paste("Rows ",i,"and",j,"are identical",sep=" "))
     if (isTRUE(all.equal(mat[i,],mat[j,])))  
       print(paste("Rows ",i,"and",j,"are identical",sep=" "))
    }
  
  
}

#Find how many 0 eigenvalues for fixed grid size and varying number of transducers
howMany0 <- function(sm,lg,gridsz){
  if (sm>lg)stop("smallest must be smaller or less than largest")
  for (j in sm:lg){
    out <- handySetup(j,15,gridsz)
    S <- createSij(out)
    STS <- t(S) %*% S
    eig <- eigen(STS,symm=TRUE)
    print(paste(j,gridsz,sum(eig$values<1e-7),sep="  "))
  }
}

fixS <- function(out){
  S <- createSij(out)
  cl <- dim(S)[2]
  qr <- qr(S)
  if (qr$rank<cl){
    rmcl <- tail(qr$pivot,cl-qr$rank)
    cutS <- S[,-rmcl]
  }
  cutS
}

fixImage <- function(noisy){
  nrw <- nrow(noisy)
  ncl <- ncol(noisy)
  lessnoisy <- matrix(0,nrw,ncl)
  tot <- 0
  for (j in 1:(nrw-1)){
    ad <- .5*(max(noisy[j,]-noisy[j+1,])+min(noisy[j,]-noisy[j+1,]))
    tot = tot + ad
    if (ad>0) lessnoisy[j,] <- noisy[j,]-ad
    else lessnoisy[j,] <- noisy[j,]+ad
    print(ad)
  }
  if (tot>0) lessnoisy[nrw,] <- noisy[nrw,]-ad
  else lessnoisy[nrw,] <- noisy[nrw,]+ad
  lessnoisy
}

getY <- function(setup){
  yt <- recalcForScan(setup$phantom,setup$xmitr,setup$rcvr)$transmitters[2,]
  yr <- recalcForScan(setup$phantom,setup$xmitr,setup$rcvr)$receivers[2,]
  list(yt=yt,yr=yr)
}

xmitrcvrDist <- function(setup){
  dist <- matrix(0,setup$n-1,setup$n)
  for (i in 1:setup$n-1)
    for (j in 1:setup$n)
      dist[i,j] <- sqrt( (setup$align$transmitters[1,i]-setup$align$receivers[1,j])^2 +
                           (setup$align$transmitters[2,i]-setup$align$receivers[2,j])^2 +
                           (setup$align$transmitters[3,i]-setup$align$receivers[3,j])^2 
        )
  dist
}

getXR_fromRow <- function(rownumber,setup){
  X <- (rownumber-1)%%(setup$n-1) + 1
  R <- floor((rownumber-1)/(setup$n-1))+1
  ans <- cbind(rownumber,X,R) 
  ans
}


#given pixel (x,y) find all the xmitr/rcvr pairs that go through it
getXR_fromPixel <- function(x,y,setup,myS=S){
  XR <- matrix(0,0,3)
  if (x<1 || x>setup$npix) stop("x is out of bounds ",x)
  if (y<1 || y>setup$npix) stop("y is out of bounds ",y)
  if (is.null(myS)) S <- createSij(setup)
  else S <- myS
  #find which column of S corresponds to pixel (x,y)
  idx <- (y-1) * setup$npix + x
  for (i in 1:(setup$n*(setup$n-1)))
    if (S[i,idx]>0){
      cl <- floor((i-1)/(setup$n-1))+1
      rw <- (i-1)%%(setup$n-1) + 1
      XR <- rbind(XR,c(i,rw,cl))
    }
  XR
  }

tryAll <- function(tof){
  len <- length(tof)
  maxtry <- 100
  put <- numeric()
  stor <- rep(1000,len)
  old <- numeric()
  new <- numeric()
  aidx <- matrix(0,len,3)
  tof <- tof*maxtry
  c_silicone <- 1000.0 # m/sec
  c_alcohol <- 1180.0 # m/sec (ethyl alcohol)
  c_water <- 1480.0 # m/sec
  for (i in 0:maxtry){
    for (j in 0:(maxtry-i)){
      k <- maxtry-i-j
      put <- rep(c_silicone*i + c_alcohol*j + c_water*k,len)
      old <- abs(tof-stor)
      new <- abs(tof-put )
      for (idx in 1:len) 
        if (new[idx]<old[idx]) {
        stor[idx] <- put[idx]
        aidx[idx,1] <- i
        aidx[idx,2] <- j
        aidx[idx,3] <- k
    }
  }
  }
  list(best=stor,indx=aidx)
}