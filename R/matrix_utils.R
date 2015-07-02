#find identical rows in nr x nc matrix mat
findIdentiRows <- function(mat,nr,nc){
  for (i in 1:(nr-1))
    for (j in (i+1):nr){
      ctr <- sum(mat[i,]==mat[j,])
      if (ctr==nc) print(paste("Rows ",i,"and",j,"are identical",sep=" "))
      
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
#given pixel (x,y) find all the xmitr/rcvr pairs that go through it
getXR <- function(x,y,setup,myS=S){
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