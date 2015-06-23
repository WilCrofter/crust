#find identical rows in nr x nc matrix mat
findIdentiRows <- function(mat,nr,nc){
  for (i in 1:(nr-1))
    for (j in (i+1):nr){
      ctr <- sum(mat[i,]==mat[j,])
      if (ctr==nc) print(paste("Rows ",i,"and",j,"are identical",sep=" "))
      
    }
  
  
}
getY <- function(setup){
  yt <- recalcForScan(setup$phantom,setup$xmitr,setup$rcvr)$transmitters[2,]
  yr <- recalcForScan(setup$phantom,setup$xmitr,setup$rcvr)$receivers[2,]
  list(yt=yt,yr=yr)
}