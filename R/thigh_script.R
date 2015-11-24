
# Utility to calculate gridsize in mm which will give a certain height, h, in pixels
# for a given rectangle, r.
calcGridSize <- function(r, h){
  hv <- r[4,] - r[1,]
  lhv <- sqrt(sum(hv^2))
  mean(c(lhv/(h-.5), lhv/(h-1.5)))
}

generatePixImages <- function(datadir="data"){
  source(paste0(datadir, "../R/inkscape_utils.R"))
  phantom <- importIPhantom(paste0(datadir,"/thigh_inkscape.csv"))
  ap <- read.csv("data/acoustic_properties_thigh.csv", as.is=TRUE)
  scan_1 <- computeRectangle(phantom, c(100,100), pi/4, 127*.23)
  scan_2 <- computeRectangle(phantom, c(100,100), pi/4 + pi/8, 127*.23)
  scan_3 <- computeRectangle(phantom, c(90, 80), 0, 127*.23)
  scan_4 <- computeRectangle(phantom, c(180, 125), pi/8, 127*.23)
  scan_5 <- computeRectangle(phantom, c(125,140), pi/2+pi/4+pi/16, 127*.23)
  scan_6 <- computeRectangle(phantom, c(130,125), pi/2+pi/4-pi/8, 127*.23)
  scan_7 <- computeRectangle(phantom, c(130,125), pi/2-pi/16, 127*.23)
  scan_8 <- computeRectangle(phantom, c(180,125), pi/4+pi/16, 127*.23)
  scans <- paste0("scan_",1:8)
  for(n in 1:8){
    scan <- get(scans[n])
    for(res in 2^(4:7)){
      gridsize <- calcGridSize(scan, res)
      img <- pixelate(phantom, scan, gridsize)
      colnames(img) <- as.character(signif(gridsize*(1:ncol(img)-1),4))
      outf <- paste0(datadir,"/",scans[n],"_",res,".csv")
      write.csv(as.data.frame(img),outf, row.names=FALSE)
      writeLines(outf)
    }
  }
}

readScan <- function(scan, resolution, datadir="data"){
  temp <- as.matrix(read.csv(paste0(datadir,"/scan_",scan,"_",resolution,".csv"), 
                             as.is=TRUE, comment.char = "#"))
  colnames(temp) <- gsub("X","",colnames(temp))
  temp
}

readWideScan <- function(scan, resolution, datadir="data"){
  temp <- as.matrix(read.csv(paste0(datadir,"/wide_scan_",scan,"_",resolution,".csv"), 
                             as.is=TRUE, comment.char = "#"))
  colnames(temp) <- gsub("X","",colnames(temp))
  temp
}

tissue2Slowness <- function(img, datadir="data"){
  ap <- read.csv(paste0(datadir,"/acoustic_properties_thigh.csv"), 
                        as.is=TRUE, comment.char = "#")
  for(id in unique(as.vector(img))){
    idx <- img==id 
    img[idx] <- ap[ap$ID==id, "speed"]
  }
  1/img
}