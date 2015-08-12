# Utilities for phantoms constructed using Inkscape

if(!require(sp))stop(paste("This script requires the sp package,",
                           "which does not appear to be installed.",
                           "Please install it and source this script again."))

### Converting Inkscape paths to scaled, R polygons.

ink2R <- function(xml_string){
  code <- substr(xml_string,1,1)
  if(code != "m" & code != "M")stop(paste("Unrecognized code:", code))
  xml_string <- substr(xml_string, 3, nchar(xml_string)-2)
  xml_string <- paste0("c(", gsub(" ", ",", xml_string), ")")
  pts <- t(matrix( eval(parse(text=xml_string)), nrow=2))
  if(code == "m"){
    pts[,1] <- cumsum(pts[,1])
    pts[,2] <- cumsum(pts[,2])
  }
  pts <- rbind(pts, pts[1,])
  # Inkscape's y coordinate is in the direction opposite to R's.
  pts[,2] <- -pts[,2]
  pts
}

deriveXform <- function(pgon, max_coordinate=200){
  r1 <- range(pgon[,1])
  r2 <- range(pgon[,2])
  # translation of minimum coordinates to 0
  shift <- -c(r1[1], r2[1])
  # scale maximum coordinate
  scale <- max_coordinate/(max(diff(r1), diff(r2)))
  list(shift=shift, scale=scale)
}

applyXform <- function(pgon, xform){
  pgon[,1] <- pgon[,1] + xform$shift[1]
  pgon[,2] <- pgon[,2] + xform$shift[2]
  pgon*xform$scale
}
