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

importIPhantom <- function(file, max_coordinate=200){
  paths <- read.table(file, header=TRUE, as.is = TRUE, sep=",")
  k <- which(paths$layer == 1)[1]
  xform <- deriveXform(ink2R(paths[k,"path"]), max_coordinate=max_coordinate)
  ans <- list()
  for(n in 1:nrow(paths)){
    temp <- list()
    for(name in names(paths)[-ncol(paths)]){
      temp[[name]] <- paths[n, name]
    }
    shape <- applyXform(ink2R(paths[n,"path"]), xform)
    temp[["shape"]] <- shape
    ans[[paths[n,"description"]]] <- temp
  }
  ans
}

plotIPhantom <- function(phantom,
                         ccode=c("gray", "lightgray", "tan", "blue",
                                 "red", "yellow", "darkgray","brown"),
                         col="lightblue", asp=1, xlab="x (mm)", ylab="y (mm)",
                         main="IPhantom", ...){
  layers <- sapply(phantom, function(x)x$layer)
  k <- which(layers==1)
  rx <- range(phantom[[k]]$shape[,1])
  ry <- range(phantom[[k]]$shape[,2])
  plot(rx, ry, type='n', asp=asp, xlab=xlab, ylab=ylab, main=main,...)
  rect(rx[1], ry[1], rx[2], ry[2], col=col)
  for(n in order(layers)){
    x <- phantom[[n]]
    polygon(x$shape[,1], x$shape[,2], col=ccode[x$layer])
  }
}

computeRectangle <- function(phantom, center, angle_in_radians, height, ...){
  a <- angle_in_radians # for brevity
  # rotation by a (as applied from the left)
  rot <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  # find the phantom's boundary 
  k <- which(1 == sapply(phantom, function(x)x$layer))
  bdry <- phantom[[k]]$shape
  # calculate a precision with which to determine rectangle boundaries
  eps <- abs(diff(range(bdry)))/100
  delta <- rot %*% c(eps, 0)
  # find first points outside phantom, working from center,
  # returning the associated multiplier of the increment.
  fct <- function(start, increment){
    pt <- start
    k <- 0
    while(point.in.polygon(pt[1], pt[2], bdry[,1], bdry[,2])){
      k <- k+1
      pt <- pt + increment
    }
    k
  }
  # half the height, properly rotated
  hh <- rot %*% (height*c(0,.5))
  # upper left increment
  kul <- fct(center+hh, -delta)
  # upper right
  kur <- fct(center+hh, delta)
  # lower left
  kll <- fct(center-hh, -delta)
  # lower right
  klr <- fct(center-hh, delta)
  # these increment multiples define a trapezoid, not a rectangle,
  # so define a miminal rectangle which includes them
  kl <- max(kul, kll)
  ul <- center+hh-kl*delta
  ll <- center-hh-kl*delta
  kr <- max(kur, klr)
  ur <- center+hh+kr*delta
  lr <- center-hh+kr*delta
  # return corners
  matrix(c(ul,ur,lr,ll), 4, 2, byrow = TRUE)
}

overlayRectangle <- function(r,...){
  polygon(rbind(r,r[1,]),...)
}

pixelate <- function(phantom, rectangle, gridsize){
  # assuming rectangle corners are arranged in clockwise order,
  # compute oriented base and height
  bv <- r[2,]-r[1,]
  hv <- r[4,]-r[1,]
  # rectangle dimensions in integer multiples of gridsize
  b <- ceiling(.5+sqrt(sum(bv^2))/gridsize)
  h <- ceiling(.5+sqrt(sum(hv^2))/gridsize)
  # gridsize increments for base and height
  bi <- gridsize*bv/(sqrt(sum(bv^2)))
  hi <- gridsize*hv/(sqrt(sum(hv^2)))
  # order the phantom shapes by layer
  ord <- order(sapply(phantom, function(x)x$layer))
  # function to determine tissue type of a point
  tissue <- function(pt){
    ID <- 0
    for(n in ord){
      shape <- phantom[[n]]$shape
      if(point.in.polygon(pt[1], pt[2], shape[,1], shape[,2])){
        ID <- phantom[[n]]$ID
      }
    }
    ID
  }
  img <- matrix(0,b,h)
  # assuming r[1,] is a midpoint, compute remaining midpoints and
  # determine their tissue types
  for(i in 1:b){
    b0 <- r[1,] + (i-1)*bi
    for(j in 1:h){
      pt <- b0 + (j-1)*hi # midpoint for jth column in row
      img[i,j] <- tissue(pt)
    }
  }
  img <- img[,ncol(img):1] # was constructed upside down
  attr(img, "base") <- (b+1)*gridsize
  attr(img, "height") <- (h+1)*gridsize
  img
}

pixelatedImage <- function(img, col=NULL, xlab="x (mm)", ylab="y (mm)", ...){
  if(is.null(col))col <- gray.colors(length(unique(as.vector(img))))
  image(x=seq(0,attr(img,"base"), length.out=dim(img)[1]),
        y=seq(0,attr(img,"height"), length.out=dim(img)[2]), z=img, col=col, asp=1,
        xlab=xlab, ylab=ylab)
}