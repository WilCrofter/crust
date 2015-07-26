#' Utilities related to the inner product <u,v> = t(Su)Sv on the orthogonal
#' complement of the kernel of S.
#' 
#' Sets of vectors which are orthonormal with respect to this inner product may be
#' useful for analysis of system matrices and the potential for image recovery in
#' co-robotic ultrasound tomography. If U is such an orthonormal set, then the
#' least squares solution to Sb ~ tau in the span of U is just the sum of the 
#' projections of t(S)tau on the members of U, where the projection on u is just 
#' <t(S)tau,u>*u. Since it is a subspace of all possible images, the minimum over
#' the span of U is not likely to be a global minimum. However, it will be unique. 
#' In fact, it is easily shown that the span of U is in the orthogonal complement
#' of the kernel of S, hence the minimum on the span of U is the projection of all
#' global minima. 
#' 
#' The kernel of S includes, though is not necessarily limited to, the so-called stripe space.
#' The stripe space consists of images (i.e., their associated vectors) such that
#' each row of pixels is identical, and the sum of values in a row is zero. The
#' orthogonal complement of the stripe space consists of images such that each
#' column of pixels sums to the same value as every other.

#' inner product with respect to S
inner <- function(u, S, v){
  if(is(u, "test image"))u <- as.vector(u)
  if(is(v, "test image"))v <- as.vector(v)
  if(is.vector(u))u <- matrix(u, length(u), 1)
  if(is.vector(v))v <- matrix(v, length(v), 1)
  if(!is.matrix(u) & !is(u, "Matrix"))stop("u is not a test image, matrix, Matrix, or vector")
  if(!is.matrix(v) & !is(v, "Matrix"))stop("v is not a test image, matrix, Matrix, or vector")
  t(S %*% u) %*% (S %*% v)
}

#' Euclidean norm with respect to S
eunorm <- function(u, S)sqrt(sum((S %*% as.vector(u))^2))

#' Given a matrix, U, whose columns are orthonormal in the inner product
#' associated with S, return a vector v, orthonormal to the columns of U,
#' such that the span of v and the columns of U includes b.
extendONS <- function(b, S, U=NULL){
  if(is.null(U))return(as.vector(b)/eunorm(b, S))
  if(is.vector(U))dim(U)<- c(length(U), 1)
  # find the residual of the projection of b on the column space of U.
  ans <- as.vector(b) - U %*% inner(U, S, b)
  # return its normalized value
  ans/eunorm(ans, S)
}