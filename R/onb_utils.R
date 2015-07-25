#' Utilities related to the inner product <u, v> = t(Su)Sv on the orthogonal
#' complement of the kernel of S, which we abbreviate here as S-perp.
#' 
#' Sets of vectors which are orthonormal with respect to this inner product may be
#' useful for analysis of system matrices and the potential for image recovery in
#' co-robotic ultrasound tomography. If U is such an orthonormal set, then the
#' least squares solution to Sb ~ tau in U is just the sum of the projections of
#' t(S)tau on the members of U, where the projection on u is just <t(S)tau, u>u.
#' 
#' Ker(S) includes, though is not necessarily limited to, the so-called stripe space.
#' The stripe space consists of images (i.e., their associated vectors) such that
#' each row of pixels is identical, and the sum of values in a row is zero. The
#' orthogonal complement of the stripe space consists of images such that each
#' column of pixels sums to the same value as every other. We'll abbreviate this
#' space as stripe-perp.
#' 
#' It appears true, though has not been proven, that if the transducer-to-pixel
#' ratio is large enough, ker(S) consists of the stripe space alone. Thus, it
#' seems to make sense to form orthonormal sequences, u1, u2, u3, ..., in
#' stripe-perp such that shorter subsequences divide the image into cruder
#' blocks of pixels than do longer subsequences. These utilities are meant to
#' facilitate that.

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
extend <- function(b, S, U=NULL){
  if(is.null(U))return(as.vector(b)/eunorm(b))
  # find the residual of the projection of b on the column space of U.
  ans <- b - U %*% inner(U, S, b)
  # return its normalized value
  ans/eunorm(ans)
}