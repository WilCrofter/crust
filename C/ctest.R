# MAKE SURE HEIGHT and WEiGHT are properly defined in genS.h file
# BEFORE you compile with R CMD SHFILE command

# CgenS is the name of the C function contained
# in the C file genS.c. CgenS takes 2 arguments
# which follow it in the .C line.
# The first is the address of the S matrix which
# has height^2 rows and height*width columns. It is 
# of type double.
# The second argument, also of type double, is the gridsize.
# It looks like passing by value works, since the C function expects
# a value not an address.
# The C program is compiled from a terminal window with the 
# command R CMD SHLIB genS.c. This invokes the C compiler and 
# creates the shared object file genS.so.
# This last must be dynamically loaded from R using
# dyn.load("genS.so"). Make sure you specify the correct path to 
# this .so file, e.g. dyn.load("C/genS.so")

#dyn.load("C/genS.old.so")
dyn.load("C/genS.so")
genStest <- function(height,width,gridsize){
  # SfromR <- genS(height,width,gridsize)
  # The following shouldn't have worked since R calls C by passing addresses, not values
  # and oldCgenS expected a value for gridsize, not a pointer!!
  # SfromC <- .C("oldCgenS", S=double(height^3*width),gridsize)
  
  SfromC <- .C("CgenS",S=double(height^3*width),as.integer(height),as.integer(width),as.double(gridsize))
  matrix((SfromC["S"])[[1]],height^2,height*width,byrow=TRUE)
}