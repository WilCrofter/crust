#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG1
#define HEIGHT 8  //y dim of image
#define WIDTH  12  //x dim of image

#define HEIGHT1 (HEIGHT+1)
#define WIDTH1  (WIDTH+1)
#define TOTAL   (HEIGHT1+WIDTH1)
#define SQR(A)  ((A)*(A))
struct point{double x; double y; double z;};
struct ipt{int x; int y; int z;};
struct segment{int prow, pcol; double length;};
