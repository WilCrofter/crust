#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG1

#define SQR(A)  ((A)*(A))

struct point{double x; double y; double z;};
struct segment{int prow, pcol; double length;};
