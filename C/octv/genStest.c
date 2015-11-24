#include "mex.h"
#include "genS.c"

void
mexFunction (int nlhs, mxArray* plhs[],
             int nrhs, const mxArray* prhs[])
{
  
  int *pheight, *pwidth, height, width;
  double *pgridsize, gridsize, *S;

  if (nrhs != 3 )
    mexErrMsgTxt ("3 args please");

  height = mxGetScalar(prhs[0]);
  width =  mxGetScalar(prhs[1]);
  pgridsize = mxGetPr(prhs[2]);
  
  pheight = &height;
  pwidth = &width;
  gridsize = *pgridsize;
  
  mexPrintf ("I have %d inputs and %d outputs\n", nrhs, nlhs);
 
  mexPrintf("value of height is %x value of width is %x gridsize is %f\n",
             height,width,gridsize);
  mexPrintf ("height adddress is %x holding %d width is %x grid is %x %f\n",
              pheight,*pheight,*pwidth,pgridsize,gridsize);
  S = mxCreateDoubleMatrix (height^2,height*width,mxREAL);
  CgenS(S,pheight,pwidth,pgridsize);
  plhs[0]=S;
}
