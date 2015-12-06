#include "mex.h"
//#include "genS.c"
void CgenS(double* S, int* pheight, int* pwidth, double* gridsize);
void mexFunction (int nlhs, mxArray* plhs[],
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
  plhs[0] = mxCreateDoubleMatrix (height^3*width,1,mxREAL);
  S=mxGetPr(plhs[0]);
  
  mexPrintf("S address is %x plhs[0] is %x\n",S,plhs[0]);
  
  CgenS(S,pheight,pwidth,pgridsize);
 
}
