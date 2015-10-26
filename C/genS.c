#include "genS.h"

void CgenS(S,pheight,pwidth,pgridsize)
double *pgridsize,S[];
int *pheight,*pwidth;
{
  int i,j,k,ridx,cidx,retval;
  int total,npix;
  int height,width;
  int segLengths();
  double gridsize;
  struct point u,v;
  struct segment *segs;

#ifdef DEBUG
  printf("Hello height is %d width is %d grid is %f\n",height,width,gridsize);
#endif

  height = *pheight;
  width = *pwidth;
  gridsize = *pgridsize;

  npix = height*width;
  total = (height+width+2);
  segs = calloc(total,sizeof(segs[0]));

  u.x = -1e-6;
  v.x = (width*gridsize)+1e-6;
  u.z = v.z = 0;
  for (i=0;i<(height*height);i++)
    for (j=0;j<(npix);j++)      S[i*npix+j]=0.0;

  for (i=1;i<=height;i++){ //for each transmitter
    u.y = (i-0.5)*gridsize;
    for (j=1;j<=height;j++){ //for each receiver
      v.y = (j-0.5)*gridsize;
      retval = segLengths(u,v,height,width,gridsize,segs); 
#ifdef DEBUG
      printf("retval is %d\n",retval);
#endif
      ridx=(j-1)*height + i - 1;
      for (k=0;k<retval;k++){
	cidx=(segs[k].pcol-1)*width + segs[k].prow -1;
	S[(ridx*npix)+cidx]=segs[k].length;
      }//

    }//j
  }//i
  free(segs);
}

int segLengths(u,v,height,width,gridsize,segs)
     struct point u,v;
     double gridsize;
     struct segment segs[];
     int    height,width;
{
  double        i1[2],i2[2],*kx,*ky;
  double        xdiff,ydiff,*xlambda,*ylambda;
  double        *xylambda,*lengths;
  int           i,xylen,xlen,ylen,retval;
  int           total;
  int           wCrossings(),removeInvalid(),merge();
  struct point uprime,vprime,*crossings,*interiors;
  int          *validCross;
  void         computeCrossings(),swap();

#ifdef DEBUG
  printf("in segLengths\nu is %f %f %f\nv is %f %f %f\n",
	 u.x,u.y,u.z,v.x,v.y,v.z);
#endif

  total = (height+width+2);

  kx = calloc(width+1,sizeof(kx[0]));
  ky = calloc(height+1,sizeof(ky[0]));
  xlambda = calloc(width+1,sizeof(xlambda[0]));
  ylambda = calloc(width+1,sizeof(ylambda[0]));
  xylambda = calloc(total,sizeof(xylambda[0]));
  lengths = calloc(total,sizeof(lengths[0]));
  crossings = calloc(total,sizeof(crossings[0]));
  interiors = calloc(total,sizeof(interiors[0]));
  validCross = calloc(total,sizeof(validCross[0]));

  xdiff = v.x-u.x;
  ydiff = v.y-u.y;

  //make sure points are distinct
  if ((u.x==v.x)&&(u.y==v.y)) {return(1);}

  //make sure line connecting points crosses grid interior
  i1[0] = (-u.x)/xdiff;
  i1[1] = (width*gridsize-u.x)/xdiff;
  if (i1[1]<i1[0]) swap(&i1[0],&i1[1]);

  i2[0] = (-u.y)/ydiff;
  i2[1] = (height*gridsize-u.y)/ydiff;
  if (i2[1]<i2[0]) swap(&i2[0],&i2[1]);

#ifdef DEBUG
  printf("i1 is %f %f   i2 is %f %f\n",i1[0],i1[1],i2[0],i2[1]);
#endif

  if ((i1[0]>=i2[1]) || (i2[0]>=i1[1])) {return(2);}

  for (i=0;i<total;i++) validCross[i]=0;

  //set up grid coordinates
  kx[0]=ky[0]=0;
  for (i=1;i<=width;i++) kx[i]=kx[i-1]+gridsize;
  for (i=1;i<=height;i++) ky[i]=ky[i-1]+gridsize;

  if (u.x!=v.x){
    uprime.x = 0; 
    uprime.y = u.y+ydiff*(-u.x/xdiff); 
    uprime.z = u.z+(v.z-u.z)*(-u.z)/xdiff;
    vprime.x = width*gridsize;
    vprime.y = u.y+ydiff*(width*gridsize-u.x)/xdiff;
    vprime.z = u.z+(v.z-u.z)*(width*gridsize-u.x)/xdiff;
  }
  else {//x coords are same so y coords are diff
    uprime.x = u.x + xdiff*(-u.y/ydiff);
    uprime.y = 0;
    uprime.z = u.z+(v.z-u.z)*(-u.y)/ydiff;
    vprime.x = u.x + xdiff*(height*gridsize-u.y)/ydiff; 
    vprime.y = height*gridsize; 
    vprime.z = u.z+(v.z-u.z)*(height*gridsize-u.y)/ydiff;
  }


  retval=wCrossings(uprime,vprime,0,kx,width+1,xlambda);//rgr
  if (retval==(-2)) 
    xlen=0;
  else 
    xlen=width+1;

  retval=wCrossings(uprime,vprime,1,ky,height+1,ylambda);//rgr
  if (retval==(-3)) 
    ylen=0;
  else 
    ylen=height+1;
  
#ifdef DEBUG
  if (xlen>0) for (i=0;i<=width;i++) printf("xlamba %.10f\n",xlambda[i]);
  if (ylen>0) for (i=0;i<=height;i++) printf("ylamba %.10f\n",ylambda[i]);
#endif

  xylen = merge(xlambda,ylambda,xlen,ylen,xylambda);

#ifdef DEBUG
  printf("merge len is %d\n",xylen);
  for (i=0;i<xylen;i++) printf("xylamba %.10f\n",xylambda[i]);
#endif

  computeCrossings(uprime,vprime,xylambda,xylen,crossings);
  //check for valid crossings
  for (i=0;i<xylen;i++)
    if ((crossings[i].x>=0)&&(crossings[i].x<=(width*gridsize))&&
	(crossings[i].y>=0)&&(crossings[i].y<=(height*gridsize)))
      validCross[i]=1;
    else
      validCross[i]=0;
#ifdef DEBUG
  for (i=0;i<xylen;i++) printf("crossing %.10f %.10f %.10f   * %d\n",
			       crossings[i].x,crossings[i].y,
			       crossings[i].z,validCross[i]);
#endif

  //remove invalid crossings
  retval=removeInvalid(crossings,validCross,height,width,xylen);
  xylen=retval;

  for (i=1;i<xylen;i++) 
    lengths[i-1]=sqrt( SQR(crossings[i].x-crossings[i-1].x)+
		       SQR(crossings[i].y-crossings[i-1].y)+
		       SQR(crossings[i].z-crossings[i-1].z));

  for (i=0;i<(xylen-1);i++)
    //    if (lengths[i]==0) validCross[i]=0;
      if (lengths[i]<1e-10) validCross[i]=0;

#ifdef DEBUG
  printf("Before culling number of  lengths is %d\n",xylen-1);
  for (i=0;i<(xylen-1);i++)
    printf("length %d is %e validity is %d\n",i,lengths[i],validCross[i]);
#endif


  retval=removeInvalid(crossings,validCross,height,width,xylen);
  xylen=retval;

  //now recompute lengths; all should be greater than epsilon
  for (i=1;i<xylen;i++) 
    lengths[i-1]=sqrt( SQR(crossings[i].x-crossings[i-1].x)+
		       SQR(crossings[i].y-crossings[i-1].y)+
		       SQR(crossings[i].z-crossings[i-1].z));

#ifdef DEBUG
  printf("After culling length is %d\n",xylen);
  for (i=0;i<xylen-1;i++) printf("crossing %.10f %.10f %.10f  %e  * %d\n",
			       crossings[i].x,crossings[i].y,
    crossings[i].z,lengths[i],validCross[i]);
  printf("crossing %.10f %.10f %.10f  * %d\n", crossings[i].x,crossings[i].y,crossings[i].z,validCross[i]);
#endif

  for (i=0;i<(xylen-1);i++){
    interiors[i].x=(crossings[i].x+crossings[i+1].x)/2;
    interiors[i].y=(crossings[i].y+crossings[i+1].y)/2;
    interiors[i].z=(crossings[i].z+crossings[i+1].z)/2;
  }
#ifdef DEBUG
  for (i=0;i<(xylen-1);i++){
    printf("%d   %f %f  %f\n",i,interiors[i].x,interiors[i].y,interiors[i].z);
  }
#endif
  for (i=0;i<(xylen-1);i++){
    segs[i].prow=ceil(interiors[i].x/gridsize);
    segs[i].pcol=ceil(interiors[i].y/gridsize);
    segs[i].length=lengths[i];
  }

#ifdef DEBUG
  printf("Segments\n");
  for (i=0;i<(xylen-1);i++){
    printf("%d   %d %d  %f\n",i,segs[i].prow,segs[i].pcol,segs[i].length);
  }
#endif  

  free(kx);
  free(ky);
  free(xlambda);
  free(ylambda);
  free(xylambda);
  free(lengths);
  free(crossings);
  free(interiors);
  free(validCross);

  return(xylen-1);
}//segLengths

void computeCrossings(u,v,lamb,len,cross)
struct point u,v,cross[];
double lamb[];
int len;
{
  int i;
  for (i=0; i<len; i++) {
    cross[i].x = u.x + lamb[i]*(v.x-u.x);
    cross[i].y = u.y + lamb[i]*(v.y-u.y);
    cross[i].z = u.z + lamb[i]*(v.z-u.z);
  }
}//computeCrossings

int removeInvalid(cross,valid,height,width,xylen)
struct point cross[];
int valid[],height,width,xylen;
{
  int i,idx,sum,total;
  struct point *temp;
#ifdef DEBUG
  printf("in REMOVE len is %d\n",xylen);
#endif

  total = height + width + 2;
  temp = calloc(total,sizeof(temp[0]));

  sum=0;
  for (i=0;i<xylen;i++) sum+=valid[i];
  if (sum==xylen) return(xylen);  //valid so return
#ifdef DEBUG
  printf("in REMOVE about to check %d\n",xylen);
#endif

  idx=0;
  for (i=0;i<(xylen);i++)
    if (valid[i]) {
      temp[idx].x=cross[i].x;
      temp[idx].y=cross[i].y;
      temp[idx].z=cross[i].z;
      idx++;
    }//valid
  for (i=0;i<idx;i++){
    cross[i].x=temp[i].x;
    cross[i].y=temp[i].y;
    cross[i].z=temp[i].z;
    valid[i]=1;
  }//i

  free(temp);

  return(idx);
}


int wCrossings(u,v,which,scale,len,lamb)
struct point u,v;
int    which,len;
double  scale[],lamb[];
{
  int i;
  int sz;
  int cmpfunc(); 

  if ((which>2)||(which<0)){
    printf("Calling wCrossings with incorrect coordinate %d\n",which);
    exit(2);
  }
  if ((v.x==u.x)&&(which==0)) return(-2);
  if ((v.y==u.y)&&(which==1)) return(-3);
  sz=sizeof(u.x);

  if (which==0){
    for (i=0;i<len;i++)      lamb[i] = (scale[i]-u.x)/(v.x-u.x);
    qsort(&lamb[0],len,sz,cmpfunc);
    return(0);
  }
  if (which==1){
    for (i=0;i<len;i++)      lamb[i] = (scale[i]-u.y)/(v.y-u.y);
    qsort(&lamb[0],len,sz,cmpfunc);
    return(0);
  }
  if (which==2){
    for (i=0;i<len;i++)      lamb[i] = (scale[i]-u.z)/(v.z-u.z);
    qsort(&lamb[0],len,sz,cmpfunc);
    return(0);
  }
}

void swap(a,b)
double *a,*b;
{
  double tmp;
  tmp = *b;
  (*b)=*a;
  (*a)=tmp;
}

int merge(a,b,m,n,sorted)
double a[],b[],sorted[];
int m,n;
{
  int i,j,k;
  j=k=0;

  if (m==0) {
    for (i=0;i<n;i++) sorted[i]=b[i];
    return(i);
  }
  if (n==0) {
    for (i=0;i<m;i++) sorted[i]=a[i];
    return(i);
  }

  i=0;
  while((j<m)||(k<n))
  {
    if (j<m && k<n){
      if (a[j]<b[k]){
	sorted[i]=a[j];
	j++;
      }//if a element < b element
      else{
	  sorted[i]=b[k];
	  if (a[j]==b[k])
	    j++;
	  k++;
      }//else  choosing b
      i++;

    }//if a and b both have guys to put in sorted
    else if ((j==m)&&(k<n)){
      for (;k<n;){
	sorted[i]=b[k];
	k++;
	i++;
      }
    }// else j==m
    else{
      if ((j<m)&&(k==n)){
	for (; j<m;){
	  sorted[i]=a[j];
	  j++;
	  i++;
	}
      }
    }//last else
  }//while
  return(i);
}

int cmpfunc( const void  *a, const void *b){
  const double *fa= (const double *)a;
  const double *fb= (const double *)b;
  if ( (*fa) > *fb) return(1);
  if ( (*fb) > *fa) return(-1);
  if ( (*fb)==(*fa)) return(0);
}
