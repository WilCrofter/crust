#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG1
#define HEIGHT 8  //y dim of image
#define WIDTH  12  //x dim of image

#define SQHT (HEIGHT*HEIGHT)
#define NPIX (HEIGHT*WIDTH)
#define HEIGHT1 (HEIGHT+1)
#define WIDTH1  (WIDTH+1)
#define TOTAL   (HEIGHT1+WIDTH1)
#define SQR(A)  ((A)*(A))
struct point{double x; double y; double z;};
struct ipt{int x; int y; int z;};
struct segment{int prow, pcol; double length;};


void oldCgenS(S,gridsize)
double gridsize,S[][NPIX];
{
  int i,j,k,ridx,cidx,retval;
  int segLengths();
  struct point u,v;
  struct segment segs[TOTAL];
#ifdef DEBUG
  printf("Hello height is %d width is %d grid is %f\n",HEIGHT,WIDTH,gridsize);
#endif
  u.x = -1e-6;
  v.x = (WIDTH*gridsize)+1e-6;
  u.z = v.z = 0;
  for (i=0;i<SQHT;i++)
    for (j=0;j<NPIX;j++)      S[i][j]=0.0;

  for (i=1;i<=HEIGHT;i++){ //for each transmitter
    u.y = (i-0.5)*gridsize;
    for (j=1;j<=HEIGHT;j++){ //for each receiver
      v.y = (j-0.5)*gridsize;
      retval = segLengths(u,v,gridsize,&segs);
#ifdef DEBUG
      printf("retval is %d\n",retval);
#endif
      ridx=(j-1)*HEIGHT + i - 1;
      for (k=0;k<retval;k++){
	cidx=(segs[k].pcol-1)*WIDTH + segs[k].prow -1;
	S[ridx][cidx]=segs[k].length;
      }//

    }//j
  }//i
}

int segLengths(u,v,gridsize,segs)
     struct point u,v;
     double gridsize;
     struct segment segs[];
{
  double        i1[2],i2[2],kx[WIDTH1],ky[HEIGHT1];
  double        xdiff,ydiff,xlambda[WIDTH1],ylambda[HEIGHT1];
  double        xylambda[TOTAL],lengths[TOTAL];
  int           i,xylen,xlen,ylen,retval;
  int           wCrossings(),removeInvalid(),merge();
  struct point uprime,vprime,crossings[TOTAL],interiors[TOTAL];
  int          validCross[TOTAL];
  void         computeCrossings(),swap();

#ifdef DEBUG
  printf("in segLengths\nu is %f %f %f\nv is %f %f %f\n",
	 u.x,u.y,u.z,v.x,v.y,v.z);
#endif

  xdiff = v.x-u.x;
  ydiff = v.y-u.y;

  //make sure points are distinct
  if ((u.x==v.x)&&(u.y==v.y)) {return(1);}

  //make sure line connecting points crosses grid interior
  i1[0] = (-u.x)/xdiff;
  i1[1] = (WIDTH*gridsize-u.x)/xdiff;
  if (i1[1]<i1[0]) swap(&i1[0],&i1[1]);

  i2[0] = (-u.y)/ydiff;
  i2[1] = (HEIGHT*gridsize-u.y)/ydiff;
  if (i2[1]<i2[0]) swap(&i2[0],&i2[1]);

#ifdef DEBUG
  printf("i1 is %f %f   i2 is %f %f\n",i1[0],i1[1],i2[0],i2[1]);
#endif

  if ((i1[0]>=i2[1]) || (i2[0]>=i1[1])) {return(2);}

  for (i=0;i<TOTAL;i++) validCross[i]=0;

  //set up grid coordinates
  kx[0]=ky[0]=0;
  for (i=1;i<=WIDTH;i++) kx[i]=kx[i-1]+gridsize;
  for (i=1;i<=HEIGHT;i++) ky[i]=ky[i-1]+gridsize;

  if (u.x!=v.x){
    uprime.x = 0; 
    uprime.y = u.y+ydiff*(-u.x/xdiff); 
    uprime.z = u.z+(v.z-u.z)*(-u.z)/xdiff;
    vprime.x = WIDTH*gridsize;
    vprime.y = u.y+ydiff*(WIDTH*gridsize-u.x)/xdiff;
    vprime.z = u.z+(v.z-u.z)*(WIDTH*gridsize-u.x)/xdiff;
  }
  else {//x coords are same so y coords are diff
    uprime.x = u.x + xdiff*(-u.y/ydiff);
    uprime.y = 0;
    uprime.z = u.z+(v.z-u.z)*(-u.y)/ydiff;
    vprime.x = u.x + xdiff*(HEIGHT*gridsize-u.y)/ydiff; 
    vprime.y = HEIGHT*gridsize; 
    vprime.z = u.z+(v.z-u.z)*(HEIGHT*gridsize-u.y)/ydiff;
  }


  retval=wCrossings(uprime,vprime,0,kx,WIDTH1,&xlambda);
  if (retval==(-2)) 
    xlen=0;
  else 
    xlen=WIDTH1;

  retval=wCrossings(uprime,vprime,1,ky,HEIGHT1,&ylambda);
  if (retval==(-3)) 
    ylen=0;
  else 
    ylen=HEIGHT1;
  
#ifdef DEBUG
  if (xlen>0) for (i=0;i<=WIDTH;i++) printf("xlamba %.10f\n",xlambda[i]);
  if (ylen>0) for (i=0;i<=HEIGHT;i++) printf("ylamba %.10f\n",ylambda[i]);
#endif

  xylen = merge(xlambda,ylambda,xlen,ylen,xylambda);

#ifdef DEBUG
  printf("merge len is %d\n",xylen);
  for (i=0;i<xylen;i++) printf("xylamba %.10f\n",xylambda[i]);
#endif

  computeCrossings(uprime,vprime,xylambda,xylen,&crossings[0]);
  //check for valid crossings
  for (i=0;i<xylen;i++)
    if ((crossings[i].x>=0)&&(crossings[i].x<=(WIDTH*gridsize))&&
	(crossings[i].y>=0)&&(crossings[i].y<=(HEIGHT*gridsize)))
      validCross[i]=1;
    else
      validCross[i]=0;
#ifdef DEBUG
  for (i=0;i<xylen;i++) printf("crossing %.10f %.10f %.10f   * %d\n",
			       crossings[i].x,crossings[i].y,
			       crossings[i].z,validCross[i]);
#endif

  //remove invalid crossings
  retval=removeInvalid(&crossings[0],&validCross[0],xylen);
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


  retval=removeInvalid(&crossings[0],&validCross[0],xylen);
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
  return(xylen-1);
}//gens

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

int removeInvalid(cross,valid,xylen)
struct point cross[];
int valid[],xylen;
{
  int i,idx,sum;
  struct point temp[TOTAL];
#ifdef DEBUG
  printf("in REMOVE len is %d\n",xylen);
#endif

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
/*
void main(){
  double S[64][96];
  int i,j;
  oldCgenS(S,8,12,.5);
  for (i=0;i<64;i++){
    printf("Row %d\n",i);
    for (j=0;j<96;j++) 
      if (S[i][j]!=0)
	printf("%d %f\n",j,S[i][j]);
    printf("\n");
  }

}
*/
