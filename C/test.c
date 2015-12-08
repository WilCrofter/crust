#include <stdio.h>
#include <stdlib.h>

main(){
  /*  int arr[4][6];
  int i,j;
  fillit(arr,4,6);
  for (i=0;i<4;i++){
    for (j=0;j<6;j++) printf("%d ",arr[i][j]);
    printf("\n");
  }
  */

  //nrow=ht^2
  //ncol=ht*wd;

  double S[16][64];  //row is ht^2, col is ht*wid, see assignments below
  int i,j,ht,wid,temp[2];
  double grid;
  double nonzero[163840];  //max number of nonzero entries in sparse matrix
  int  colind[641]; //colind is 1+number of columns in matrix
  int rowind[163840];   //max number of nonzero entries in sparse matrix
  int rowct[256],nz;

  ht=4;
  wid=16;
  grid=1.84;

  /* test qsort - make sure compare function compares correct type!!!
  temp[0]=3;
  temp[1]=0;
  for (i=0;i<2;i++) printf("%d %d\n",i,temp[i]);
  qsort(&temp[0],2,sizeof(temp[0]),icmpfunc);
  for (i=0;i<2;i++) printf("%d %d\n",i,temp[i]);
  exit(1);
  */

  CgenS(S,&ht,&wid,&grid);
  for (i=0;i<(ht*ht);i++){
    printf("Row %d\n",i);
    for (j=0;j<(ht*wid);j++) 
      if (S[i][j]!=0)
	printf("%d %f\n",j,S[i][j]);
    printf("\n");
  }

  /*
   CgenS_sparse(&nonzero,&colind,&rowind,&ht,&wid,&grid);
   i=0;
   while ((nonzero[i]>0) && (i<163840)) {
     rowct[rowind[i]]++;
     if (rowind[i]==55) 
       printf("%f\n",nonzero[i]);
     i++;
   }//while
   if (i>163840) {
     printf("Uh Oh - too many nonzero guys!!!\n");
     exit(1);
   }
   nz=i;
   printf("Number of nonzero guys is %d\n",nz);
   for (i=0;i<256;i++) printf("%3d %d\n",i,rowct[i]);
  */
}

fillit(arr,nrow,ncol)
int arr[],nrow,ncol;
{ 
  int i,j;
  for (i=0;i<nrow;i++){
    for (j=0;j<ncol;j++)      arr[(i*ncol)+j]=(i*ncol)+j;
    }
}
