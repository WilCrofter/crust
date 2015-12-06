#include <stdio.h>
#include <stdlib.h>

main(){

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

  CgenS(S,&ht,&wid,&grid);
  for (i=0;i<(ht*ht);i++){
    printf("Row %d\n",i);
    for (j=0;j<(ht*wid);j++) 
      if (S[i][j]!=0)
	printf("%d %f\n",j,S[i][j]);
    printf("\n");
  }
}

fillit(arr,nrow,ncol)
int arr[],nrow,ncol;
{ 
  int i,j;
  for (i=0;i<nrow;i++){
    for (j=0;j<ncol;j++)      arr[(i*ncol)+j]=(i*ncol)+j;
    }
}
