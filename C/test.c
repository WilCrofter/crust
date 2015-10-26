#include <stdio.h>

main(){
  /*  int arr[4][6];
  int i,j;
  fillit(arr,4,6);
  for (i=0;i<4;i++){
    for (j=0;j<6;j++) printf("%d ",arr[i][j]);
    printf("\n");
  }
  */
  double S[64][96];
  int i,j;
  CgenS(S,8,12,1.);
  for (i=0;i<64;i++){
    printf("Row %d\n",i);
    for (j=0;j<96;j++) 
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
