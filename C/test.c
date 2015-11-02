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

  //nrow=ht^2
  //ncol=ht*wd;

  double S[4][8];
  int i,j,ht,wid,temp[2];
  double grid;
  double nonzero[32];
  int  rowind[32],colind[9];
  int icmpfunc();

  ht=2;
  wid=4;
  grid=1.;

  /* test qsort - make sure compare function compares correct type!!!
  temp[0]=3;
  temp[1]=0;
  for (i=0;i<2;i++) printf("%d %d\n",i,temp[i]);
  qsort(&temp[0],2,sizeof(temp[0]),icmpfunc);
  for (i=0;i<2;i++) printf("%d %d\n",i,temp[i]);
  exit(1);
  */
  /*
  CgenS(S,&ht,&wid,&grid);
  for (i=0;i<(ht*ht);i++){
    printf("Row %d\n",i);
    for (j=0;j<(ht*wid);j++) 
      if (S[i][j]!=0)
	printf("%d %f\n",j,S[i][j]);
    printf("\n");
  }
  */

  CgenS_sparse(&nonzero,&colind,&rowind,&ht,&wid,&grid);
}

fillit(arr,nrow,ncol)
int arr[],nrow,ncol;
{ 
  int i,j;
  for (i=0;i<nrow;i++){
    for (j=0;j<ncol;j++)      arr[(i*ncol)+j]=(i*ncol)+j;
    }
}
