#include <stdlib.h>
#include <stdio.h>
#include "/home/wang/Labwork/Software/mldaetlib/allocate.c"
#include "/home/wang/Labwork/Software/mldaetlib/lapack/lapacke.h"
#define N 200


void geteigen(double **matrix,int dim,double *eigenvalue_real,double *eigenvalue_image,double **leftmatrix)
 {
 int i,j;
 char JOBVL[1],JOBVR[1];
 int DIM[1],LDA[1],LDB[1],LDVL[1],LDVR[1],LWORK[1],INFO[1];
 double *A,*B,*ALPHAR,*ALPHAI,*BETA,*VL,*VR,*WORK;

 DIM[0]=dim;
 JOBVL[0]='V';
 JOBVR[0]='N';
 A=doublearray(DIM[0]*DIM[0]);
 LDA[0]=DIM[0];
 B=doublearray(DIM[0]*DIM[0]);
 LDB[0]=DIM[0];
 ALPHAR=doublearray(DIM[0]);
 ALPHAI=doublearray(DIM[0]);
 BETA=doublearray(DIM[0]);
 VL=doublearray(DIM[0]*DIM[0]);
 LDVL[0]=DIM[0];
 VR=doublearray(DIM[0]*DIM[0]);
 LDVR[0]=DIM[0];
 LWORK[0]=DIM[0]*8;
 WORK=doublearray(LWORK[0]);

 for(i=0;i<DIM[0];i++)
  {
  for(j=0;j<DIM[0];j++)
   {
   A[j*DIM[0]+i] = matrix[i][j];
   if(j==i)
    B[j*DIM[0]+i] = 1;
   else
    B[j*DIM[0]+i] = 0;
   }
  }

 LAPACK_dggev(JOBVL,JOBVR,DIM,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
 
// printf("%d\n",INFO[0]);
 for(i=0;i<DIM[0];i++)
  {
  eigenvalue_real[i] = ALPHAR[i];
  eigenvalue_image[i] = ALPHAI[i];
  for(j=0;j<DIM[0];j++)
   {
   leftmatrix[i][j] = VL[j+i*DIM[0]];
   }
  }
 

 free(A);
 free(B);
 free(ALPHAI);
 free(ALPHAR);
 free(BETA);
 free(VL);
 free(VR);
 free(WORK);


 }





void geteigen_sort_dec(double *eigenvalue,int dim,double **leftmatrix)
{
int i,j,maxid;
double tempdata,maxvalue;


//tempvector = doublearray(dim);

for(i=0;i<dim-1;i++)
 {
 maxvalue = eigenvalue[i];
 maxid = i;
//find the next maximum eigenvalue
 for(j=i+1;j<dim;j++)
  {
  if(eigenvalue[j] > maxvalue)
   {
   maxvalue = eigenvalue[j];
   maxid = j;
   }
  }
//swap eigenvalues and eigenvectors
 if(maxid != i)
  {
  eigenvalue[maxid] = eigenvalue[i];
  eigenvalue[i] = maxvalue;
  for(j=0;j<dim;j++)
   {
   tempdata = leftmatrix[i][j];
   leftmatrix[i][j] = leftmatrix[maxid][j];
   leftmatrix[maxid][j] = tempdata;
   }
  }
 }
}













void geteigen_sort_inc(double *eigenvalue,int dim,double **leftmatrix)
{
int i,j,minid;
double tempdata,minvalue;


//tempvector = doublearray(dim);

for(i=0;i<dim-1;i++)
 {
 minvalue = eigenvalue[i];
 minid = i;
//find the next maximum eigenvalue
 for(j=i+1;j<dim;j++)
  {
  if(eigenvalue[j] < minvalue)
   {
   minvalue = eigenvalue[j];
   minid = j;
   }
  }
//swap eigenvalues and eigenvectors
 if(minid != i)
  {
  eigenvalue[minid] = eigenvalue[i];
  eigenvalue[i] = minvalue;
  for(j=0;j<dim;j++)
   {
   tempdata = leftmatrix[i][j];
   leftmatrix[i][j] = leftmatrix[minid][j];
   leftmatrix[minid][j] = tempdata;
   }
  }
 }
}













