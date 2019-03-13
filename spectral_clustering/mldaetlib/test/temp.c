#include<stdio.h>
#include<stdlib.h>
#include"/home/group/code/c/mldaetlib/geteigen.c"
#include"/home/group/code/c/matrix/inverse/matrixinverse.c"
//#include"/home/group/code/c/mldaetlib/allocate.c"

void main()
{
int i,j,dim;
double **matrix,*eigenr,*eigeni,**leftmatrix;

dim = 2;
matrix = doublematrix(dim,dim);
eigenr = doublearray(dim);
eigeni = doublearray(dim);
leftmatrix = doublematrix(dim,dim);


matrix[0][0] = 0.8;
matrix[0][1] = 0.2;
matrix[1][0] = 0.4;
matrix[1][1] = 0.6;

geteigen(matrix,dim,eigenr,eigeni,leftmatrix);
matrixinverse(leftmatrix,dim,matrix);

for(i=0;i<dim;i++)
 printf("%lf  %lf\n",eigenr[i],eigeni[i]);

for(i=0;i<dim;i++)
 printf("%lf  %lf\n",leftmatrix[i][0],leftmatrix[i][1]);

for(i=0;i<dim;i++)
 printf("%lf  %lf\n",matrix[i][0],matrix[i][1]);

for(i=0;i<dim;i++)
 {
 free(matrix[i]);
 free(leftmatrix[i]);
 }
free(matrix);
free(leftmatrix);
free(eigenr);
free(eigeni);

}

