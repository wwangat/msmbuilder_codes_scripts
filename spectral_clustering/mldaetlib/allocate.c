#ifndef ALLOCATE_H
#define ALLOCATE_H

#include <string.h>
#include <stdlib.h>
#include "/home/wang/Labwork/Software/mldaetlib/arguement.c"

float *floatarray(int length)
{
float *array;
if(!(array = (float *)malloc(sizeof(float)*length)))
 {
 printf("Failed to allocate the memory for the float array:\n");
 exit(0);
 }
return array;
};




float **floatmatrix(int row,int col)
{
float **matrix;
int i;
if(!(matrix = (float **)malloc(sizeof(float *)*row)))
 {
 printf("Failed to allocate the memory for the matrix\n");
 exit(0);
 }

for(i=0;i<row;i++)
 {
 if(!(matrix[i] = (float *)malloc(sizeof(float)*col)))
  {
  printf("Failed to allocate the memory for the matrix\n");
  exit(0);
  }
 }
return matrix;
}




float ***floatmatrixarray(int row1, int row2, int row3)
{
float ***matrixarray;
int i;
if(!(matrixarray = (float ***)malloc(sizeof(float **)*row1)))
 {
 printf("Failed to allocate the memory for the float matrixarray\n");
 exit(0);
 }
for(i=0;i<row1;i++)
 matrixarray[i] = floatmatrix(row2,row3);
return matrixarray;
}




char **namearray(int filenum,int length)
{
 char **namelist;
 int i;
 if(!(namelist = (char **)malloc(sizeof(char *)*filenum)))
  {
  printf("Failed to allocate the memory for the file name list1\n");
  exit(0);
  }

 for(i=0;i<filenum;i++)
  {
  if(!(namelist[i] = (char *)malloc(sizeof(char)*length)))
   {
   printf("Failed to allocate the memory for the file name list2\n");
   exit(0);
   }
  } 
return namelist;
}



char *chararray(int length)
{
char *string;
if(!(string=(char *)malloc(sizeof(char)*length)))
 {
 printf("Failed to allocate the memory for the string\n");
 exit(0);
 }
return string;
}



FILE **filearray(int num)
{
FILE **temp;
if(!(temp=(FILE **)malloc(sizeof(FILE *)*num)))
 {
 printf("Failed to allocate the memory for the filearray\n");
 exit(0);
 }
return temp;
}




int *intarray(int length)
{
int *array;
if(!(array = (int *)malloc(sizeof(int)*(length+1))))
 {
 printf("Failed to allocate the memory for the int array:\n");
 exit(0);
 }
return array;
}



int **intmatrix(int row,int col)
{
int **matrix;
int i;
if(!(matrix = (int **)malloc((row+1)*sizeof(int *))))
 {
 printf("Failed to allocate the memory for the int\n");
 exit(0);
 }
matrix[row] = NULL;

for(i=0;i<row;i++)
 {
 if(!(matrix[i] = (int *)malloc(sizeof(int)*(col+1))))
  {
  printf("Failed to allocate the memory for the int matrix\n");
  exit(0);
  }
 }
return matrix;
}




double *doublearray(int length)
{
double *array;
if(!(array = (double *)malloc(sizeof(double)*length)))
 {
 printf("Failed to allocate the memory for the double array:\n");
 exit(0);
 }
return array;
}




double **doublematrix(int row,int col)
{
double **matrix;
int i;
if(!(matrix = (double **)malloc(sizeof(double *)*row)))
 {
 printf("Failed to allocate the memory for the double\n");
 exit(0);
 }
for(i=0;i<row;i++)
 {
 if(!(matrix[i] = (double *)malloc(sizeof(double)*col)))
  {
  printf("Failed to allocate the memory for the double matrix\n");
  exit(0);
  }
 }
return matrix;

}




double ***doublematrixarray(int row1, int row2, int row3)
{
double ***matrixarray;
int i;
if(!(matrixarray = (double ***)malloc(sizeof(double **)*row1)))
 {
 printf("Failed to allocate the memory for the double matrix array");
 exit(0);
 }

for(i=0;i<row1;i++)
 matrixarray[i] = doublematrix(row2,row3);

return matrixarray;
}





int ***intmatrixarray(int row1, int row2, int row3)
{
int ***matrixarray;
int i;
if(!(matrixarray = (int ***)malloc(sizeof(int **)*row1)))
 {
 printf("Failed to allocate the memory for the int matrixarray\n");
 exit(0);
 }
for(i=0;i<row1;i++)
 matrixarray[i] = intmatrix(row2,row3);
return matrixarray;
}





FILE **filematrix(int num)
{
FILE **matrix;
int i;
if(!(matrix = (FILE **)malloc(sizeof(FILE *)*num)))
 {
 printf("Failed to allocate the memory for the FILE\n");
 exit(0);
 }

for(i=0;i<num;i++)
 {
 if(!(matrix[i] = (FILE *)malloc(sizeof(FILE))))
  {
  printf("Failed to allocate the memory for the FILE matrix\n");
  exit(0);
  }
 }
return matrix;
};



struct arguement *arguearray(int num)
{
int i,j;
struct arguement *arguelist;
arguelist=(struct arguement *)malloc(sizeof(struct arguement)*num);
for(i=0;i<num;i++)
 {
 arguelist[i].arguename=chararray(200);
 arguelist[i].assigned=0;
 arguelist[i].parameter.s=chararray(200);
 }
return arguelist;
};


double ****doublematrixmatrix(int row1, int row2, int row3, int row4)
{
double ****matrixmatrix;
int i;
if(!(matrixmatrix = (double ****)malloc(sizeof(double ***)*row1)))
 {
 printf("Failed to allocate the memory for the double matrix array");
 exit(0);
 }

for(i=0;i<row1;i++)
 matrixmatrix[i] = doublematrixarray(row2,row3,row4);

return matrixmatrix;
}






#endif

