#include <stdio.h>
#include <stdlib.h>

FILE *openfile(char *filename,char type)
{
FILE *file;
if(type == 'r')
 {
 if((file = fopen(filename,"r")) == NULL)
  {
  printf("Failed to open file %s\n",filename);
  exit(0);
  }
 }

if(type == 'w')
 {
 if((file = fopen(filename,"w")) == NULL)
  {
  printf("Failed to build file %s\n",filename);
  exit(0);
  }
 }

return(file);
}

/*void test(char type)
{
if(type == 'r')
 printf("it is %c\n",type);
else
 printf("it is not %c\n",type);
}
*/



void readintarray(char *filename,int *array,int dim)
 {
 int i;
 FILE *file;
 file = openfile(filename,'r');
 for(i=0;i<dim;i++)
  fscanf(file,"%d",&array[i]);
 fclose(file);
 }


void readintmatrix(char *filename,int **matrix,int rownum,int colnum)
 {
 int i,j;
 FILE *file;
 file = openfile(filename,'r');
 for(i=0;i<rownum;i++)
  for(j=0;j<colnum;j++)
   fscanf(file,"%d",&matrix[i][j]);
 fclose(file);
 }



void readdoublearray(char *filename,double *array,int dim)
 {
 int i;
 FILE *file;
 file = openfile(filename,'r');
 for(i=0;i<dim;i++)
  fscanf(file,"%lf",&array[i]);
 fclose(file);
 }


void readdoublematrix(char *filename,double **matrix,int rownum,int colnum)
 {
 int i,j;
 FILE *file;
 file = openfile(filename,'r');
 for(i=0;i<rownum;i++)
  for(j=0;j<colnum;j++)
   fscanf(file,"%lf",&matrix[i][j]);
 fclose(file);
 }
