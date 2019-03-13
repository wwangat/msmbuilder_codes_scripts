#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include <stdio.h>
#include <stdlib.h>
#define N 100

void main()
{
char name='r';
float *array;
int n=5,i;

array = floatarray(n);

for(i=0;i<n;i++)
 {
 array[i] = i*i;
 printf("%d:  %f\n",i,array[i]);
 }
}
