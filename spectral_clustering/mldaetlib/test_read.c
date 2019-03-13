#include "read.c"
#include <stdio.h>
#include <stdlib.h>
#define N 100

void main()
{
char filename[N];
int i,endfile;
FILE *file;

printf("Input the name for the input file:\n");
scanf("%s",filename);

file = openfile(filename,'r');

for(;!feof(file);)
 {
 endfile = fscanf(file,"%d",&i);
 if(endfile != 1)
  break;
 printf("%d\n",i);
 }


fclose(file);
}

