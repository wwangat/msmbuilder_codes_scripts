#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/read.c"
#define N 200

void main()
{
char line[N];
int linenum;
FILE *file;

file = openfile("test.dat",'r');

for(linenum = 0; ; linenum ++)
 {
 if(fgets(line,N,file) == NULL)
  break;
 printf("%d:%s",linenum,line);
 }

printf("%d\n",linenum);

fclose(file);
}



