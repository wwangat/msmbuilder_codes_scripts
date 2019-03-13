#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "/home/group/code/c/mldaetlib/allocate.c"

char *getstrpart(char *string,int begin, int end)
{
int i,len;
char *output;

len = strlen(string);
if(len < end - begin)
 {
 printf("The string %s is not long enough to decompose\n");
 exit(0);
 }

output = chararray(end - begin);
for(i=0;i<end - begin;i++)
 output[i] = string[begin+i];

return(output);
}

