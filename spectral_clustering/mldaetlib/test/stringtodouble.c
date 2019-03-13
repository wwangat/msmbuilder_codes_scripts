#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200

void main()
{
char string[N];
double f;

printf("Input a number:\n");
scanf("%s",string);

f = getdoublefromstring(string,0,2);

printf("%lf\n",f);

//printf("%lf\n",stringtodouble(string));


}
