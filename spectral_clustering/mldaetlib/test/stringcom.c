#include <stdio.h>
#include <string.h>
#include "../allocate.c"

void main()
{
char *string1,*string2;

string1=chararray(100);
string2=chararray(100);

scanf("%s",string1);
scanf("%s",string2);

//printf("%d\n",strcmp(string1,string2));
strcpy(string1,string2);
printf("%s\n%s\n",string1,string2);

string2="test";

printf("%s\n%s\n",string1,string2);

}
