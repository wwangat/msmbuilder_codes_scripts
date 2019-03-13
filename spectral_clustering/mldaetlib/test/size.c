#include <stdio.h>
#include <string.h>
#define N 200


void main()
{
char string[N];
int len;

printf("Input a string:\n");
scanf("%s",string);

len = strlen(string);

//printf("The lengthe of the string:%d\n",length(string));
printf("The lengthe of the string:%d\n",strlen(string));
printf("The lengthe of the string:%d\n",len);
}
