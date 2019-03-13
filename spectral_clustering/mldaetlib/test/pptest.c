#include <stdio.h>

void main()
{
int i,j;

for(i=0,j=0;i<10;)
 {
 j=i++;
 printf("%d\n",j);
 }

printf("\n");
for(i=0,j=0;i<10;)
 {
 j=++i;
 printf("%d\n",j);
 }
}
