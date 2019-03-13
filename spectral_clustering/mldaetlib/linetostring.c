#include <stdio.h>
#define N 200
#include <string.h>

int linetostrings(char *inputline,char *firststring,char *secondstring)
{
//The return value is one or zero. One means that this line is not a comment line, ie not begins with ";". Zero means it is a comment line

int i,j,k;
int check;
check =0;
//firststring="";
//secondstring="";

for(i=0,j=0,k=0;i<N;i++)
 {
 if((inputline[i] == ' ')||(inputline[i] == '\t'))
  continue;
 else if((inputline[i] == '\0')||(inputline[i] == ';')||(inputline[i] == '\n'))
  break;
 else
 {
 check = 1;
 if(inputline[i] == '=')
  {
   k=1;
   j=0;
  }
  else
  {
  if(k==0)
   {
   firststring[j]=inputline[i];
   j++;
   firststring[j]='\0';
   }
  if(k==1)
   {
   secondstring[j]=inputline[i];
   j++;
   secondstring[j]='\0';
   }
  }
 }
 }
 
return check;
}
