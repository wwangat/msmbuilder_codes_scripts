#include <string.h>
#include <stdlib.h>
#include "/home/wang/Labwork/Software/mldaetlib/stringdecompose.c"
#define INT 0
#define FLOAT 1
#define DOUBLE 2
#define STRING 3




//Define the data type
enum data_type
{
int_type,float_type,double_type,string_type
};



union data
{
int i;
float f;
double d;
char *s;
};



struct arguement
{
char *arguename;
enum data_type type;
union data parameter;
int assigned;
};

void argueinitialize(struct arguement *arguearray,int seq,char *arguename,int datatype)
{
arguearray[seq].arguename=arguename;

arguearray[seq].type=datatype;
arguearray[seq].assigned=0;
};






//Do the assignment for the arguement parameters
void argueassign(struct arguement *arguearray,int argnum,char *paraname,char *parameter)
{
int i,seq,check=0;
int stringtoint(char *string);
float stringtofloat(char *string);
double stringtodouble(char *string);
//printf("%s=%s,%d\n",arguearray[0].arguename,arguearray[0].parameter.s,arguearray[0].type);
for(i=0;i<argnum;i++)
 {
 if(strcmp(arguearray[i].arguename,paraname) == 0)
  {
  check = 1;
  seq = i;
  break;
  }
 }


if(check == 0)
 {
 printf("%s not founded in the parameter list\n",paraname);
 exit(0);
 }

//printf("%s=%s,%d\n",arguearray[0].arguename,arguearray[0].parameter.s,arguearray[0].type);

    if(arguearray[seq].type == int_type)
     arguearray[seq].parameter.i = stringtoint(parameter);
  else  if(arguearray[seq].type == float_type)
     arguearray[seq].parameter.f = stringtofloat(parameter);
  else  if(arguearray[seq].type == double_type)
     arguearray[seq].parameter.d = stringtodouble(parameter);
  else  if(arguearray[seq].type == string_type)
     strcpy(arguearray[seq].parameter.s,parameter);
    else
     {
     printf("The type of the parameter %s is not pre-determined\n",paraname);
     exit(0);
     }
/*
for(i=0;i<argnum;i++)
{
if(arguearray[i].type==string_type)
printf("%s=%s,%d\n",arguearray[i].arguename,arguearray[i].parameter.s,arguearray[i].type);
if(arguearray[i].type==int_type)
printf("%s=%d,%d\n",arguearray[i].arguename,arguearray[i].parameter.i,arguearray[i].type);
}
printf("\n");
*/
arguearray[seq].assigned=1;
};




void argueallassigned(struct arguement *arguearray,int argnum)
{
int i;
for(i=0;i<argnum;i++)
 {
 if(arguearray[i].assigned == 0)
  {
  printf("The %s has not been assigned\n",arguearray[i].arguename);
  exit(0);
  }
 }
};



//Decompose the lines that from the file to the parameter that we need
int linetostrings(char *inputline,char *firststring,char *secondstring)
{
//The return value is one or zero. One means that this line is not a comment line, ie not begins with ";". Zero means it is a comment line

int i,j,k;
int check;
check =0;
//firststring="";
//secondstring="";

for(i=0,j=0,k=0;i<200;i++)
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
};


