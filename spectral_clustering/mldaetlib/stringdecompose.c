
int stringlen(char *string)
{
int i;
for(i=0;string[i] != '\0';i++);
return i;

}



void stringcopy(char *oristring,char *tarstring)
{
int i;

for(i=0;oristring[i] != '\0';i++)
 {
 tarstring[i] = oristring[i];
 }

tarstring[i] = oristring[i];

}


int stringtoint(char *string)
{
int len,i;
int value;

for(len=0;string[len]!='\0';len++);

//preval=0;
value=0;
for(i=0;i<len;i++)
 {
 if((string[i] > 57) || (string[i] < 48))
  continue;
 value =value*10+(int)(string[i])-48;
 }
return value;
};


double stringtodouble(char *string)
{
double num,power;
int intlen,declen,label,i,sign,begin;
//printf("%s\n",string);

for(intlen=0,declen=0,label=0,i=0,sign=1,begin=0;string[i]!='\0';i++)
 {
 if(string[i] == ' ')
  {
  begin ++;
  continue;
  }
 if(string[i] == '-')
  {
  sign = -1;
  begin ++;
  continue;
  }
 if(string[i] == '.')
  {
  label=1;
  continue;
  }
 if(label==0)
  intlen++;
 if(label==1)
  declen++;
// printf("i=%d,%c,%d,%d\n",i,string[i],intlen,declen);
 }

num=0;

for(i=0;i<intlen;i++)
 {
 num=num*10+(int)(string[begin+i])-48;
// printf("%c\n",string[begin+i]);
// printf("%lf\n",num);
 }
//printf("%lf\n",num);

for(i=0,power=0.1;i<declen;i++)
 {
 num += power*((int)string[begin+intlen+1+i]-48);
 power*=0.1;
// printf("%lf\n",num);
 }
//printf("%s\n",string);

return num*sign;

};



float stringtofloat(char *string)
{
float num,power;
int intlen,declen,label,i,begin,sign;

sign = 1;	//for the positive or negative sign of data
begin = 0;	//begin = 1 means the it is not the space

for(intlen=0,declen=0,label=0,i=0;string[i]!='\0';i++)
 {
 if(string[i] == ' ')		//Get rid of the space before the sign label and number
  {
  begin ++;
  continue;
  }

 if(string[i] == '-')		//checking the sign of the data
  {
  sign = -1;
  begin ++;
  continue;
  }

  if(string[i] == '.')
   {
   label=1;
   continue;
   }
  if(label==0)
   intlen++;
  if(label==1)
   declen++;

 }

num=0;


for(i=0;i<intlen;i++)
 num=num*10+(int)(string[begin+i])-48;
//printf("%lf\n",num);

for(i=0,power=0.1;i<declen;i++)
 {
 num += power*((int)string[begin+intlen+1+i]-48);
 power*=0.1;
// printf("%lf\n",num);
 }

num = num*sign;

return num;

};


char *getstrpart(char *string,int begin, int end)
{
int i,len;
char *output;
char *chararray(int length);

len = stringlen(string);
if(len < end - begin)
 {
 printf("The string %s is not long enough to decompose\n",string);
 exit(0);
 }

output = chararray(end - begin+1);
for(i=0;i<end - begin;i++)
 output[i] = string[begin+i];
output[end-begin+1] = '\0';

return(output);
}


void getstrpartin(char *template,int begin,int end, char *target)
 {
 int i,j;
 for(i=0;i<end-begin;i++)
  target[i] = template[i+begin];
 target[i] = '\0';
 }





void exstrpart(char *string, int begin, int end, char *output)
{
int i,len;

len = stringlen(string);
if(len < end - begin)
 {
 printf("The string %s is not long enough to decompose\n",string);
 exit(0);
 }

//output = chararray(end - begin+1);
for(i=0;i<end - begin;i++)
 output[i] = string[begin+i];
output[end-begin] = '\0';

//return(output);
}




float getfloatfromstring(char *string, int begin, int end)
 {
 int i,len;
 char *output;
 float f;
 char *chararray(int length);

 len = stringlen(string);
 if(len < end - begin)
  {
  printf("The string %s is not long enough to decompose\n",string);
  exit(0);
  }
 
 output = chararray(end - begin+1);
 for(i=0;i<end - begin;i++)
  output[i] = string[begin+i];
 output[end-begin+1] = '\0';

 f = stringtofloat(output);
 free(output);
 return(f);
 }





int getlastintfromstring(char *string)
 {
int len,i,begin,end;
int value;

len = stringlen(string);
//printf("len:%d\n",len);
 for(i= len -1;(string[i] > 57) || (string[i] < 48);i--);
end = i+1;
//printf("end:%d\n",end);

 for(i = end-1;(string[i] <= 57) && (string[i] >= 48);i--);
begin = i+1;
//printf("begin:%d\n",begin);
//preval=0;
value=0;
for(i=begin;i<end;i++)
 {
 if((string[i] > 57) || (string[i] < 48))
  continue;
 value =value*10+(int)(string[i])-48;
 }
return value;
 }






double getdoublefromstring(char *string, int begin, int end)
 {
 int i,len;
 char *output;
 double f;
// double num,power;
// int i,len,intlen,declen,label,sign,beginpoint;
//
// len = stringlen(string);
////printf("%d\n",len);
// if(len < end)
//  {
//  printf("The string %s is not long enough to decompose\n",string);
//  exit(0);
//  }
////printf("%d\n",len);
////printf("%s\n",string);
//
//for(intlen=0,declen=0,label=0,i=0,sign=1,beginpoint=begin;i<end-begin;i++)
// {
// if(string[begin+i] == ' ')
//  {
//  beginpoint ++;
//  continue;
//  }
// if(string[begin+i] == '-')
//  {
//  sign = -1;
//  begin ++;
//  continue;
//  }
// if(string[begin+i] == '.')
//  {
//  label=1;
//  continue;
//  }
// if(label==0)
//  intlen++;
// if(label==1)
//  declen++;
//// printf("i=%d,%c,%d,%d\n",i,string[i],intlen,declen);
// }
//
//num=0;
//
//for(i=0;i<intlen;i++)
// {
// num=num*10+(int)(string[beginpoint+i])-48;
//// printf("%c\n",string[begin+i]);
//// printf("%lf\n",num);
// }
////printf("%lf\n",num);
//
//for(i=0,power=0.1;i<declen;i++)
// {
// num += power*((int)string[begin+intlen+1+i]-48);
// power*=0.1;
//// printf("%lf\n",num);
// }
////printf("%s\n",string);
//
//return num*sign;


output = (char *)malloc(sizeof(char)*(end-begin+1));
 
// output = chararray(end - begin+1);
//printf("%s\n",output);
 for(i=0;i<end - begin;i++)
  output[i] = string[begin+i];
 output[end-begin] = '\0';

//printf("%s\n",output);

 f = stringtodouble(output);
 free(output);
 return(f);
 }




int getintfromstring(char *string, int begin, int end)
 {
 int i,len,value;
 char *output;
 char *chararray(int length);

 len = stringlen(string);
 if(len < end - begin)
  {
  printf("The string %s is not long enough to decompose\n",string);
  exit(0);
  }

 output = chararray(end - begin+1);
 for(i=0;i<end - begin;i++)
  output[i] = string[begin+i];
// printf("%s\n",output);
 value = stringtoint(output);
 free(output);
 return(value);

 }




