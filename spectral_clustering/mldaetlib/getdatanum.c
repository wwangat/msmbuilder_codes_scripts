#include <stdlib.h>
#include <stdio.h>
#define M 3000

int getintnum(char *popfilename) //Get the number of integer in the file
 {
 int statenum,eof;
 FILE *popfile;
 int data;

 popfile = openfile(popfilename,'r');
 for(statenum=0;;statenum++)
  if((eof = fscanf(popfile,"%d",&data)) == -1)
   break;

  fclose(popfile);
 return(statenum);
 }



int getnamenum(char *listfilename)      //Get the number of file names in the filelist
 {
 int filenum,eof;
 FILE *listfile;
 char filename[M];

 listfile = openfile(listfilename,'r');
 for(filenum=0;;filenum++)
  if((eof = fscanf(listfile,"%s",filename)) == -1)
   break;

 fclose(listfile);
 return(filenum);
 }




int getfloatnum(char *popfilename) //Get the number of float data in the file
 {
 int statenum,eof;
 FILE *popfile;
 float data;

 popfile = openfile(popfilename,'r');
 for(statenum=0;;statenum++)
  if((eof = fscanf(popfile,"%f",&data)) == -1)
   break;

  fclose(popfile);
 return(statenum);
 }




int getdoublenum(char *popfilename) //Get the number of float data in the file
 {
 int statenum,eof;
 FILE *popfile;
 double data;

 popfile = openfile(popfilename,'r');
 for(statenum=0;;statenum++)
  if((eof = fscanf(popfile,"%lf",&data)) == -1)
   break;

  fclose(popfile);
 return(statenum);
 }




int getlinenum(char *inputfilename)
{
char line[M];
int linenum;
FILE *inputfile;

inputfile = openfile(inputfilename,'r');

for(linenum = 0; ; linenum ++)
 {
 if(fgets(line,M,inputfile) == NULL)
  break;
// printf("%d:%s",linenum,line);
 }

//printf("%d\n",linenum);

fclose(inputfile);
return(linenum);
}



int getmaxint(char *inputfilename)
 {
 int i,datanum,max,data;
 FILE *inputfile;

 datanum = getintnum(inputfilename);
 inputfile = openfile(inputfilename,'r');
 fscanf(inputfile,"%d",&data);
// printf("%d\n",datanum);
 max = data;
 for(i=0;i<datanum-1;i++)
  {
  fscanf(inputfile,"%d",&data);
  if(data > max)
   max = data;
  }
 fclose(inputfile);

 return(max);
 }




