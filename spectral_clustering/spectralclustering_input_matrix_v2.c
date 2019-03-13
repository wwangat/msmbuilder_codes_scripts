#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mldaetlib/allocate.c"
#include "mldaetlib/read.c"
#include "mldaetlib/geteigen.c"
#include <time.h>
#define N 200

//the code is written by Meng Luming

//ww:spectral clustering program, input matrix should be in the format of: i   j   A[i, j], where A[i, j] is all the non-zero elements in the matrix, i can be less or bigger than j
//ww:Sep30,2017: input matrix: transition count matrix, can be symmetry or assymmetry

void rownormalize(double **leftmatrix,int dim)
 {
 int i,j;
 double rowsum;

 for(i=0;i<dim;i++)
  {
  rowsum = 0;
  for(j=0;j<dim;j++)
   rowsum += leftmatrix[i][j]*leftmatrix[i][j];
  rowsum = sqrt(rowsum);
  for(j=0;j<dim;j++)
   leftmatrix[i][j] /= rowsum;
  }
 }

double getdist(double *coor,int colnum,double *centercoor)
 {
 double dist;
 int i;

 dist = 0;
 for(i=0;i<colnum;i++)
  dist += (coor[i] - centercoor[i])*(coor[i] - centercoor[i]);
 dist = sqrt(dist);
 return(dist);
 }

double getdistsq(double *coor,int colnum,double *centercoor)
 {
 double dist;
 int i;

 dist = 0;
 for(i=0;i<colnum;i++)
  dist += (coor[i] - centercoor[i])*(coor[i] - centercoor[i]);
 return(dist);
 }

void getneardist(double **data,int rownum,int colnum,int *centeridlist,int centernum,double *centerdist)
 {
 int i,j,k;
 double dist,mindist;

 for(i=0;i<rownum;i++)
  {
  mindist = 10000000;
  for(j=0;j<centernum;j++)
   {
   dist = getdist(data[i],colnum,data[centeridlist[j]]);
   if(mindist > dist)
     mindist = dist;
   }
  centerdist[i] = mindist;
  }
 }


void getnewassign(double **data,int rownum,int colnum,double **centeraverage,int clusternum,int *newassign,int *assigncount)
 {
 int i,j,k,nearid;
 double dist,mindist;
 
 for(i=0;i<clusternum;i++)
  assigncount[i] = 0;

 for(i=0;i<rownum;i++)
  {
  mindist = 10000000;
  nearid = -1;
  for(j=0;j<clusternum;j++)
   {
   dist = getdist(data[i],colnum,centeraverage[j]);
   if(dist < mindist)
    {
    nearid = j;
    mindist = dist;
    }
   newassign[i] = nearid;
   assigncount[nearid]  ++;
   }
  }
 }

void getnewaverage(double **data,int rownum,int colnum,double **centeraverage,int clusternum,int *assign)
 {
 int i,j,*count;

 count = intarray(clusternum);
 for(i=0;i<clusternum;i++)
  {
  count[i] = 0;
  for(j=0;j<colnum;j++)
   centeraverage[i][j] = 0;
  }

 for(i=0;i<rownum;i++)
  {
  for(j=0;j<colnum;j++)
   {
   centeraverage[assign[i]][j] += data[i][j];
   }
   count[assign[i]] ++;
  }

 for(i=0;i<clusternum;i++)
  for(j=0;j<colnum;j++)
   centeraverage[i][j] /= count[i];

//for(i=0;i<clusternum;i++)
// printf("count:%d\n",count[i]);

 free(count);

 }






int compareassign(int *assign,int rownum,int *newassign)
 {
 int i,j,samebool;

 samebool = 1;
 for(i=0;(i<rownum) && (samebool == 1);i++)
  if(assign[i] != newassign[i])
   {
   samebool = 0;
//   printf("assign:%d newassign:%d\n",assign[i],newassign[i]);
   }
 return(samebool);
 }








double getdistsum(double **data,int rownum,int colnum,double **centeraverage,int *assign)
 {
 int i,j;
 double distsum,dist;

 distsum = 0;
 for(i=0;i<rownum;i++)
  distsum += getdistsq(data[i],colnum,centeraverage[assign[i]]);
 return(distsum);
 }



int checkzeroassign(int *count,int len)
 {
 int i,result;
 result = 0;
 for(i=0;i<len;i++)
  {
  if(count[i] == 0)
   {
   result = 1;
   break;
   }
  }
 return(result);
 }



void main()
{
char inputfilename[N],outputfilename[N];
int i,j,k,clusternum,*assign,*newassign,initseed,*centeridlist,maxid,convergebool,count,*assigncount;
double **data,**centeraverage,*centerdist,maxdist,distsum,mindistsum,test;
int statenum,bestseed;
double **matrix,*eigenvalue,*eigenvalue_i,**leftmatrix,*rowsum,tempsum;
FILE *inputfile,*outputfile;
double temp;
int row,column;
clock_t start, end;
double cpu_time_used;


start=clock();
printf("Input filename for your transition count matrix (dense matrix):\n");
scanf("%s",inputfilename);

printf("Input the dimension of your transition count matrix:\n");
scanf("%d",&statenum);

printf("Input the number of macrostates:\n");
scanf("%d",&clusternum);

printf("Input filename for output assignment:\n");
scanf("%s",outputfilename);

matrix = doublematrix(statenum,statenum);
eigenvalue = doublearray(statenum);
eigenvalue_i = doublearray(statenum);
leftmatrix = doublematrix(statenum,statenum);
rowsum = doublearray(statenum);
data = doublematrix(statenum,clusternum);
centeraverage = doublematrix(clusternum,clusternum);
centerdist = doublearray(statenum);
assign = intarray(statenum);
newassign = intarray(statenum);
centeridlist = intarray(clusternum);
assigncount = intarray(clusternum);

//read file
//annotate by Wei
  inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 {
 rowsum[i] = 0;
 for(j=0;j<statenum;j++)
  {
  fscanf(inputfile,"%lf",&matrix[i][j]);
//  rowsum[i] += matrix[i][j];
  if(matrix[i][j] < 0)
   {
   printf("Error for negative value in original matrix: %d %d\n",i,j);
   exit(0);
   }
  }
 }
fclose(inputfile);

/*
//add by Wei
inputfile = openfile(inputfilename,'r');
for (i=0;i<statenum;i++)
{
    for (j=0;j<statenum;j++)
    {
        matrix[i][j] = 0;
    }
} 
for (i=0;i<length;i++)
{
    for (j=0;j<3;j++)
    {  
        fscanf(inputfile, "%d %d %lf", &row, &column, &temp);
        matrix[row-1][column-1] = temp;
    }
}
fclose(inputfile);
*/
printf("finish reading input matrix, matrix[1][1] and matrix[%d][%d] is %lf%lf respectively\n", statenum, statenum, matrix[0][0], matrix[statenum-1][statenum-1]);
//symmetry
for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  {
  if(i < j)
   matrix[i][j] = (matrix[i][j] + matrix[j][i])/2;
  else
   matrix[i][j] = matrix[j][i];
  rowsum[i] += matrix[i][j];
  }

//get Laplacian
for(i=0;i<statenum;i++)
 {
 matrix[i][i] = rowsum[i] - matrix[i][i];
 for(j=0;j<statenum;j++)
  if(i != j)
   matrix[i][j] = - matrix[i][j];
 rowsum[i] = sqrt(rowsum[i]);
 }

for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  matrix[i][j] /= rowsum[i]*rowsum[j];

//printf("matrix:\n");
//for(i=0;i<statenum;i++)
// {
// for(j=0;j<statenum;j++)
//  printf("%lf ",matrix[i][j]);
// printf("\n");
// }

geteigen(matrix,statenum,eigenvalue,eigenvalue_i,leftmatrix);
//geteigen_sort_dec(eigenvalue,statenum,leftmatrix);
printf("Eigenvalue:\n");
for(i=0;i<statenum;i++)
 printf("%lf %lf\n",eigenvalue[i],eigenvalue_i[i]);
//printf("Eigenmatrix:\n");
//for(i=0;i<statenum;i++)
// {
// for(j=0;j<statenum;j++)
//  printf("%lf ",leftmatrix[i][j]);
// printf("\n");
// }
geteigen_sort_inc(eigenvalue,statenum,leftmatrix);
rownormalize(leftmatrix,statenum);


//get data
for(i=0;i<statenum;i++)
 for(j=0;j<clusternum;j++)
  data[i][j] = leftmatrix[j][i];


//printf("data:\n");
//for(i=0;i<statenum;i++)
// {
// for(j=0;j<clusternum;j++)
//  printf("%lf  ",data[i][j]);
// printf("\n");
// }
//printf("check:\n");
//for(i=0;i<clusternum;i++)
// {
// test = 0;
// for(j=0;j<statenum;j++)
//  test += data[j][i]*data[j][i];
// printf("%e\n",test);
// }

//normalize data
for(i=0;i<statenum;i++)
 {
 tempsum = 0;
 for(j=0;j<clusternum;j++)
  tempsum += data[i][j]*data[i][j];
 tempsum = sqrt(tempsum);
 for(j=0;j<clusternum;j++)
  data[i][j] /= tempsum;
 }

//printf("Eigenvalue:\n");
//for(i=0;i<statenum;i++)
// printf("%e\n",eigenvalue[i]);
//
printf("data:\n");
for(i=0;i<statenum;i++)
 {
 for(j=0;j<clusternum;j++)
  printf("%lf  ",data[i][j]);
 printf("\n");
 }

//loop all the seed
bestseed = -1;
mindistsum = 10000000000;
for(initseed = 0;initseed < statenum; initseed ++)
 {
 //getinitcenter
 for(i=0;i<clusternum;i++)
  centeridlist[i] = -1;	//-1 for unassigned
 centeridlist[0] = initseed;
 
// printf("centerIDList:%d\n",centeridlist[0]);
 //get initial center id
 for(i=1;i<clusternum;i++)
  {
  getneardist(data,statenum,clusternum,centeridlist,i,centerdist);
 // for(j=0;j<rownum;j++)
 //  printf("centerdist:%lf\n",centerdist[j]);
  maxdist = centerdist[0];
  maxid = 0;
  for(j=1;j<statenum;j++)
   if(maxdist < centerdist[j])
    {
    maxid = j;
    maxdist = centerdist[j];
    }
  centeridlist[i] = maxid;	//next id
//  printf("centerIDList:%d\n",centeridlist[i]);
  }
 
 getneardist(data,statenum,clusternum,centeridlist,4,centerdist);
// for(i=0;i<rownum;i++)
//  printf("centerdist:%lf\n",centerdist[i]);

 
 //get initial coordinate
 for(i=0;i<clusternum;i++)
  for(j=0;j<clusternum;j++)
   centeraverage[i][j] = data[centeridlist[i]][j];

 
 //do the loop for converage
 getnewassign(data,statenum,clusternum,centeraverage,clusternum,assign,assigncount);
 for(count=0,convergebool = 0; (convergebool != 1) && (count <= 1000);count ++)
  {
  getnewaverage(data,statenum,clusternum,centeraverage,clusternum,assign);
//  printf("Center coor:\n");
//  for(j=0;j<clusternum;j++)
//   {
//   for(k=0;k<colnum;k++)
//    printf("%lf  ",centeraverage[j][k]);
//   printf("\n");
//   }
  getnewassign(data,statenum,clusternum,centeraverage,clusternum,newassign,assigncount);
  convergebool = compareassign(assign,statenum,newassign);
 // printf("%d\n",convergebool);
  if(convergebool == 0)	//copy assign
   {
   for(i=0;i<statenum;i++)
    assign[i] = newassign[i];
   }
  }
 
 if(checkzeroassign(assigncount,clusternum) == 1)
  {
  printf("Empty cluster for initial seed %d\n",initseed);
  continue;
  }

 distsum = getdistsum(data,statenum,clusternum,centeraverage,assign);
 if(distsum < mindistsum)
  {
  mindistsum = distsum;
  bestseed = initseed;
  }
 printf("distsum:%lf\n",distsum);
 }


printf("mindistsum:%lf\n",mindistsum);

//get the best result
initseed = bestseed;
 //getinitcenter
 for(i=0;i<clusternum;i++)
  centeridlist[i] = -1;	//-1 for unassigned
 centeridlist[0] = initseed;
 
 //get initial center id
 for(i=1;i<clusternum;i++)
  {
  getneardist(data,statenum,clusternum,centeridlist,i,centerdist);
 // for(j=0;j<rownum;j++)
 //  printf("centerdist:%lf\n",centerdist[j]);
  maxdist = centerdist[0];
  maxid = 0;
  for(j=1;j<statenum;j++)
   if(maxdist < centerdist[j])
    {
    maxid = j;
    maxdist = centerdist[j];
    }
  centeridlist[i] = maxid;	//next id
//  printf("centerIDList:%d\n",centeridlist[i]);
  }
 
 getneardist(data,statenum,clusternum,centeridlist,4,centerdist);
// for(i=0;i<rownum;i++)
//  printf("centerdist:%lf\n",centerdist[i]);

 
 //get initial coordinate
 for(i=0;i<clusternum;i++)
  for(j=0;j<clusternum;j++)
   centeraverage[i][j] = data[centeridlist[i]][j];

 
 //do the loop for converage
 getnewassign(data,statenum,clusternum,centeraverage,clusternum,assign,assigncount);
 for(count=0,convergebool = 0; (convergebool != 1) && (count <= 1000);count ++)
  {
  getnewaverage(data,statenum,clusternum,centeraverage,clusternum,assign);
//  printf("Center coor:\n");
//  for(j=0;j<clusternum;j++)
//   {
//   for(k=0;k<colnum;k++)
//    printf("%lf  ",centeraverage[j][k]);
//   printf("\n");
//   }
  getnewassign(data,statenum,clusternum,centeraverage,clusternum,newassign,assigncount);
  convergebool = compareassign(assign,statenum,newassign);
 // printf("%d\n",convergebool);
  if(convergebool == 0)	//copy assign
   {
   for(i=0;i<statenum;i++)
    assign[i] = newassign[i];
   }
  }
 
 if(checkzeroassign(assigncount,clusternum) == 1)
  {
  printf("Error:empty cluster exists %d\n",initseed);
  exit(0);
  }


 outputfile = openfile(outputfilename,'w');
 for(i=0;i<statenum;i++)
  fprintf(outputfile,"%d\n",assign[i]);
 fclose(outputfile);
 



for(i=0;i<statenum;i++)
 {
 free(matrix[i]);
 free(leftmatrix[i]);
 free(data[i]);
 }

for(i=0;i<clusternum;i++)
 {
 free(centeraverage[i]);
 }

free(matrix);
free(eigenvalue);
free(eigenvalue_i);
free(leftmatrix);
free(rowsum);
free(data);
free(centeraverage);
free(assign);
free(newassign);
free(centeridlist);
free(assigncount);
end = clock();
cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;
printf("use a total of %lf seconds in spectral clustering\n", cpu_time_used);
}



