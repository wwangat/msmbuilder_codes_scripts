#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void sym(double **matrix,int dim, int r1, int r2)
 {
 double ct,t,rss,check,**tempmatrix,sin,cos,temp1,temp2;
 int i,sign,rr;

 int j;

 if((matrix[r1][r1] + matrix[r2][r2]) == 0)
  {
  tempmatrix = doublematrix(2,2);
  ct = 1/(matrix[r2][r1] - matrix[r1][r2])*(matrix[r1][r1]+matrix[r2][r2]);
  rss = sqrt(1+ct*ct);
  sign = 1;

  for(rr=0;rr<dim;rr++)
   if((rr != r1) && (rr !=r2))
    break;
  
//  check = ((1-rss)*(matrix[r1][r1]+matrix[r2][r2])*(matrix[r1][r1]+matrix[r2][r2]) + (matrix[r2][r1]-matrix[r1][r2])*(matrix[r2][r1]-matrix[r1][r2]))/(matrix[r1][r1]+matrix[r2][r2]);
  check = matrix[r1][r1]+matrix[r2][r2];
  if(check < 0)
   sign = -1;
  sin = sign/rss;
  cos = ct/rss;
  

  tempmatrix[0][0] = matrix[r1][r1]*cos-matrix[r1][r2]*sin;
  tempmatrix[0][1] = matrix[r1][r1]*sin+matrix[r1][r2]*cos;
  tempmatrix[1][0] = matrix[r2][r1]*cos-matrix[r2][r2]*sin;
  tempmatrix[1][1] = matrix[r2][r1]*sin+matrix[r2][r2]*cos;
  temp1 = matrix[rr][r1]*cos-matrix[rr][r2]*sin;
  temp2 = matrix[rr][r1]*sin+matrix[rr][r2]*cos;


  
  matrix[r1][r1] = tempmatrix[0][0];
  matrix[r1][r2] = tempmatrix[0][1];
  matrix[r2][r1] = tempmatrix[1][0];
  matrix[r2][r2] = tempmatrix[1][1];
  matrix[rr][r1] = temp1;
  matrix[rr][r2] = temp2;

//debug
//printf("%d %d:\n",r1,r2);
//
//  for(i=0;i<dim;i++)
//   {
//   for(j=0;j<dim;j++)
//    printf("%lf   ",matrix[i][j]);
//   printf("\n");
//   }
//   printf("\n");

  for(i=0;i<2;i++)
   free(tempmatrix[i]);
  free(tempmatrix);

  }
 else
  {
  tempmatrix = doublematrix(2,2);
  t = (matrix[r2][r1] - matrix[r1][r2])/(matrix[r1][r1]+matrix[r2][r2]);
  rss = sqrt(1+t*t);
  sign = 1;

  for(rr=0;rr<dim;rr++)
   if((rr != r1) && (rr !=r2))
    break;
  
//  check = ((1-rss)*(matrix[r1][r1]+matrix[r2][r2])*(matrix[r1][r1]+matrix[r2][r2]) + (matrix[r2][r1]-matrix[r1][r2])*(matrix[r2][r1]-matrix[r1][r2]))/(matrix[r1][r1]+matrix[r2][r2]);
  check = matrix[r1][r1]+matrix[r2][r2];
  if(check < 0)
   sign = -1;
  t = sign*t;
  sin = t/rss;
  cos = 1/rss;
  

  tempmatrix[0][0] = matrix[r1][r1]*cos-matrix[r1][r2]*sin;
  tempmatrix[0][1] = matrix[r1][r1]*sin+matrix[r1][r2]*cos;
  tempmatrix[1][0] = matrix[r2][r1]*cos-matrix[r2][r2]*sin;
  tempmatrix[1][1] = matrix[r2][r1]*sin+matrix[r2][r2]*cos;
  temp1 = matrix[rr][r1]*cos-matrix[rr][r2]*sin;
  temp2 = matrix[rr][r1]*sin+matrix[rr][r2]*cos;


  
  matrix[r1][r1] = tempmatrix[0][0];
  matrix[r1][r2] = tempmatrix[0][1];
  matrix[r2][r1] = tempmatrix[1][0];
  matrix[r2][r2] = tempmatrix[1][1];
  matrix[rr][r1] = temp1;
  matrix[rr][r2] = temp2;

//debug
//printf("%d %d:\n",r1,r2);
//
//  for(i=0;i<dim;i++)
//   {
//   for(j=0;j<dim;j++)
//    printf("%lf   ",matrix[i][j]);
//   printf("\n");
//   }
//   printf("\n");

  for(i=0;i<2;i++)
   free(tempmatrix[i]);
  free(tempmatrix);
  }
 }


double doubleabs(double d)
 {
 if(d < 0)
  d = -d;
 return(d);
 }


double getunsymdiff(double **matrix,int dim)
 {
 double diff,tempdiff;
 int i,j;

 diff = 0;
 for(i=0;i<dim;i++)
  for(j=i+1;j<dim;j++)
   {
   tempdiff = doubleabs(matrix[i][j]-matrix[j][i]);   
   if(tempdiff > diff)
    diff = tempdiff;
//   printf("testdiff:%lf   %lf   %lf\n",tempdiff,matrix[i][j],matrix[j][i]);
   }
/* diff = abs(matrix[0][1]-matrix[1][0]);
 if(diff < abs(matrix[1][2]-matrix[2][1]))
   diff = abs(matrix[1][2]-matrix[2][1]) ;
 if(diff < abs(matrix[0][2]-matrix[2][0]))
   diff = abs(matrix[0][2]-matrix[2][0]) ;
*/
 return(diff);
 }


//double gettrace(double **matrix, int dim)
// {
// int i;
// double trace=0;
//
// for(i=0;i<dim;i++)
//  trace += matrix[i][i];
// return(trace);
// }



double getrmsd(double **coor1,double **coor2,int len)	//Here COMs of coor1 and coor2 are all set to 0
 {
 int i,j,k;
 double g1,g2,trace,maxtrace,**matrix,diff,rmsdsq;
 double *center1,*center2;

 matrix = doublematrix(3,3);
 center1 = doublearray(3);
 center2 = doublearray(3);

 g1=0;
 g2=0;
 for(i=0;i<3;i++)
  {
  center1[i] = 0;
  center2[i] = 0;
  }

 for(i=0;i<len;i++)
  {
  for(j=0;j<3;j++)
   {
   center1[j] += coor1[i][j]/len;
   center2[j] += coor2[i][j]/len;
   }
  }

 for(i=0;i<len;i++)
  {
  for(j=0;j<3;j++)
   {
   g1+=(coor1[i][j]-center1[j])*(coor1[i][j]-center1[j]);
   g2+=(coor2[i][j]-center2[j])*(coor2[i][j]-center2[j]);
   }
  }

maxtrace = 0;


  for(i=0;i<3;i++)
   for(j=0;j<3;j++)
    matrix[i][j] = 0;

  for(i=0;i<len;i++)
   for(j=0;j<3;j++)
    for(k=0;k<3;k++)
     matrix[j][k] += (coor1[i][j]-center1[j])*(coor2[i][k]-center2[k]);


 
//  diff = 100;
  diff = getunsymdiff(matrix,3);
//  printf("check1,%lf\n",diff);
  if(diff > 0.00001)
   {
   for(;diff > 0.00001;)
    {
//printf("test1\n");
    sym(matrix,3,0,1);
//printf("test2\n");
    sym(matrix,3,1,2);
//printf("test3\n");
    sym(matrix,3,2,0);
    diff = getunsymdiff(matrix,3);
    }
   }
//  printf("check2\n");
   trace = matrix[0][0] + matrix[1][1] + matrix[2][2];
   if(trace > maxtrace)
    maxtrace = trace;


 for(i=0;i<3;i++)
  free(matrix[i]);
 free(matrix);

 rmsdsq = (g1+g2-2*maxtrace)/len;

 if(rmsdsq < 0)
  rmsdsq = 0;
 return(sqrt(rmsdsq));

 }





