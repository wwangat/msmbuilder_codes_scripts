#include <math.h>



void vectorcrosstime(double *array1, double *array2, double *array3)
{
int i;
array3[0] = array1[1]*array2[2]-array1[2]*array2[1];
array3[1] = array1[2]*array2[0]-array1[0]*array2[2];
array3[2] = array1[0]*array2[1]-array1[1]*array2[0];
}




double vectordottime(double *array1,double *array2,int dim)
{
double result;
int i;
result = 0;
for(i=0;i<dim;i++)
 result += array1[i]*array2[i];
return(result);
}



void vectornormalize(double *array,int dim)
 {
 int i;
 double mode;
 mode = 0;

 for(i=0;i<dim;i++)
  mode += array[i]*array[i];
 mode = sqrt(mode);
// printf("%6f   %6f   %6f    %f\n",array[0],array[1],array[2],mode);//test

 for(i=0;i<dim;i++)
  array[i] = array[i]/mode;
 }
