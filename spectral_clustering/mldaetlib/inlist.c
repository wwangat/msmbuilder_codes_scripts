int doublevalueequal(double value1,double value2,double diff)
 {
 int equal;
 if((value1 > value2-diff) && (value1 < value2+diff))
  equal = 1;
 else
  equal = 0;

 return(equal);
 }







int doublevectorequal(double *vector1,double *vector2,int dim, double diff)
 {
 int equal,i;

 for(i=0,equal=1;(i<dim)&&equal;i++)
  if(doublevalueequal(vector1[i],vector2[i],diff) == 0)
   equal = 0;

 return(equal);
 }






int indoublevectorlist(double **vectorlist,int listlen,double *vector,int dim,double diff)
 {
 int i,found;

 found = -1;
 for(i=0;i<listlen;i++)
  {
  if(doublevectorequal(vectorlist[i],vector,dim,diff) == 1)
   {
   found = i;
   break;
   }
  }
 return(found);
 }






int inintlist(int *array,int len,int value)                //if the value is in the array, return the position of the value, else return -1
 {
 int found,i;

 found = -1;

 if(len == 0)
  return(found);
 for(i=0;i<len;i++)
  {
  if(array[i] == value)
   {
   found = i;
   break;
   }
  }

 return(found);
 }







int incharlist(char *array,char value)                //if the value is in the array, return the position of the value, else return -1
 {
 int found,i;

 found = -1;

// if(len == 0)
//  return(found);
 for(i=0;;i++)
  {
  if(array[i] == '\0')
   break;

  if(array[i] == value)
   found = i;
  }

 return(found);
 }



int instringlist(char **strarray,int len, char *string)
 {
 int found,i;
 found = -1;


 if(len == 0)
  return(found);

 for(i=0;i<len;i++)
  {
  if(strcmp(strarray[i],string) == 0)
   {
   found = i;
   break;
   }
  }

 return(found);
 }

