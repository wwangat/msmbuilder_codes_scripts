#include <stdlib.h>
#define INT 0
#define FLOAT 1
#define DOUBLE 2
#define STRING 3



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
