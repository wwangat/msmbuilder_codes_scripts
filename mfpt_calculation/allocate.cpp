/*
 * ============================================================================
 *       Filename:  allocate.cpp
 *		Description: dynamic allocate
 		Since could not overload functions according to parameter type, thus should define using different names
 *		Functions: int *alloarray_int(int dim), double *alloarray_double(int dim), int **allomatrix_int(int row, int col),double **allomatrix_double(int row, int col)
 *        Created:  2015-11-08 16:01
 *         Author:  Wei WANG        (wwangat@gmail.com)
 * ============================================================================
 */
 int *alloarray_int(int dim){
 	int *array;
	array = (int *)malloc(dim*sizeof(int));
	return array;
 }
 double *alloarray_double(int dim){
 	double *array;
	array = (double *)malloc(dim*sizeof(double));
	return array;
 }
 int **allomatrix_int(int row, int col){
 	int **matrix;
 	matrix = (int **)malloc(row*sizeof(int *));
 	int i;
 	for (i=0; i<row;i++){
 		matrix[i] = (int *)malloc(col*sizeof(int));
 	}
 	return matrix;
 }
 double **allomatrix_double(int row, int col){
 	double **matrix;
 	matrix = (double **)malloc(row*sizeof(double *));
 	int i;
 	for (i=0; i<row;i++){
 		matrix[i] = (double *)malloc(col*sizeof(double));
 	}
 	return matrix;
 }
/*  
double ***allo3dmatrix_double(int row, int col, int height){
    double ***matrix;
    matrix = (double ***)malloc(row*sizeof(double **));
    int i, j;
    for (i=0; i<row;i++){
        matrix[i] = (double **)malloc(col*sizeof(double *)){
            for (j = 0; j<col; j++){
                matrix[i][j] = (double *)malloc(height*sizeof(double));
            }
        }
    }
    return matrix;
}
int ***allo3dmatrix_int(int row, int col, int height){
    int ***matrix;
    matrix = (int ***)malloc(row*sizeof(int **));
    int i, j;
    for (i=0; i<row;i++){
        matrix[i] = (int **)malloc(col*sizeof(int *)){
            for (j = 0; j<col; j++){
                matrix[i][j] = (int *)malloc(height*sizeof(int));
            }
        }
    }
    return matrix;
}
*/
