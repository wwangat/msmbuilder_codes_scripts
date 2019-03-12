/*
 * ============================================================================
 *       Filename:  operator.cpp
 *	Description: min, max, sum, quicksort, dot
 Functions:int min(int *array, int dim), double min(double *array, int dim),
 			int max(int *array, int dim), double max(double *array, int dim)
 			double sum(int *array, int dim)
 			void quicksort(double arr[], int len, int indices[]), void quicksort(double arr[], int len)
			double dotProduct(double *vector1, double *vector2, int dim), double dotProduct(double *vector1, double *vector2, int dim, double *micropop)
			void matrix_array(double **matrix, double *array, int dim1, int dim2)
			void array_matrix(double *array, double **matrix, int dim1, int dim2)
			-------------------------------
 *        Created:  2015-11-08 16:01
 *         Author:  Wei WANG        (wwangat@gmail.com)
 * ============================================================================
 */
//////////////////////////min, max, sum of array///////////////////////////////////
int min(int *array, int dim){
    int temp;
    temp = 100000;
    for (int i=0; i<dim; i++){
        if (array[i]<temp){
            temp = array[i];
        }
    }
    return temp;
 }
 double min_double(double *array, int dim){
    double temp;
    temp = 100000;
    for (int i=0; i<dim; i++){
        if (array[i]<temp){
            temp = array[i];
        }
    }
    return temp;
 }

int max(int *array, int dim){
    int temp = -10000;
    for (int i = 0; i<dim;i++){
        if(array[i]>temp){
            temp = array[i];
        }
    }
    return temp;
}
double max_double(double *array, int dim){
    double temp = -10000;
    for (int i = 0; i<dim;i++){
        if(array[i]>temp){
            temp = array[i];
        }
    }
    return temp;
}
double sum(double *array, int dim){
    double temp= 0.0;
    int i;
    for (i = 0; i<dim; i++){
        temp = temp+array[i];
    }
    return temp;
}
double cumsum(double *array, int dim){
  int i;
  for (i=1; i<dim;i++){
    array[i] = array[i-1]+array[i];
  }
}

int sum_int(int *array, int dim){
    int temp= 0;
    int i;
    for (i = 0; i<dim; i++){
        temp = temp+array[i];
    }
    return temp;
}
////////////////////quick sort in 'descend order', revised based on web: http://stackoverflow.com/questions/9558804/quick-sort-in-c//////////////
void swap(double arr[], int i, int j) {
    double temp = arr[j];
    arr[j] = arr[i];
    arr[i] = temp;
}
void swap_int(int arr[], int i, int j) {
    int temp = arr[j];
    arr[j] = arr[i];
    arr[i] = temp;
}
void quicksort0(double arr[], int a, int b, int indices[]) {
    if (a >= b)
        return;
    double key = arr[a];
    int i = a + 1, j = b;
    while (i < j) {
        while (i < j && arr[j] <= key)
            --j;
        while (i < j && arr[i] >= key)
            ++i;
        if (i < j){
            swap_int(indices, i, j);
            swap(arr, i, j);
        }
    }
    if (arr[a] < arr[i]) {
        swap_int(indices, a, i);
        swap(arr, a, i);
        quicksort0(arr, a, i - 1, indices);
        quicksort0(arr, i + 1, b, indices);
    } else { // there is no left-hand-side
        quicksort0(arr, a + 1, b, indices);
    }
}
void quicksort0(double arr[], int a, int b) {
    if (a >= b)
        return;
    double key = arr[a];
    int i = a + 1, j = b;
    while (i < j) {
        while (i < j && arr[j] <= key)
            --j;
        while (i < j && arr[i] >= key)
            ++i;
        if (i < j){
            swap(arr, i, j);
        }
    }
    if (arr[a] < arr[i]) {
        swap(arr, a, i);
        quicksort0(arr, a, i - 1);
        quicksort0(arr, i + 1, b);
    } else { // there is no left-hand-side
        quicksort0(arr, a + 1, b);
    }
}
///return sorted array and sort indices
void quicksort(double arr[], int len, int indices[]) {
    quicksort0(arr, 0, len - 1, indices);
}
///only return sorted array
void quicksort(double arr[], int len) {
    quicksort0(arr, 0, len - 1);
}
///deal with integer type,sort in 'descend', for ascend, plz use inverse after sorting
void quicksort0_int(int arr[], int a, int b) {
    if (a >= b)
        return;
    int key = arr[a];
    int i = a + 1, j = b;
    while (i < j) {
        while (i < j && arr[j] <= key)
            --j;
        while (i < j && arr[i] >= key)
            ++i;
        if (i < j){
            swap_int(arr, i, j);
        }
    }
    if (arr[a] < arr[i]) {
        swap_int(arr, a, i);
        quicksort0_int(arr, a, i - 1);
        quicksort0_int(arr, i + 1, b);
    } else { // there is no left-hand-side
        quicksort0_int(arr, a + 1, b);
    }
}
void quicksort0_int(int arr[], int a, int b, int indices[]) {
    if (a >= b)
        return;
    double key = arr[a];
    int i = a + 1, j = b;
    while (i < j) {
        while (i < j && arr[j] <= key)
            --j;
        while (i < j && arr[i] >= key)
            ++i;
        if (i < j){
            swap_int(indices, i, j);
            swap_int(arr, i, j);
        }
    }
    if (arr[a] < arr[i]) {
        swap_int(indices, a, i);
        swap_int(arr, a, i);
        quicksort0_int(arr, a, i - 1, indices);
        quicksort0_int(arr, i + 1, b, indices);
    } else { // there is no left-hand-side
        quicksort0_int(arr, a + 1, b, indices);
    }
}
///only return sorted array
void quicksort_int(int arr[], int len) {
    quicksort0_int(arr, 0, len - 1);
}
void quicksort_int(int arr[], int len, int indices[]) {
    quicksort0_int(arr, 0, len - 1, indices);
}

///////inner product of vectors//////////
double dotProduct(double *vector1, double *vector2, int dim){
	double temp=0.0;
	for(int i=0; i<dim; i++){
	    temp += vector1[i]*vector2[i];
	}
	return temp;
}

void reverse(int *array, int dim){
  int temp;
  if(dim%2==0){
    for(int j=0;j<dim/2;j++){
      temp = array[j];
      array[j] = array[dim-j-1];
      array[dim-j-1] = temp;
    }
  }
  else{
    for(int j=0;j<dim/2;j++){
      temp = array[j];
      array[j] = array[dim-j-1];
      array[dim-j-1] = temp;
    }
    array[dim/2] = array[dim/2];
  }
}


double dotProduct(double *vector1, double *vector2, int dim, double *micropop){
//in microstate population norm space
    double temp=0.0;
    for(int i=0; i<dim; i++){
        temp += vector1[i]*vector2[i]/micropop[i];
    }
    return temp;
}

//////transform between matrix and array////////////
void matrix_array(double **matrix, double *array, int dim1, int dim2){
	for (int i=0; i<dim1; i++){
	    for (int j = 0; j<dim2; j++){
	        array[i*dim1+j] = matrix[i][j];
	    }
	}
}
void array_matrix(double *array, double **matrix, int dim1, int dim2){
	for(int i=0; i<dim1; i++){
	    for (int j=0; j<dim2; j++){
	        matrix[i][j] = array[i*dim1+j];
	    }
	}
}
