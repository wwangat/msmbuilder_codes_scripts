/*
 * ============================================================================
 *       Filename:  read.cpp
 *       Functions: count line of file, read file, micro2macro, mapping_scheme
 *        Created:  2015-11-08 16:01
 *         Author:  Wei WANG        (wwangat@gmail.com)
 * ============================================================================
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <iostream>
 using namespace std;
 //////////////////////count number of lines of input files//////////////////////////////////
 int countline(char *filename){
 	int nline = 0, i;
 	int ch;
 	FILE *file = fopen(filename, "r");
	if (!file) {
		perror("\033[31mfail to open this file\33[0m\n"); 	
		exit(1);
	}
    char input[40960];
    fseek(file, 0, SEEK_SET); nline = 0;
    while (fgets(input, sizeof(input), file)) nline ++;
	fclose(file);
	printf("This file has %d lines\n", nline);
	return nline;
 }


//////////////////////////////// //read files, inputing dimension of array or matrix/////////////////////////////////
 void readarray(char *filename, int *array){
     int i, dim;
     FILE *file;
     dim = countline(filename);
     file = fopen(filename, "r");
     if (!file){
         printf("\033[31mfail to open and read file\33[0m\n");
         exit(1);
     }
     for (i=0; i<dim;i++){
         fscanf(file, "%d", &array[i]);
//         cout<<array[i]<<endl;
     }
     fclose(file);
 }
 void readarray(char *filename, double *array){
     int i, dim;
     FILE *file;
     dim = countline(filename);
     file = fopen(filename, "r");
     if (!file){
         printf("\033[31mfail to open and read file\33[0m\n");
         exit(1);
     }
     for (i=0; i<dim;i++){
         fscanf(file, "%lf", &array[i]);
     }
     fclose(file);
 }
 void readarray(char *filename, float *array){
     int i, dim;
     FILE *file;
     dim = countline(filename);
     file = fopen(filename, "r");
     if (!file){
         printf("\033[31mfail to open and read file\33[0m\n");
         exit(1);
     }
     for (i=0; i<dim;i++){
         fscanf(file, "%f", &array[i]);
     }
     fclose(file);
 }
//below revised on July 4, 2016; To make it more general
 void readmatrix(char *filename, int **matrix, int col){
     int i, j, row;
     FILE *file;
     row = countline(filename);
     file = fopen(filename, "r");
     if (!file){
         printf("\033[31mfail to open and read file\33[0m\n");
         exit(1);
     }
     for (i=0; i<row;i++){
     	for (j = 0; j<col; j++){
	         fscanf(file, "%d", &matrix[i][j]);
     	}
     }
     fclose(file);
 }
 void readmatrix(char *filename, double **matrix, int col){
     int i, j, row;
     FILE *file;
     row = countline(filename);
     file = fopen(filename, "r");
     if (!file){
         printf("\033[31mfail to open and read file\33[0m\n");
         exit(1);
     }
     for (i=0; i<row;i++){
     	for (j = 0; j<col; j++){
	         fscanf(file, "%lf", &matrix[i][j]);
     	}
     }
     fclose(file);
 }


 ///////////////////special: reading coordinate, microassignment, macroassignment::for 1D potential output//////////////////////
 void read_microANDmacro_2files(char *filename1, char *filename2, int *micro, int *macro, int row){
     int i;
     FILE *file1, *file2;
     file1 = fopen(filename1, "r");
     file2 = fopen(filename2, "r");
     if (!file1 || !file2){
         printf("\033[31mfail to open and read these files\33[0m\n");
         exit(1);
     }
     for (i=0; i<row;i++){
       if (i%1000==0) fprintf(stderr, "main : %d lines read\r", i+1);
        fscanf(file1, "%d", &micro[i]);
        fscanf(file2, "%d", &macro[i]);
        macro[i] = macro[i];
     }
     fclose(file1);
     fclose(file2);
     //adjust index to make both micro and macro assignment start from 0
 }
void read_microANDmacro_1file(char *filename1, int *micro, int *macro, int row){
    int i;
    FILE *file1;
    double temp;
    file1 = fopen(filename1, "r");
    if (!file1){
        printf("\033[31mfail to open and read these files\33[0m\n");
        exit(1);
    }
    for (i=0; i<row;i++){
        if (i%1000==0) fprintf(stderr, "main : %d lines read\r", i+1);
        fscanf(file1, "%lf %d %d", &temp, &micro[i], &macro[i]);
    }
    fprintf(stderr, "\n");
    fclose(file1);
}
void read_microANDmapping_2files(char *filename1, char *filename2, int *micro, int *mapping, int row1, int row2){
    int i;
    FILE *file1, *file2;
    file1 = fopen(filename1, "r");
    file2 = fopen(filename2, "r");
    if (!file1 || !file2){
        printf("\033[31mfail to open and read these files\33[0m\n");
        exit(1);
    }
    for (i=0; i<row1;i++){
        if (i%1000==0) fprintf(stderr, "main : %d lines read\r", i+1);   
        fscanf(file1, "%d", &micro[i]);
    }
    for(i = 0; i<row2; i++){
        if (i%1000==0) fprintf(stderr, "main : %d lines read\r", i+1);
        fscanf(file2, "%d", &mapping[i]);
    }
    fclose(file1);
    fclose(file2);
}

//function: micro2macro: map microstate sequence to macrostate sequence, used in hard clustering
void micro2macro(int *micro, int *mapping, int *macro, int nline){
    for (int i=0; i<nline; i++){
        if(micro[i]>=0){
            macro[i] = mapping[micro[i]];
        }
        else{
            macro[i] = -1;
        }
    }
}
void mapping_scheme(int *micro, int *mapping, int *macro, int nline){
    for (int i=0; i<nline; i++){
        mapping[micro[i]] = macro[i];
    }
}



//function: micro2macro_membership: map microstate sequence to macrostate sequence with given probability (membership)
//if you have the hard assignment, then just use micro2macro
void micro2macro_membership(int *micro, double **membership, int *macro, int nMacro, int nline){
    double *array;
    for (int i=0;i<nline; i++){
        if (micro[i]>=0){
            array = membership[micro[i]];//get the head address for that row
            macro[i] = MonteCarlo(array, nMacro);
        }
        else{
            macro[i] = -1;
        }
    }
}
///////////////color/////////////////
/*
    \33[31m: red, \33[32m: green, \33[33m:yellow, \33[34m:blue, \33[35m:pink, \33[36m:cyan, \33[37m:grey
example: for ((i=30;i<=38;i++)); do printf "%d " $i; done | awk '{for(i=1;i<=NF;i++)printf("\33[%dm%d\33[0m ",$i,$i);printf("\n")}'
        for ((i=30;i<=38;i++)); do printf "%d " $i; done | awk '{for(i=1;i<=NF;i++)printf("\33[1;%dm%d\33[0m ",$i,$i);printf("\n")}'

*/
