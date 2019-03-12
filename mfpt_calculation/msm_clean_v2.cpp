/*
 * ============================================================================
 *       Filename:  msm_clean.cpp
 *    Description:  Count transition matrix, implied timescale, recover of long
 *    time kinetics, mfpt by counting
 *		functions:	void transmatrix(double **ProbMatrix, int dim, int lagtime, int nline, int *MicroAssign, int *traj_len, int traj_num)
 *			void ImpliedTimescale(int *assignment, int nline, int start, int interval, int end, int dim, int *traj_len, int traj_num)
 *		double Overlap_eigR(int *micro, int *macro, int *mapping, int lagtime, int nline, int nMicro, int nMacro, int *traj_len, int traj_num)
 *		double real_TPM(char *filename, int nMacro, int lagtime, int terminate, int nline, int *macro, int *traj_len, int traj_num)
 *		double GOE_TPM(char *trad_file, char *GOE_file, int nMacro, int nMicro, int lagtime, int terminate, int nline, int *micro, int *macro, int *traj_len, int traj_num)
 *      double mfpt_direct_count(int *macro, int *traj_len, int traj_num, int nMacro, double **mfpt_matrix)
 *        Created:  2015-11-09 01:45
 *         Author:  Wei WANG        (wwangat@gmail.com)
 * ============================================================================
 Reference: Da, L., Sheong, F.K., Silva, D., Huang, X.
 "Application of Markov State Models to Simulate Long Timescale Dynamics of Biological Macromolecules",
  Book Chapter, Protein Conformational Dynamics, Advances in Experimental Medicine and Biology. 805, 29-66, Springer, (2014)
 */

#include <iostream>
#include <math.h>
//#include "dgeev.cpp"
//#include "read.cpp"
//#include "allocate.cpp"
using namespace std;

/*Part I: generate transition probability matrix with column normalized*/
//transCount: default: 7 args, sliding window, sliding step=lagtime
void transCount(double **CountMatrix, int dim, int lagtime, int nline, int *MicroAssign, int *traj_len, int traj_num){
	int i, j;
	for (i=0; i<dim; i++){
		for (j=0; j<dim;j++){
			CountMatrix[i][j] = 0.0;
		}
	}
	int curr, prev;
	if (traj_num==1){
		for(i=0;i<nline-lagtime;i++){
			curr = MicroAssign[i+lagtime];
	 		prev = MicroAssign[i];
			if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
		}
	}
	else{
		for(i=0; i<traj_len[0]-lagtime; i++){
			curr = MicroAssign[i+lagtime];
			prev = MicroAssign[i];
			if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
		}
		int temp = traj_len[0];
		for(i=0; i<traj_num-1;i++){
			for (j=temp; j<temp+traj_len[i+1]-lagtime; j++){
				curr = MicroAssign[j+lagtime];
				prev = MicroAssign[j];
				if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
			}
			temp += traj_len[i+1];
		}
	}
	//symmetrize of transition count matrix
	for (i = 0; i<dim; i++){
		for (j = i+1; j<dim; j++){
			CountMatrix[i][j] = (CountMatrix[i][j]+CountMatrix[j][i])/2.0;
			CountMatrix[j][i] = CountMatrix[i][j];
		}
	}
}
//new function: 8 args, jumping window with jumping step=jump_step
void transCount(double **CountMatrix, int dim, int lagtime, int nline, int *MicroAssign, int *traj_len, int traj_num, int jump_step){
	int i, j;
	for (i=0; i<dim; i++){
		for (j=0; j<dim;j++){
			CountMatrix[i][j] = 0.0;
		}
	}
	int curr, prev;
	if (traj_num==1){
		for(i=0;i<nline-lagtime;i=i+jump_step){
			curr = MicroAssign[i+lagtime];
	 		prev = MicroAssign[i];
			if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
		}
	}
	else{
        //for the 1st trajectory
		for(i=0; i<traj_len[0]-lagtime; i=i+jump_step){
			curr = MicroAssign[i+lagtime];
			prev = MicroAssign[i];
			if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
		}
		int temp = traj_len[0];
        //for other trajectories
		for(i=0; i<traj_num-1;i=i+1){
			for (j=temp; j<temp+traj_len[i+1]-lagtime; j=j+jump_step){
				curr = MicroAssign[j+lagtime];
				prev = MicroAssign[j];
				if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
			}
			temp += traj_len[i+1];
		}
	}
	//symmetrize of transition count matrix
	for (i = 0; i<dim; i++){
		for (j = i+1; j<dim; j++){
			CountMatrix[i][j] = (CountMatrix[i][j]+CountMatrix[j][i])/2.0;
			CountMatrix[j][i] = CountMatrix[i][j];
		}
	}
}

//the transition count matrix is not symmetrized
void transCount(double **CountMatrix, int dim, int lagtime, int nline, int *MicroAssign, int *traj_len, int traj_num, int jump_step, int sym){
	//sym is either 0 or 1, 1:do symmetrize, 0:not symmetrize
	int i, j;
	for (i=0; i<dim; i++){
		for (j=0; j<dim;j++){
			CountMatrix[i][j] = 0.0;
		}
	}
	int curr, prev;
	if (traj_num==1){
		for(i=0;i<nline-lagtime;i=i+jump_step){
			curr = MicroAssign[i+lagtime];
	 		prev = MicroAssign[i];
			if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
		}
	}
	else{
        //for the 1st trajectory
		for(i=0; i<traj_len[0]-lagtime; i=i+jump_step){
			curr = MicroAssign[i+lagtime];
			prev = MicroAssign[i];
			if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
		}
		int temp = traj_len[0];
        //for other trajectories
		for(i=0; i<traj_num-1;i=i+1){
			for (j=temp; j<temp+traj_len[i+1]-lagtime; j=j+jump_step){
				curr = MicroAssign[j+lagtime];
				prev = MicroAssign[j];
				if(prev>=0 && curr>=0)	CountMatrix[curr][prev]++;
			}
			temp += traj_len[i+1];
		}
	}
	if (sym==1){	
		//which means we need to symmetrize the matrix
		//symmetrize of transition count matrix
		for (i = 0; i<dim; i++){
			for (j = i+1; j<dim; j++){
				CountMatrix[i][j] = (CountMatrix[i][j]+CountMatrix[j][i])/2.0;
				CountMatrix[j][i] = CountMatrix[i][j];
			}
		}
	}
}


//currently transition probability matrix only for sliding window
void transmatrix(double **ProbMatrix, int dim, int lagtime, int nline, int *MicroAssign, int *traj_len, int traj_num){
	int i, j;
	transCount(ProbMatrix, dim, lagtime, nline, MicroAssign, traj_len, traj_num);
	//column normalized
	double temp;
	for (i = 0; i<dim; i++){
		temp = 0.0;
		for (j = 0; j<dim; j++){
			temp += ProbMatrix[j][i];
		}
		if(temp==0){
			printf("\33[31mThis count matrix has zero rowsum term(state %d), need to remove these states first\33[0m\n", i);
		}
		//for each column
		for (j = 0; j<dim; j++){
			ProbMatrix[j][i] /= temp;
		}
	}
}
//return both transCount and transProb
void transmatrix(double **CountMatrix, double **ProbMatrix, int dim, int lagtime, int nline, int *MicroAssign, int *traj_len, int traj_num){
	int i, j;
	transCount(CountMatrix, dim, lagtime, nline, MicroAssign, traj_len, traj_num);
	//column normalized
	double temp;
	for (i = 0; i<dim; i++){
		temp = 0.0;
		for (j = 0; j<dim; j++){
			temp += CountMatrix[j][i];
		}
		if(temp==0){
			printf("\33[31mThis count matrix has zero rowsum term(state %d), need to remove these states first\33[0m\n", i);
		}
		//for each column
		for (j = 0; j<dim; j++){
			ProbMatrix[j][i] = CountMatrix[j][i]/temp;
		}
	}
}


/*calculate implied timescale in different timescale*/
/*
void ImpliedTimescale(int *assignment, int nline, int start, int interval, int end, int dim, int *traj_len, int traj_num){
	//Just need to record eigenvalues, no need other parts
	double **ProbMatrix, *eigenVal_r, *eigenVal_i;
	ProbMatrix = allomatrix_double(dim, dim); //column normalized
	int numRec;
	if (dim>=10){
		numRec = 10;
	}
	else{
		numRec = dim;
	}
	eigenVal_r = alloarray_double(dim);
	eigenVal_i = alloarray_double(dim);
	for (int k = start; k<=end; k=k+interval){
		transmatrix(ProbMatrix, dim, k, nline, assignment, traj_len, traj_num);
	 	eig(ProbMatrix, eigenVal_r, eigenVal_i, dim);
	 	printf("major implied timescale for lagtime %d are as follows:\n", k);
	 	for (int i = 1; i<numRec;i++){
	 		printf("%lf\t", -k/log(eigenVal_r[i]));//from the second
	 	}
	 	printf("\n");
	}
	free(eigenVal_r);
	free(eigenVal_i);
	for (int i = 0; i<dim; i++){
		free(ProbMatrix[i]);
	}
	free(ProbMatrix);
}
void angle(double **matrix1, double **matrix2, double **angleMatrix, int nMicro, int nMacro, double *micropop){
    double temp1, temp2, temp;
    double *vector1 = alloarray_double(nMicro), *vector2 = alloarray_double(nMicro);
    for (int i=0; i<nMacro; i++){
        for (int k=0; k<nMicro; k++){
            vector1[k] = matrix1[i][k];
        }
        temp1 = sqrt(dotProduct(vector1, vector1, nMicro, micropop));
        for (int j=0; j<nMacro; j++){
            for (int k=0; k<nMicro; k++){
                vector2[k] = matrix2[j][k];
            }
            temp2 = sqrt(dotProduct(vector2, vector2, nMicro, micropop));
            temp = dotProduct(vector1, vector2, nMicro, micropop);
            angleMatrix[i][j] = fabs(temp)/(temp1*temp2);
        }
    }
    free(vector1);
    free(vector2);
}

double Overlap_eigR(int *micro, int *macro, int *mapping, int lagtime, int nline, int nMicro, int nMacro, int *traj_len, int traj_num){
	////////////microstate eigs calculation, record both eigenvalues and right eigenvector overlap//////////////////
	printf("At lagtime %d, leading eigenvalues and eigenvector overlap:\n", lagtime);
	double **micro_Matrix, *micro_eigenVal_r, **micro_eigenVec, *micro_eigenVal_i;
	micro_Matrix = allomatrix_double(nMicro, nMicro); //column normalized
	micro_eigenVal_r = alloarray_double(nMicro);
	micro_eigenVal_i = alloarray_double(nMicro);
	micro_eigenVec = allomatrix_double(nMicro, nMicro);
	transmatrix(micro_Matrix, nMicro, lagtime, nline, micro, traj_len, traj_num);
	eig(micro_Matrix, micro_eigenVal_r, micro_eigenVal_i, micro_eigenVec, nMicro);///sorted ones, each line an eigenvector
	cout << "Major MicroState EigenValue:\t";
	for (int i = 0; i<nMacro;i++){
		printf("%lf\t", micro_eigenVal_r[i]);
	}
	printf("\n");
	double **macro_Matrix, *macro_eigenVal_r, **macro_eigenVec, *macro_eigenVal_i;
	macro_Matrix = allomatrix_double(nMacro, nMacro); //column normalized
	macro_eigenVal_r = alloarray_double(nMacro);
	macro_eigenVal_i = alloarray_double(nMacro);
	macro_eigenVec = allomatrix_double(nMacro, nMacro);
	transmatrix(macro_Matrix, nMacro, lagtime, nline, macro, traj_len, traj_num);
	eig(macro_Matrix, macro_eigenVal_r, macro_eigenVal_i, macro_eigenVec, nMacro);///sorted ones, each line an eigenvector
	cout << "MacroState EigenValue:\t";
	for (int i = 0; i<nMacro;i++){
		printf("%lf\t", macro_eigenVal_r[i]);
	}
	printf("\n");
	///////////////////////////Eigenvector overlap/////////////////////
	double *micropop = alloarray_double(nMicro), *invmacropop = alloarray_double(nMacro);
	for (int i=0; i<nMacro; i++){
		invmacropop[i] = 0.0;
	}
	double temp = 0.0;
	for (int i=0; i<nMicro; i++){
		temp += micro_eigenVec[0][i];
	}
	for (int i=0; i<nMicro; i++){
		micropop[i] = micro_eigenVec[0][i]/temp;
		invmacropop[mapping[i]]+=micropop[i];
	}
	for (int i=0; i<nMacro; i++){
		invmacropop[i] = 1.0/invmacropop[i];
	}
	/////////////getting Dn, A, inv(DN), now we can calculating projected back eigenvectors///////////////////
	double **right_PK = allomatrix_double(nMacro, nMicro);///projected back one
	int temp1; //membership function
	for (int i=0; i<nMicro; i++){
		temp1 = mapping[i];//mapping should start from 0
		temp = micropop[i]*invmacropop[temp1];
		for (int j = 0; j<nMacro; j++){
			right_PK[j][i] = temp*macro_eigenVec[j][temp1];
		}
	}
	/////////////Having got all backward projected eigenvectors, now calculating overlap dot products////////
	double **angleMatrix=allomatrix_double(nMacro, nMacro);
	angle(right_PK, micro_eigenVec, angleMatrix, nMicro, nMacro, micropop);
	cout << "Eigenvector overlap matrix is:\n";
	for (int i = 0; i<nMacro;i++){
		for (int j = 0; j<nMacro; j++){
			printf("%lf\t", angleMatrix[i][j]);
		}
 	printf("\n");
	}
	cout<<"***********************************************************"<<endl;
	///////////////begin to calculate time dilation factor under this lagtime////////////////////
	cout <<"Major time dilation factor except the first equilibrium mode"<<endl;
	double time_dilation, temp_sum=0.0;
	for (int i = 1; i<nMacro;i++){
	temp_sum = temp_sum+log(micro_eigenVal_r[i])/log(macro_eigenVal_r[i]);
			printf("%lf\t", log(micro_eigenVal_r[i])/log(macro_eigenVal_r[i]));
		}
	cout<<endl;
	time_dilation = temp_sum/(nMacro-1);
	cout <<"Time scaling factor under this lagtime (average of time dilation factor) is " <<time_dilation <<endl;
	////////////////////////////
	for (int i = 0; i<nMicro; i++){
		free(micro_eigenVec[i]);
		free(micro_Matrix[i]);
	}
	for (int i = 0; i<nMacro; i++){
		free(macro_eigenVec[i]);
		free(macro_Matrix[i]);
		free(angleMatrix[i]);
		free(right_PK[i]);
	}
	free(micro_eigenVal_r);
	free(micro_eigenVal_i);
	free(right_PK);
	free(macro_eigenVal_r);
	free(macro_eigenVal_i);
	free(micro_eigenVec);
	free(micro_Matrix);
	free(macro_eigenVec);
	free(macro_Matrix);
	free(angleMatrix);
	return time_dilation;
}


///If we know eigenvectors of microstate TPM and mapping scheme in advance, then using the below function
//still writing, WW on Feb 25, 2016
double Overlap_eigR(double **micro_eigenVec,double *micropop, double **macro_Matrix, int *mapping, int nMicro, int nMacro){
	////////////microstate eigs calculation, record both eigenvalues and right eigenvector overlap//////////////////
	double *macro_eigenVal_r, **macro_eigenVec, *macro_eigenVal_i;
	macro_eigenVal_r = alloarray_double(nMacro);
	macro_eigenVal_i = alloarray_double(nMacro);
	macro_eigenVec = allomatrix_double(nMacro, nMacro);
	//macro_Matrix is the converged tCount in main program
	double temp;
	for (int i = 0; i<nMacro; i++){
		temp = 0.0;
		for (int j = 0; j<nMacro; j++){
			temp += macro_Matrix[j][i];
		}
		//for each column
		for (int j = 0; j<nMacro; j++){
			macro_Matrix[j][i] /= temp;
		}
	}
	eig(macro_Matrix, macro_eigenVal_r, macro_eigenVal_i, macro_eigenVec, nMacro);///sorted ones, each line an eigenvector
	///////////////////////////Eigenvector overlap/////////////////////
	double *invmacropop = alloarray_double(nMacro);
	for (int i=0; i<nMacro; i++){
		invmacropop[i] = 0.0;
	}
	for (int i=0; i<nMicro; i++){
		invmacropop[mapping[i]]+=micropop[i];
	}
	for (int i=0; i<nMacro; i++){
		invmacropop[i] = 1.0/invmacropop[i];
	}
	/////////////getting Dn, A, inv(DN), now we can calculating projected back eigenvectors///////////////////
	double **right_PK = allomatrix_double(nMacro, nMicro);///projected back one
	int temp1; //membership function
	for (int i=0; i<nMicro; i++){
		temp1 = mapping[i];//mapping should start from 0
		temp = micropop[i]*invmacropop[temp1];
		for (int j = 0; j<nMacro; j++){
			right_PK[j][i] = temp*macro_eigenVec[j][temp1];
		}
	}
	/////////////Having got all backward projected eigenvectors, now calculating overlap dot products////////
	double **angleMatrix=allomatrix_double(nMacro, nMacro);
	angle(right_PK, micro_eigenVec, angleMatrix, nMicro, nMacro, micropop);
	double evaluation = angleMatrix[0][0];
	for (int j=1;j<nMacro;j++){
		evaluation += angleMatrix[j][j];
	}
	for (int i = 0; i<nMacro; i++){
		free(macro_eigenVec[i]);
		free(angleMatrix[i]);
		free(right_PK[i]);
	}
	free(right_PK);
	free(macro_eigenVal_r);
	free(macro_eigenVal_i);
	free(macro_eigenVec);
	free(angleMatrix);
	free(invmacropop);
	return evaluation;
}


 double real_TPM(char *filename, int nMacro, int lagtime, int terminate, int nline, int *macro, int *traj_len, int traj_num){
	FILE *file=fopen(filename, "w");
 	double **ProbMatrix;
 	ProbMatrix = allomatrix_double(nMacro, nMacro); //column normalized
	//first element is identity matrix
	cout<<"Now recording real TPM under different lagtime in P[i][j] form"<<endl;
	fprintf(file, "%d\t", 0);
	for (int i=0; i<nMacro; i++){
		for(int j=0; j<nMacro; j++){
			if(j==i){
				fprintf(file, "%lf\t", 1.0);
			}
			else{
				fprintf(file, "%lf\t", 0.0);
			}
		}
	}
	fprintf(file, "\n");
 	for (int k = lagtime; k<=terminate; k=k+lagtime){
   		transmatrix(ProbMatrix, nMacro, k, nline, macro, traj_len, traj_num);
     	fprintf(file, "%d\t", k);
   		for (int i = 0; i<nMacro;i++){
       		for (int j=0;j<nMacro; j++){
           		fprintf(file, "%lf\t", ProbMatrix[i][j]);
       		}
   		}
   		fprintf(file, "\n");
 	}
 	fclose(file);
}


double GOE_TPM(char *trad_file, char *GOE_file, int nMacro, int nMicro, int lagtime, int terminate, int nline, int *micro, int *macro, int *traj_len, int traj_num){
  	///Directly propagate  matrix in a certain lagtime
	int i,j, m;
	double **macro_Matrix;
	macro_Matrix = allomatrix_double(nMacro, nMacro); //column normalized
	transmatrix(macro_Matrix, nMacro, lagtime, nline, macro, traj_len, traj_num); //
	cout<<"recording directly propagated TPM in various lagtime in P[i][j] form:"<<endl;
	FILE *file1=fopen(trad_file, "w");
	fprintf(file1, "%d\t", 0);
	for (i=0; i<nMacro; i++){
		for(j=0; j<nMacro; j++){
			if(j==i){
				fprintf(file1, "%lf\t", 1.0);
			}
			else{
				fprintf(file1, "%lf\t", 0.0);
			}
		}
	}
	fprintf(file1, "\n");
 	double **ProbMatrix;
 	ProbMatrix = allomatrix_double(nMacro, nMacro); //column normalized
	transmatrix(ProbMatrix, nMacro, lagtime, nline, macro, traj_len, traj_num);
	double **temp_matrix=allomatrix_double(nMacro, nMacro);
	for (i=0; i<nMacro; i++){
		for (j=0; j<nMacro; j++){
			macro_Matrix[i][j]= ProbMatrix[i][j];
		}
	}
 	for (int k = lagtime; k<=terminate; k=k+lagtime){
     	fprintf(file1, "%d\t", k);
   		for (i = 0; i<nMacro;i++){
       		for (j=0;j<nMacro; j++){
           		fprintf(file1, "%lf\t", macro_Matrix[i][j]);
       		}
   		}
   		fprintf(file1, "\n");
		//begin to propagate
		for (i=0; i<nMacro; i++){
			for (j=0; j<nMacro; j++){
				temp_matrix[i][j] = 0.0;
 				for(m=0; m<nMacro; m++){
					temp_matrix[i][j]+=macro_Matrix[i][m]*ProbMatrix[m][j];
				}
			}
		}
		for(i=0; i<nMacro; i++){
			for(j=0; j<nMacro; j++){
				macro_Matrix[i][j] = temp_matrix[i][j];
			}
		}
	}
 	fclose(file1);

	FILE *file2 = fopen(GOE_file, "w");
 	cout<<"Now calculating GOE propagate matrix in P[i][j] form"<<endl;
 	double *macro_eigenVal_r, **macro_eigenVec, *macro_eigenVal_i; //only need right eigenvectors
 	double **micro_Matrix, *micro_eigenVal_r, *micro_eigenVal_i;
	macro_eigenVal_r = alloarray_double(nMacro);
	macro_eigenVal_i = alloarray_double(nMacro);
	macro_eigenVec = allomatrix_double(nMacro, nMacro);
	micro_Matrix = allomatrix_double(nMicro, nMicro); //column normalized
	micro_eigenVal_r = alloarray_double(nMicro);
	micro_eigenVal_i = alloarray_double(nMicro);
	eig(macro_Matrix, macro_eigenVal_r, macro_eigenVal_i, macro_eigenVec, nMacro);///sorted ones, each line an eigenvector
	transmatrix(micro_Matrix, nMicro, lagtime, nline, micro, traj_len, traj_num);
 	eig(micro_Matrix, micro_eigenVal_r, micro_eigenVal_i, nMicro);///sorted ones, each line an eigenvector
	double **macro_leftVec=allomatrix_double(nMacro, nMacro);
 	double *macropop = alloarray_double(nMacro);
 	double temp = 0.0;
 	for (int i=0; i<nMacro; i++){
 		temp += macro_eigenVec[0][i];
 	}
 	for (int i=0; i<nMacro; i++){
		macropop[i]=macro_eigenVec[0][i]/temp;
 	}
	double *array=alloarray_double(nMacro);
	//calculation of left eigenvector of macro TPM
	for (i=0; i<nMacro; i++){
		for(j=0;j<nMacro;j++){
			array[j] = macro_eigenVec[i][j];
		}
		temp = dotProduct(array, array, nMacro, macropop);//inner product in population norm space
		temp=sqrt(temp);
		for(j=0;j<nMacro;j++){
			macro_eigenVec[i][j]=macro_eigenVec[i][j]/temp;
		}
	}
	fprintf(file2, "%d\t", 0);
	for (i=0; i<nMacro; i++){
		for(j=0; j<nMacro; j++){
			if(j==i){
				fprintf(file2, "%lf\t", 1.0);
			}
			else{
				fprintf(file2, "%lf\t", 0.0);
			}
		}
	}
	fprintf(file2, "\n");
 	for (int k = lagtime; k<=terminate; k=k+lagtime){
     	fprintf(file1, "%d\t", k);
   		for (i = 0; i<nMacro;i++){
       		for (j=0;j<nMacro; j++){
				macro_Matrix[i][j]=0.0;
				for(int m=0;m<nMacro;m++){
					macro_Matrix[i][j]+=macro_eigenVec[m][i]*pow(micro_eigenVal_r[m],k/lagtime)*macro_eigenVec[m][j]/macropop[j];
				}
           		fprintf(file2, "%lf\t", macro_Matrix[i][j]);
       		}
   		}
   		fprintf(file2, "\n");
     	if (k%(100*lagtime)==0) fprintf(stderr, "write data : %d timestep finished\r", k%lagtime+1);
	}
	fclose(file2);
	for (int i = 0; i<nMicro; i++){
		free(micro_Matrix[i]);
	}
	for (int i = 0; i<nMacro; i++){
		free(macro_eigenVec[i]);
		free(macro_Matrix[i]);
	}
	free(micro_eigenVal_r);
	free(micro_eigenVal_i);
	free(macropop);
	free(array);
	free(macro_eigenVal_r);
	free(macro_eigenVal_i);
	free(micro_Matrix);
	free(macro_eigenVec);
	free(macro_Matrix);
 }
*/ 
 /////////////calculate mean first passage time by directly counting
 double mfpt_direct_count(int *macro, int *traj_len, int traj_num, int nMacro, double **mfpt_matrix){
    double **mfpt_count = allomatrix_double(nMacro, nMacro);
    double **mfpt_time = allomatrix_double(nMacro, nMacro);
    int *mflag = alloarray_int(nMacro);
    int i, j, k, old_index, new_index;
    for(i=0;i<nMacro;i++){
        for(j=0;j<nMacro;j++){
            mfpt_count[i][j] = 0.0;//how long it takes to reach another macrostate for the first time
            mfpt_time[i][j] = 0.0;//how many times of this kind of transition
        }
    }
    /////now splitting the macrostate chains into pieces according to trajectory lengths
    if(traj_num != 1){
        int temp = 0;
        //allocate intermediate dynamic array to store splitted trajectories
        for (i = 0;i<traj_num;i++){
            int *chain = alloarray_int(traj_len[i]);
            for (j = temp; j<traj_len[i]+temp; j++){
                chain[j-temp] = macro[j];
            }
            temp += traj_len[i];
            //starting counting for this trajectory
            for(j = 0; j<traj_len[i]-1;j++){
                for(k=0;k<nMacro;k++){
                    mflag[k] = 0;//set flag
                }
                old_index = chain[j];
                for(k = 1; k<traj_len[i]-j; k++){
                    if(sum_int(mflag, nMacro) == nMacro){
                        break;
                    }//mean first passage time for each macrostate is considered
                    else{
                        new_index = chain[j+k];
                        if(mflag[new_index] == 0){
                            mfpt_count[new_index][old_index] += 1;
                            mfpt_time[new_index][old_index] += k;
                            mflag[new_index] = 1;
                        }
                    }

                }
            }
            ////free chain
            free(chain);
        }
    }
    else{
        //only have 1 trajectory, repeat previous steps
        //starting counting for this trajectory
        for(j = 0; j<traj_len[0]-1;j++){
            for(k=0;k<nMacro;k++){
                mflag[k] = 0;//set flag
            }
            old_index = macro[j];
            for(k = 0; k<traj_len[0]-j; k++){
                if(sum_int(mflag, nMacro) == nMacro){
                    break;
                }//mean first passage time for each macrostate is considered
                else{
                    new_index = macro[j+k];
                    if(mflag[new_index] == 0){
                        mfpt_count[new_index][old_index] += 1;
                        mfpt_time[new_index][old_index] += k;
                        mflag[new_index] = 1;
                    }
                }
            }
        }
    }
    for(j = 0; j<nMacro; j++){
        for(k=0; k<nMacro;k++){
            mfpt_matrix[j][k] = mfpt_time[j][k]/mfpt_count[j][k];
//            mfpt_std[j][k] = sqrt(mfpt_std[j][k]/mfpt_count[j][k]-mfpt_matrix[j][k]*mfpt_matrix[j][k]);
        }
    }
    //thus we have updated the mean first passage time
    //free mflag, mfpt_count, mfpt_time
    free(mflag);
    for(j = 0;j<nMacro;j++){
    	free(mfpt_count[j]);
    	free(mfpt_time[j]);
    }
 }
 double mfpt_direct_count_multi(int *macro, int *traj_len, int traj_num, int nMacro, double **mfpt_count, double **mfpt_time){
    int *mflag = alloarray_int(nMacro);
    int i, j, k, old_index, new_index;
    /////now splitting the macrostate chains into pieces according to trajectory lengths
    if (traj_num == 1){
        //only have 1 trajectory, repeat previous steps
        //starting counting for this trajectory
        for (k=0; k<nMacro; k++){
            mflag[k] = 0;
        }
        old_index = macro[0];
        for(j = 0; j<traj_len[0];j++){
            if(sum_int(mflag, nMacro) == nMacro){
                break;
            }//mean first passage time for each macrostate is considered
            else{
                new_index = macro[j];
                if(mflag[new_index] == 0){
                    mfpt_count[new_index][old_index] += 1;
                    mfpt_time[new_index][old_index] += j;
                    mflag[new_index] = 1;
                }
            }
        }
    }
    //thus we have updated the mean first passage time
    //free mflag, mfpt_count, mfpt_time
    free(mflag);
}
