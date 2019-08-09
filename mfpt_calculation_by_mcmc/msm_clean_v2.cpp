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

