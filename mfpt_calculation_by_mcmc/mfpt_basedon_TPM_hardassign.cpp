/*
 * ============================================================================
 *       Filename:  mfpt_basedon_TPM_hardassign.cpp
 *    Description:  input microstate TPM & micro-to-macro mapping, output MFPT. method: markov chain monte carlo
 *          Usage:  g++ mfpt_basedon_TPM_hardassign.cpp;./a.out
 *     Parameters:  nMacro, nMicro, MCtime, syn_len, timeunit, lagtime, ProbMatrix, mapping.
 *        Created:  2016-07-01 20:05
 *         Author:  Wei WANG        (wwangat@gmail.com)
 * ============================================================================
 */
 #include <iostream>
 #include <fstream>
 #include <string>
 #include <stdlib.h>
 #include <stdio.h>
 #include <time.h>
 #include <math.h>
 #include <sys/stat.h>
 #include <sys/types.h>
 #include "allocate.cpp"
 #include "basic.cpp"
 #include "operator.cpp"
 #include "read_v2.cpp"
 #include "msm_clean_v2.cpp"
 using namespace std;
 
int main()
 {
    ///////read multiple microstate trajectories
    int i, j, k, nline;
    long a = time(NULL);
    cout << "Seed: " << a << endl;
    srand(a);
    clock_t t=clock();
    double cpu_time_used;

//parameters you may change
    int nMacro=4;  //input; number of macrostates
    int nMicro = 1483; //input: number of microstates
    int MCtime = 5000000; //input:the length of each markov chain
    int syn_len = 1; //input: number of markov chains
    int *macro, *temp_micro;
    double timeunit=1e-4;//unit:micro-second, the time unit for each saving interval
    double lagtime=800; //input:the microstate markovian lag-time
//

    double **mfpt_matrix=allomatrix_double(nMacro, nMacro);
    double **ProbMatrix = allomatrix_double(nMicro, nMicro);
    
    readmatrix("microstate_TPM_col_normalized.txt", ProbMatrix, nMicro);//input: the input TPM should be column normalized, plz make sure that. can be either from transpose or mle
    int mapping[nMicro];
    readarray("pcca_plus_4_state_mapping.txt", mapping);//input: the lumping relationship. change the filename
    cout<<"step 2: artificially generate many short microstate trajectory based on the matrix calculated above, each one has a length of "<<MCtime<<endl;
    double **mfpt_count = allomatrix_double(nMacro, nMacro);
    double **mfpt_time = allomatrix_double(nMacro, nMacro);
    for (int i=0; i<nMacro; i++){
        for(int j=0; j<nMacro; j++){
            mfpt_count[i][j]=0.0;mfpt_time[i][j] = 0.0;
        }
    }
    macro = alloarray_int(MCtime);
    temp_micro = alloarray_int(MCtime);
    int traj_len[0];
    int *count_num=alloarray_int(nMacro);
    for (int kk=0;kk<nMacro;kk++)
    {
		count_num[kk]=0;
    }
    float total_count=0;
    for (int m=0; m<syn_len; m++){
        traj_len[0] = MCtime;
        resampling(ProbMatrix, temp_micro, nMicro, MCtime);
        //trim the trajectory, remove the first 30% sampling points
        double remove_ratio=0.2; //input: remove the first 20% of samples
        int trim_traj_len[0];
        trim_traj_len[0]=(1-remove_ratio)*traj_len[0];
        int *new_temp_micro=alloarray_int(trim_traj_len[0]);
        int start_frame=remove_ratio*traj_len[0];
        ofstream myfile ("TrimmedMicroTrajectory.txt");
	if (myfile.is_open()){
        for (int j=0;j<trim_traj_len[0];j++){
			new_temp_micro[j]=temp_micro[j+start_frame];
			myfile << new_temp_micro[j] << "\n";
		}
	}
	else {
	cout << "Unable to open file"; 
        for (int j=0;j<trim_traj_len[0];j++){
                        new_temp_micro[j]=temp_micro[j+start_frame];
                }

	}
	
        int *temp_macro=alloarray_int(trim_traj_len[0]);
        micro2macro(new_temp_micro, mapping, temp_macro, trim_traj_len[0]);
        cout<<"We have removed the first 20% frames of MCMC traj, and now calculate based on "<<trim_traj_len[0]<<" frames"<<endl;
        mfpt_direct_count(temp_macro, trim_traj_len, 1, nMacro, mfpt_matrix);

        cout<<"now we are calculating the population of each macrostate in the MC chain"<<endl;
        for (int kk=0;kk<trim_traj_len[0];kk++){
            count_num[temp_macro[kk]]++; 
            total_count++;
        }
        cout<<"now outputing the count number and stationary population for each state"<<endl;
        free(temp_macro);
    }
    for (int kk=0;kk<nMacro;kk++)
        {cout<<"For state "<<kk<<":, count number: ,"<<count_num[kk]<<endl;}
        cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
    free(count_num);
    for (j = 0; j<nMicro;j++){
        free(ProbMatrix[j]);
    }
    for(j = 0; j<nMacro; j++){
        free(mfpt_count[j]);
        free(mfpt_time[j]);
    }
    free(ProbMatrix);
    free(mfpt_count);
    free(mfpt_time);

    //calculating mfpt of this macrostate chain
    cout<<"the mean first passage time is mean, unit:micro-second, from row to column"<<endl;
    for (int j = 0; j<nMacro; j++){
        for (int k=0; k<nMacro; k++){
            if(MCtime == 0) {
                cout<<mfpt_matrix[k][j]*timeunit<<'\t';
            }
            else{
                cout<<mfpt_matrix[k][j]*timeunit*lagtime<<'\t';
            }
        }
        cout<<endl;
    }
    //free matrices
    t = clock()-t;
    cpu_time_used = ((double)(t))/CLOCKS_PER_SEC;
    cout<<"****************************************************"<<endl;
    cout<<"Using a total of "<<cpu_time_used<<" seconds"<<endl;
    for(int j=0;j<nMacro;j++){
        free(mfpt_matrix[j]);
    }
    free(macro);
    free(temp_micro);
    return 0;
 }

