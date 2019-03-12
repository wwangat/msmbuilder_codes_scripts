#include<iostream>

using namespace std;
//part 1: max of an array
double max(double *array, int size){
	int i;
	double max=array[0];
	for (i=1;i<size;i++){
		if (array[i]>max)
			{max = array[i];}
	}
	return max;
}
//part 2: generate random number
int randi(int min, int max){
	return rand() % (max-min+1)+min;
} //generate random integer between min and max
double randd(){
	double randnum;
	randnum = rand()/(RAND_MAX+1.0);
	return randnum;
} //generate pseudo random number between 0 and 1
double posterior(double **tCount, int **nCount, double alpha,double beta, int sizeMacro){
//calculate gamma conjugate prior posterior probability
	double sum = 0.0, temp;
	int i, j;
	for (i = 0; i<sizeMacro; i++){
		for (j = i+1; j<sizeMacro; j++){
			//in this situation, since tCount starts from 0, we should add them by one
			temp = ceil(tCount[i][j]+alpha-1);
			sum += temp*log(temp)-temp-(tCount[i][j]+alpha)*log(nCount[i][j]+beta);
		}
	}
	return sum;
}
/*modularity */
double modularity(double *degree, double totaledge,int *array, double **matrix, int sizeMicro, int sizeMacro)
{
	int i, j, row, col;
	double sum = 0.0;
	for (i = 0; i<sizeMicro; i++){
		row = array[i];
		for (j = 0; j<sizeMicro; j++){
			col = array[j];
			if (row == col)
			{
				sum += matrix[i][j]-degree[i]*degree[j]/totaledge;
			}
		}
	}
	sum = sum/totaledge;
	return sum;
}
char* num2str(int num)
{
	static char str[10];
	int rem;
	int n, len=0, i;
	n = num;
	while (n!=0)
	{
		len++;
		n /= 10;
	}
	for (i=0;i<len;i++)
	{
		rem = num %10;
		num = num/10;
		str[len -(i+1)] = rem + '0';
	}
	str[len] = '\0';
	return str;
}
/*monte carlo*/
int MonteCarlo(double *array, int size){
	int index = 0, i;
	double randd();
	double randnum, sum = array[0];
	for (i = 1;i<size;i++){
		sum += array[i];
	}
	randnum = sum*randd();
	if (randnum<=array[0])
		index = 0;
	else
	{
		sum = 0.0;
		for (i = 0; i<size-1; i++){
			if(sum+array[i]<randnum && sum+array[i]+array[i+1]>=randnum){//bug fix on July 3, 2016. add +array[i]
				index = i+1;
				break;
			}
			else
			{
				sum = sum+array[i];
			}

		}	
	}
	return index;
}

/* resampling based on transition probability matrix
 * the transition probability matrix is column normalized
  This function is validated to be right
*/

//this one is highly correlated since using sliding window
void resampling(double **ProbMatrix, int *chain, int dim, int MCstep){
    chain[0] = randi(0, dim-1);
    //check whether the input ProbMatrix is column normalized
    int j, k;
    double temp_sum = 0.0;
    for (j = 0; j<dim; j++){
        temp_sum += ProbMatrix[j][0];
    }
    if(temp_sum > 1.1 || temp_sum < 0.9){
        cout<<"input transition probability matrix should be column normalized, exiting..."<<endl;
        exit(0);
    }
    double randd();
    double randnum;
    for (j = 1; j<MCstep; j++){
        if (j%100 == 0) fprintf(stderr, "resampling: %d states sampled\r", j);
        randnum = randd();
        if (randnum <=ProbMatrix[0][chain[j-1]]){
            chain[j] = 0;
        }
        else{
            temp_sum = 0.0;
            for (k = 0; k<dim-1; k++){
                if(temp_sum+ProbMatrix[k][chain[j-1]]<randnum && temp_sum+ProbMatrix[k][chain[j-1]]+ProbMatrix[k+1][chain[j-1]]>=randnum){
                    chain[j] = k+1;
                    break;
                }
                else{
                    temp_sum += ProbMatrix[k][chain[j-1]];
                }
            }
        }
    }
    
}

/*resampling based on transition probability matrix, but sampling would stop when traversing all the macrostates
nMicro: microstate dimension, return macrostate chain, nMacro: macrostate dimension
*/



//this one is slightly co-related
int resampling(double **ProbMatrix, int *chain, int *mapping, int nMicro, int nMacro, int MCstep){
    int tempid = randi(0, nMicro-1);
    chain[0] = mapping[tempid];
    //check whether the input ProbMatrix is column normalized
    int j, k;
    int *flag=alloarray_int(nMacro);
    for (j=0;j<nMacro;j++){
        flag[j] = 0;
    }
    int chain_length = 1;
    double temp_sum = 0.0;
    for (j = 0; j<nMicro; j++){
        temp_sum += ProbMatrix[j][tempid];
    }
    if(temp_sum > 1.1 || temp_sum < 0.9){
        cout<<"input transition probability matrix should be column normalized, exiting..."<<endl;
        exit(0);
    }
    flag[chain[0]]=1;//mark macrostate flag of initial state, then what we will do is to check other macrostate flag
    double randd();
    double randnum;
    for (j = 1; j<MCstep; j++){
//        if (j%100 == 0) fprintf(stderr, "resampling: %d states sampled\r", j);
        randnum = randd();
        if (randnum <=ProbMatrix[0][tempid]){
            tempid = 0;
        }
        else{
            temp_sum = 0.0;
            for (k = 0; k<nMicro-1; k++){
                if(temp_sum+ProbMatrix[k][tempid]<randnum && temp_sum+ProbMatrix[k][tempid]+ProbMatrix[k+1][tempid]>=randnum){
                    tempid = k+1;
                    break;
                }
                else{
                    temp_sum += ProbMatrix[k][tempid];
                }
            }
        }
        //now we have sampled one more step, now checking whether we should stop the chain
        chain_length++;
        chain[j] = mapping[tempid];
        flag[chain[j]] = 1;
        int sum_flag=0;
        for (k=0;k<nMacro;k++){
            sum_flag += flag[k];
        }
        if(sum_flag==nMacro){
            break;
        }
    }
    free(flag);
    return chain_length;
}

//this one is non-correlated, resampling and updating the matrix on the fly, repeat_time is the number of synthetic trajectories
void uncorrelated_mfpt(double **ProbMatrix, int *mapping, int nMicro, int nMacro, int MCstep, int repeat_time, double **mfpt_time, double **mfpt_count){
    int old_state, new_state, tempid, j, k;
    int loop_time = 0;
    while (loop_time<repeat_time){
        tempid = randi(0, nMicro-1);
        old_state = mapping[tempid];
        new_state = randi(0, nMacro-1);
        if(old_state == new_state){
            continue;
        }
        else{
            loop_time++;
            int chain_length = 0; //update this number to mfpt_count when stopping
            //check whether the input ProbMatrix is column normalized
            double temp_sum = 0.0;
            for (j = 0; j<nMicro; j++){
                temp_sum += ProbMatrix[j][old_state];
            }
            if(temp_sum > 1.1 || temp_sum < 0.9){
                cout<<"input transition probability matrix should be column normalized, exiting..."<<endl;
                exit(0);
            }
            double randd();
            double randnum;
            for (j = 1; j<MCstep; j++){
                randnum = randd();
                if (randnum <=ProbMatrix[0][tempid]){
                    tempid = 0;
                }
                else{
                    temp_sum = 0.0;
                    for (k = 0; k<nMicro-1; k++){
                        if(temp_sum+ProbMatrix[k][tempid]<randnum && temp_sum+ProbMatrix[k][tempid]+ProbMatrix[k+1][tempid]>=randnum){
                            tempid = k+1;
                            break;
                        }
                        else{
                            temp_sum += ProbMatrix[k][tempid];
                        }
                    }
                }
                //now we have sampled one more step, now checking whether we should stop the chain
                chain_length++;
                if(mapping[tempid] == new_state){
                    break;
                }
            }
            //mfpt_time: add intervals(traj_len), mfpt_count: add 1 
            mfpt_time[new_state][old_state] += chain_length;
            mfpt_count[new_state][old_state]++;
            if(loop_time%10000==0){
                cout<<"the "<<loop_time<<" th loop is from "<<old_state<<" to "<<new_state<<" with length "<<chain_length+1<<endl;
            }
        }
    }
}
