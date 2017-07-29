#usage: after running tICA2-tic-for_bootstrap.py and post_bootstrap_script.sh
import numpy as np
for tIC in range(1,20,1):
    a=np.loadtxt('for_stat_tIC_%d'%(tIC));
    #column 2-norm normalize
    print a.shape[0],a.shape[1]   #a.shape[1] is the # of columns
    for k in range(a.shape[1]):
        a[:,k]=a[:,k]/np.linalg.norm(a[:,k])  #frobenuis 2-norm
        a[:, k] = a[:, k]**2
    #take absolute values
    b=[]
    #now begin to calculate the min,max,mean, median, std of each row
    for k in range(a.shape[0]):
        array=a[k,:]
        b.append([np.min(array), np.max(array), np.mean(array), np.median(array), np.std(array)])
    np.savetxt('final_stat_tIC_%d_min_max_mean_median_std'%(tIC), b)
    
    
