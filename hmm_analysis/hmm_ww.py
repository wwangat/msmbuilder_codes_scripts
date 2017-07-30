"""
Author:Wei WANG (wwangat@gmail.com)
Function:applying the hidden markov model
Dependency: pyemma
Input: non-markovian microstate sequences
Output: emisson, stationary population, transition probability matrix, membership function, Implied timescale plot
"""
from pyemma import msm
import numpy as np
from msmbuilder.cluster import KCenters
import matplotlib
from msmbuilder.msm import MarkovStateModel

kcenters_sequences=[]
#input:trajlist_micro: path to all the microstate chains
traj_num=0
trajlist_name=[]

for line in open('trajlist_micro'):
    trajlist_name.append(line.strip())
    traj_num=traj_num+1
    temp = np.loadtxt(line.strip(), dtype='int')
    kcenters_sequences.append(temp.tolist())

nMacro=4
#microstate_lagtime=50
reversible_type='True'  #If true compute reversible MSM, else non-reversible MSM

connectivity_type='largest'

initial=10
ending=400
interval = ending*1.0/20;
lag_times=[]
for j in range(20):
    lag_times.append(int(initial+j*interval))

print lag_times
print "unit is 20 ps"
for j in range(20):
    print "now we are dealing with lagtime ", lag_times[j]
    mm=msm.estimate_hidden_markov_model(kcenters_sequences, nMacro, lag=lag_times[j], reversible=reversible_type,connectivity=connectivity_type)
    np.savetxt("hmm_%d_state_%d_lagtime_transition_matrix.txt"%(nMacro, lag_times[j]), mm.transition_matrix)
    np.savetxt("hmm_%d_state_%d_lagtime_stationary_pop.txt"%(nMacro, lag_times[j]), mm.stationary_distribution)
    np.savetxt("hmm_%d_state_%d_lagtime_emisson.txt"%(nMacro, lag_times[j]), mm.metastable_distributions)
    np.savetxt("hmm_%d_state_%d_lagtime_membership.txt"%(nMacro, lag_times[j]), mm.metastable_memberships)
    print mm.timescales()
