"""
Author: Wei WANG (wwangat@gmail.com)
Function: do pcca plus and project conformations onto major tICs
input: MicroAssignment list and tICA coordinates list

"""
from pyemma.msm import PCCA #actually it is pcca+
from pyemma import msm
import numpy as np
from msmbuilder.cluster import KCenters
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt 
import os
from msmbuilder.msm import MarkovStateModel
import sys 
from matplotlib.colors import LogNorm


kcenters_sequences=[]
#input:trajlist_micro: path to all the microstate chains
traj_num=0
trajlist_name=[]
for line in file('trajlist_micro'):
    trajlist_name.append(line.strip())
    traj_num=traj_num+1
    temp = np.loadtxt(line.strip(), dtype=int)
    kcenters_sequences.append(temp.tolist())

#print kcenters_sequences

tic_12 = []

#input: trajlist_tic: path to all the tic projections

for line in file('trajlist_tic'):
    temp = np.loadtxt(line.strip())
    tic_12.append(temp.tolist())
#print tic_12

#input: microstate_lagtime for buiding a markovian model
microstate_lagtime = 40
#input: nMacro: number of macrostates we use
nMacro = 6

#msm=MarkovStateModel(verbose=True, lag_time=microstate_lagtime, reversible_type='transpose', ergodic_cutoff='off')
msm=MarkovStateModel(verbose=True, lag_time=microstate_lagtime, reversible_type='mle', ergodic_cutoff='on')
msm.fit(kcenters_sequences)

np.savetxt("kCenters_stationary_population", msm.populations_)
#TCM = np.loadtxt('micro_tcm.txt')
#nMacro = 6

#nMicro = TCM.shape[0]
#print nMicro
#TPM = TCM

#for j in range(nMicro):
#    rowsum = np.sum(TCM[j,:])
#    for k in range(nMicro):
#        TPM[j][k]=TCM[j][k]/rowsum 
#print TPM
#print "having calculated the row normalized TPM"

pcca_result = PCCA(msm.transmat_, nMacro)
print "now writing the pcca+ assignment"

np.savetxt("pcca_plus_%d_state_mapping.txt"%(nMacro), pcca_result.metastable_assignment, fmt='%d')
np.savetxt("pcca_plus_%d_state_stationary_population.txt"%(nMacro), pcca_result.coarse_grained_stationary_probability)
np.savetxt("pcca_plus_%d_state_transmat.txt"%(nMacro), pcca_result.coarse_grained_transition_matrix)

#begin to save "trajectory_name, frame_number, microstate assignment, macrostate assignment"
print traj_num
f=open("pcca_plus_%d_state_for_selecting_confor.txt"%(nMacro), 'w')
f.write("trajlist_name    frame_number    microstate_assignment    macrostate_assignment\n")

print msm.mapping_
for j in range(traj_num):
 #   print kcenters_sequences[j]
    for k in range(len(kcenters_sequences[j])):
        if(msm.mapping_.has_key(kcenters_sequences[j][k])):
            stateid=msm.mapping_[kcenters_sequences[j][k]]
            f.write("%s    %d    %d    %d\n" % (trajlist_name[j], k, kcenters_sequences[j][k], pcca_result.metastable_assignment[stateid]))

f.close()

nMicro = msm.n_states_

#drawing projections on major tics
print "begin drawing figure"
tic_transformed=np.concatenate(tic_12)
micro_transformed = np.concatenate(kcenters_sequences)

macro_transformed = micro_transformed.copy()
macro_transformed[:] = -1
print nMicro
for j in range(700):  #this number may change afterwards
    if (msm.mapping_.has_key(j)):
        macro_transformed[np.where(micro_transformed == j)] = pcca_result.metastable_assignment[msm.mapping_[j]]
colors_jet = plt.cm.jet(np.linspace(0,1,nMacro))

for m in range(3):
    for n in range(m+1,3):
        print "Now we begin to deal with projection on tic ", m, " and tic ", n
        plt.figure()
        for j in range(nMacro):
            j=nMacro-j-1
            index = np.where(macro_transformed == j)
            index1=index[0][: : 10]
#    print index1
            print tic_transformed[index1,m]
            print tic_transformed[index1,n]
            plt.plot(tic_transformed[index1, m], tic_transformed[index1, n], 'x', label="%d"%(j), color=colors_jet[j])
        plt.legend(loc=1, labelspacing=0.075, prop={'size': 8.0}, scatterpoints=1, markerscale=0.8, numpoints=1)
        plt.title("pcca_plus_%d_state"%(nMacro))
        plt.xlabel("tIC%d"%(m))
        plt.ylabel("tIC%d"%(n))
        plt.savefig("pcca_plus_%d_state_projection_tic_%d%d.png"%(nMacro, m, n))
        plt.close()
