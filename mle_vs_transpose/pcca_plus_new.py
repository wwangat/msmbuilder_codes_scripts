from pyemma.msm import PCCA #actually it is pcca+
from pyemma import msm
import numpy as np
from msmbuilder.cluster import KCenters
import matplotlib
matplotlib.use('Agg')
from msmbuilder.msm import implied_timescales
from msmbuilder.msm import MarkovStateModel

import matplotlib.pyplot as plt 
import os
from msmbuilder.msm import MarkovStateModel
import sys 
from matplotlib.colors import LogNorm

def fit_and_plot(input_data, tIC_a, tIC_b):
    x = input_data[:, int(tIC_a)]
    y = input_data[:, int(tIC_b)]

    plt.axes(axisbg='w')
    plt.grid(False)
    plt.hist2d(x, y, bins=100, cmap='hot_r', norm=LogNorm())
    plt.xlabel('%s tIC'%(str(tIC_a)))
    plt.ylabel('%s tIC'%(str(tIC_b)))
    plt.title('tICA Heatmap (log color scale)')
    plt.colorbar()

kcenters_sequences=[]
#input:trajlist_micro: path to all the microstate chains
traj_num=0
trajlist_name=[]
for line in open('trajlist_micro'):
    trajlist_name.append(line.strip())
    traj_num=traj_num+1
    temp = np.loadtxt(line.strip())
    kcenters_sequences.append(temp.tolist())

#print kcenters_sequences

#tic_12 = []

#input: trajlist_tic: path to all the tic projections

#for line in open('trajlist_tic'):
#    temp = np.load(line.strip())
#    tic_12.append(temp.tolist())
#print tic_12

#input: microstate_lagtime for buiding a markovian model
#input: nMacro: number of macrostates we use
#nMacro = 4
#tic_transformed=np.concatenate(tic_12)
#micro_transformed = np.concatenate$(kcenters_sequences)

#tic1=0
#tic2=1
#print(tic_transformed)

microstate_lagtime=50
reversible='none'

initial=10
ending=400
interval=pow(ending*1.0/initial,1.0/20);
lag_times=[]
for j in range(20):
    lag_times.append(initial*pow(interval,j))
#lag_times=range(10,100,10)
msm=MarkovStateModel(verbose=True, lag_time=microstate_lagtime, reversible_type=reversible, ergodic_cutoff='on')
msm.fit(kcenters_sequences)
print msm.mapping_
print("for microstate lag time = ", microstate_lagtime, ",", msm.n_states_, " states are left")

#for microstate_lagtime in range(10,20,10):
#    msm=MarkovStateModel(verbose=True, lag_time=microstate_lagtime, reversible_type='mle', ergodic_cutoff='on')
#    msm.fit(kcenters_sequences)
#    print("for microstate lag time = ", microstate_lagtime, ",", msm.n_states_, " states are left")

#    new_micro_transformed = -1*np.ones(len(micro_transformed))
#    dict_key = list(msm.mapping_.keys())
#    dict_value = list(msm.mapping_.values())
#    for j in range(msm.n_states_):
#        new_micro_transformed[np.where(micro_transformed == dict_key[j])] = dict_value[j]

#    print(new_micro_transformed)
#    plt.figure()
#    fit_and_plot(tic_transformed, tic1, tic2)
#    index1=np.where(new_micro_transformed == -1)
#    index1=index1[0]
#    plt.plot(tic_transformed[index1, tic1], tic_transformed[index1, tic2], 'x')
#    plt.savefig("microstate_assignment_lagtime_%d_trimed_%dstate_projection_tic.png"%(microstate_lagtime, msm.n_states_))
#    plt.close()

#    plt.figure()
#    fit_and_plot(tic_transformed, tic1, tic2)
#    index = np.where(new_micro_transformed != -1)
#    index1=index[0][: : 10]
#    plt.plot(tic_transformed[index1, tic1], tic_transformed[index1, tic2], 'x')
#    plt.savefig("microstate_assignment_lagtime_%d_remains_%dstate_projection_tic.png"%(microstate_lagtime, msm.n_states_))
#    plt.close()


np.savetxt("kcenters_microstate_%s_transmat_.txt"%(reversible), msm.transmat_)
np.savetxt("kcenters_%s_stationary_population"%(reversible), msm.populations_)

#plot implied timescale
n_timescales=10
print "lagtime list is:", lag_times
msm_timescales = implied_timescales(kcenters_sequences, lag_times, n_timescales=n_timescales, msm=MarkovStateModel(verbose=True, reversible_type=reversible, ergodic_cutoff='on'))

for k in range(n_timescales):
    plt.plot(lag_times, msm_timescales[:, k], 'o-')
plt.title('Discrete-time MSM Relaxation Timescales')
plt.semilogy()
x1,x2,y1,y2 = plt.axis()
#plt.axis((x1,x2,y1,1000000))  #need to change this number
outfile_name="%s_ITS.png"%(reversible)
plt.savefig(outfile_name)
plt.close()
outfile_name="%s_ITS.dat"%(reversible)
print msm_timescales
f2=open(outfile_name, 'w')
for i in range(len(lag_times)):
    f2.write("%d    "%(lag_times[i]))
    for j in range(n_timescales):
        f2.write("%f    "%(msm_timescales[i, j]))
    f2.write('\n')
f2.close()

"""

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

print("now writing the pcca+ assignment")

np.savetxt("pcca_plus_%d_state_mapping.txt"%(nMacro), pcca_result.metastable_assignment, fmt='%d')
np.savetxt("pcca_plus_%d_state_stationary_population.txt"%(nMacro), pcca_result.coarse_grained_stationary_probability)
np.savetxt("pcca_plus_%d_state_transmat.txt"%(nMacro), pcca_result.coarse_grained_transition_matrix)

"""
"""
#begin to save "trajectory_name, frame_number, microstate assignment, macrostate assignment"
print traj_num
f=open("pcca_plus_%d_state_for_selecting_confor.txt"%(nMacro), 'w')
f.write("trajlist_name    frame_number    microstate_assignment    macrostate_assignment\n")
for j in range(traj_num):
 #   print kcenters_sequences[j]
    for k in range(len(kcenters_sequences[j])):
        f.write("%s    %d    %d    %d\n" % (trajlist_name[j], k, kcenters_sequences[j][k], pcca_result.metastable_assignment[kcenters_sequences[j][k]]))

f.close()


nMicro = msm.n_states_

#drawing projections on major tics
print "begin drawing figure"
tic_transformed=np.concatenate(tic_12)
micro_transformed = np.concatenate(kcenters_sequences)

macro_transformed = micro_transformed.copy()
print nMicro
for j in range(nMicro):
    macro_transformed[np.where(micro_transformed == j)] = pcca_result.metastable_assignment[j]
colors_jet = plt.cm.jet(np.linspace(0,1,nMacro))
plt.figure()
for j in range(nMacro):
    j=nMacro-j-1
    print j
    index = np.where(macro_transformed == j)
    index1=index[0][: : 10]
#    print index1
    print tic_transformed[index1,0]
    print tic_transformed[index1,1]
    plt.plot(tic_transformed[index1, 0], tic_transformed[index1, 1], 'x', label="State %d" % j, color=colors_jet[j])

plt.title("pcca_plus_%d_state"%(nMacro))
plt.xlabel("tIC1")
plt.ylabel("tIC2")
plt.savefig("pcca_plus_%d_state_projection_tic.png"%(nMacro))
plt.close()
"""
