from pyemma import msm
import numpy as np
import matplotlib
matplotlib.use('Agg')
from msmbuilder.msm import implied_timescales
from msmbuilder.msm import MarkovStateModel

import matplotlib.pyplot as plt 
import os
from msmbuilder.msm import MarkovStateModel

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
