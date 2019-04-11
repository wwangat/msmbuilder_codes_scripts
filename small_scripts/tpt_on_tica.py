############################################################################
###      get top paths and plot the paths on tica projection             ###
###     tpt analysis is done by pyemma,       documentation of           ###
###     tpt analysis in msmbuilder is missing                            ###
############################################################################

import os,sys,subprocess
import numpy as np
import matplotlib as mp

import pyemma
from matplotlib import cm
from matplotlib import pyplot as plt
from msmbuilder.dataset import dataset
from msmbuilder.io import load_generic

ktrajs_dir = 'ktrajs-extracted-kcenters-lag1500-2-900'
ttrajs_dir = 'ttrajs-extracted-lag1500-new'
ktrajs_pkl = '%s-mle.pickl'%ktrajs_dir
micro_pkl  = 'msm-%s-mle.pickl'%ktrajs_dir
macro_pkl  = 'msm-%s-pcca-mle.pickl'%ktrajs_dir

clusterer  = load_generic("../%s"%ktrajs_pkl)
msm_model  = load_generic("../%s"%micro_pkl)
pcca       = load_generic("../%s"%macro_pkl)

############################################################################
###      get macrostate mapping and state_label of each microstate       ###
############################################################################
mapping = pcca.microstate_mapping_
state_label = msm_model.state_labels_
msm_label   = msm_model.state_labels_
n_microstates = msm_model.n_states_
n_macrostates = pcca.n_macrostates


############################################################################
###      get first two tica and plot free energy by mle msm_model        ###
############################################################################

mlemsm_pkl  = 'msm-ktrajs-extracted-kcenters-lag1500-2-900-mle.pickl'
mlemsm_model = load_generic('../%s'%mlemsm_pkl)

ktrajs = dataset('../%s'%ktrajs_dir,mode='r',fmt='dir-npy',verbose=True)
ktrajs = [ktrajs[s].tolist() for s in range(len(ktrajs))]

ttrajs = dataset('../%s'%ttrajs_dir,mode='r',fmt='dir-npy',verbose=True)
tvalues = [ttrajs[s].tolist() for s in range(len(ttrajs)) if s!=9]
tvaluesnew = []
for i in range(108):
    ti = []
    for j in range(9001):
        if ktrajs[i][j] in msm_label:
            ti += [tvalues[i][j]]
    tvaluesnew += [np.array(ti)]
tvaluesnew = np.array(tvaluesnew)
data = np.concatenate(tvaluesnew)

assignments = mlemsm_model.fit_transform(ktrajs)
pi_0 = mlemsm_model.populations_[np.concatenate(assignments, axis=0)]
#assignments = msm_model.fit_transform(ktrajs)
#pi_0 = msm_model.populations_[np.concatenate(assignments, axis=0)]

fig, ax = plt.subplots()
pyemma.plots.plot_free_energy(*data[:,:2].T, ax=ax, weights=pi_0, logscale=True)


############################################################################
###      get top path by pyemma tpt and plot them on tica space          ###
############################################################################
source = 0
sink   = 3

set1 = np.where(mapping == source)[0]
set2 = np.where(mapping == sink)[0]
print(set1)
print(set2)

#TPM = pyemma.msm.markov_model(msm_model.transmat_)

#tpt = pyemma.msm.tpt(TPM, set1, set2)


### get path from file ###
fpath = 'paths.txt'
p = open(fpath,'r')
lines = p.readlines()
p.close()

paths = []
pathfluxes = []
for line in lines:
    pathi = line.split()
    pathfluxes += [float(pathi[0])]
    paths += [[int(s) for s in pathi[1:]]]
pathfluxes = np.array(pathfluxes)
paths = np.array(paths)

#(paths,pathfluxes) = tpt.pathways(fraction=0.95,maxiter=2000)

maxflux = np.max(pathfluxes)
minflux = np.min(pathfluxes)
n_paths = 1681
for i in range(340,n_paths,10):
     pathi = paths[i]
     #print(pathi)
     pathi_real = np.array([state_label[s] for s in pathi])
     tica_pathi = clusterer.cluster_centers_[pathi_real]  
     #print(pathi_real)
     #print(tica_pathi)
     #flux_i = (pathfluxes[i])/(maxflux)
     flux_i = (pathfluxes[i]-minflux)/(maxflux-minflux)
     #ax.plot(tica_pathi[:,0],tica_pathi[:,1],c=cm.bwr(i*1.0/n_paths)) 
     ax.plot(tica_pathi[:,0],tica_pathi[:,1],c=cm.binary(flux_i)) 

plt.show()
############################################################################
###                              END                                     ###
############################################################################


