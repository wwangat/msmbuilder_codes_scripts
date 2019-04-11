from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.cluster import KCenters
import numpy as np
from msmbuilder.decomposition import tICA
from msmbuilder.dataset import dataset
from matplotlib import pyplot as plt
import matplotlib as mp
import os
from msmbuilder.tpt import net_fluxes
import sys
#xyz = dataset('../xtc/random-150/*.xtc', topology = '~/Desktop/tica-projection/Structures/Reference-PRE.pdb')
#list1=np.loadtxt('atompairs-5pairs-5helix-P')
#featurizer = AtomPairsFeaturizer(pair_indices=list1)

#ticadist = xyz.fit_transform_with(featurizer, 'atompairsfeaturizer/', fmt='dir-npy')

#ticadist =dataset('./atompairsfeaturizer/',mode='r',fmt='dir-npy',verbose=True)
#tica_model=tICA(lag_time=400,n_components=2)
#tica_model=ticadist.fit_with(tica_model)
#tica_trajs = ticadist.transform_with(tica_model, 'tica-150/',fmt='dir-npy')
#tica_trajs=dataset('./tica',mode='r',fmt='dir-npy',verbose=True)

#txx = np.concatenate(tica_trajs)
#plt.hexbin(txx[:,0], txx[:,1],bins='log', mincnt=0.1, cmap='viridis')
#data=np.load('./tica-301/00000300.npy')
#plt.plot(data[:,0], data[:,1],color = 'red', linewidth = 3)
#data1=np.load('./tica-301/00000301.npy')
#plt.plot(data1[:,0], data1[:,1],color = 'black', linewidth = 3)
#plt.savefig('tica-150random.png')

#clusterer = KCenters(n_clusters=1000,random_state=8)
#clusterer = dataset('./ktrajs-rawpos-dp/',mode='r',fmt='dir-npy',verbose=True)
#clustered_trajs = tica_trajs.fit_transform_with(clusterer, './ktrajs-rawpos-dp/', fmt='dir-npy')
#np.savetxt('center.txt',clusterer.distances_)
#print(len(clusterer.distances_))
#for i in range(300):
#	np.savetxt('./center/traj'+str(i) +'.txt',clusterer.distances_[i])

from msmbuilder.io import load_trajs, save_trajs, save_generic
from msmbuilder.cluster import MiniBatchKMeans
#from msmbuilder.cluster import KCenters

from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import dump
from msmbuilder.dataset import dataset
import os
import numpy as np

import pyemma
from msmbuilder.lumping import PCCAPlus,PCCA
#from msmbuilder.lumping import PCCAPlus
from msmbuilder.io import load_trajs, load_generic

tlag = 400
ntrajs = 130
ktrajs_dir = 'ktrajs-extracted-kcenters-lag1500-2-1000'
ktrajs_pkl = '%s-mle.pickl'%ktrajs_dir
microtraj_dir = 'microktrajs-%s-mle8'%ktrajs_dir
microtraj_pkl = 'msm-%s-mle8.pickl'%ktrajs_dir
ttrajs_dir = 'ttrajs-atpair2-lag1500'
#unwrapbp_dir = '../msm-remove10ns-cluster-unwrapbpnew2/unwrapbpnew2'

msm = load_generic(microtraj_pkl)

#meta, k_trajs = load_trajs(ktrajs_dir)
ktrajs = dataset('%s-mle'%ktrajs_dir,mode='r',fmt='dir-npy',verbose=True)
ttrajs = dataset(ttrajs_dir,mode='r',fmt='dir-npy',verbose=True)
#utrajs = dataset(unwrapbp_dir,mode='r',fmt='dir-npy',verbose=True)
ktrajs = [ktrajs[s].tolist() for s in range(len(ktrajs))]
#tvaluesall = [ttrajs[s].tolist() for s in range(len(ttrajs)) if s!=9]
#uvaluesall = [utrajs[s].tolist() for s in range(len(utrajs)) if s!=9]
tvaluesall = [ttrajs[s].tolist() for s in range(len(ttrajs))]
#uvaluesall = [utrajs[s].tolist() for s in range(len(utrajs))]
#tvalues = [ttrajs[s] for s in range(len(ttrajs))]
#txx = np.concatenate(tvaluesnew)

msm_label = msm.state_labels_
print(len(msm_label))
#uvaluesnew = []
#for i in range(350):
#    ti = []
#    for j in range(len(ktrajs[i])):
#        if ktrajs[i][j] in msm_label:
#            ti += [uvaluesall[i][j]]
#    uvaluesnew += [np.array(ti)]
#uvaluesnew = np.array(uvaluesnew)
tvaluesnew = []
for i in range(ntrajs):
    ti = []
    for j in range(len(ktrajs[i])):
        if ktrajs[i][j] in msm_label:
            ti += [tvaluesall[i][j]]
    tvaluesnew += [np.array(ti)]
tvaluesnew = np.array(tvaluesnew)

#clusterer = dataset(ktrajs_dir,mode='r',fmt='dir-npy',verbose=True)
#clusterer = load_generic(ktrajs_pkl)

fig, ax = plt.subplots()
#import msmexplorer as msme
#assignments = clusterer.partial_transform(txx)
#kxx = np.concatenate(ktrajs, axis=0)
#assignments = msm.fit_transform(kxx)
assignments = msm.fit_transform(ktrajs)
#data = np.concatenate(assignments, axis=0)
#data = np.concatenate(ttrajs, axis=0)
#data = np.concatenate(uvaluesnew, axis=0)
data = np.concatenate(tvaluesnew, axis=0)
pi_0 = msm.populations_[np.concatenate(assignments, axis=0)]

print(len(data),len(pi_0))

#plt.plot(data[:,0],data[:,1])
#ax = msme.plot_free_energy(data,obs=(0,1), n_samples=100000,pi=pi_0,
#                      shade=True,
#                      clabel=True,
#                      clabel_kwargs={'fmt': '%.1f'},
#                      #vmax=5,
#                      #vmin=0,
#                      gridsize=31,cut=5,
#                      cbar=True,
#                      cbar_kwargs={'format': '%.1f', 'label': 'Free energy (kcal/mol)'})
#pyemma.plots.plot_free_energy(*data[:,:2].T, ax=ax, weights=pi_0, nbins=30,avoid_zero_count=True,logscale=True)
pyemma.plots.plot_free_energy(*data[:,:2].T, ax=ax, weights=pi_0, logscale=True)
#ax.set_xlim(0,28)
#ax.set_ylim(0,26)
#ax = msme.plot_free_energy(data, n_samples=10000,pi=pi_0)
                      #pi=pi_0,xlabel='tIC 1', ylabel='tIC 2')

#msme.plot_histogram(data, color='oxblood', quantiles=(0.5,),
#                    labels=['$tIC1$', '$tIC2$' ],
#                    show_titles=True)

#plt.hexbin(txx[:, 0], txx[:, 1], bins='log', mincnt=0.1, cmap="Greys")
#plt.scatter(clusterer.cluster_centers_[msm.state_labels_, 0],
#            clusterer.cluster_centers_[msm.state_labels_, 1],
#            s=5,
#            c=pcca.microstate_mapping_,
#)
#plt.hexbin(clusterer.cluster_centers_[msm.state_labels_, 0],
#            clusterer.cluster_centers_[msm.state_labels_, 1],
#            s=50,
#            cmap=pcca.microstate_mapping_,
#)
#plt.plot(data[:,0], data[:,1],color = 'red', linewidth = 1,marker='o',markersize=3)
#plt.plot(data1[:,0], data1[:,1],color = 'green', linewidth = 1,marker='o',markersize=3)
plt.xlabel('tICA 1')
plt.ylabel('tICA 2')
plt.savefig('fep-pyemma-tica.png')
plt.show()

#np.savetxt('pcca_4.txt', pcca.microstate_mapping_)
