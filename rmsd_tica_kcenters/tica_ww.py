from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.cluster import KCenters
import numpy as np
from msmbuilder.decomposition import tICA
from msmbuilder.dataset import dataset
from matplotlib import pyplot as plt
import os
import sys

atom_pairs = np.loadtxt('pairlist.txt', dtype=int)  #atom pair, starting from 0
nMicro=36
xtc_file_dir='../trajectories/'

featurizer = AtomPairsFeaturizer(pair_indices=atom_pairs)

traj_list_array=[]
for line in open("trajlist"):
    traj_list_array.append(line.strip())
print traj_list_array

ticadist=[]
for trajfile in traj_list_array:
    xyz=dataset(xtc_file_dir+trajfile, topology='../test.pdb')
    temp=featurizer.fit_transform(xyz)
    ticadist.append(temp[0])

tica_model=tICA(lag_time=20,n_components=4)

tica_trajs=tica_model.fit_transform(ticadist)  #projected tica coordinate

os.system("mkdir tica_projections/")

for j in range(len(traj_list_array)):
    np.savetxt('tica_projections/%s_tica.txt'%(traj_list_array[j][:-4]), tica_trajs[j][:,0:1])



kcenters=KCenters(n_clusters=nMicro, metric='euclidean', random_state=0)

kcenters_sequences = kcenters.fit(tica_trajs)

out_assignment_dir='Microassignment/'
os.system("mkdir %s"%(out_assignment_dir))

tmp_counter=0
for ifile in traj_list_array:
    np.savetxt("%s/%s_assignment_.txt"%(out_assignment_dir, ifile[:-4]),kcenters.labels_[tmp_counter],fmt='%d')
    tmp_counter+=1

exit()
from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import dump

#msm = MarkovStateModel(lag_time=10, reversible_type='transpose')
#msm.fit(clustered_trajs)

#print(len(msm.state_labels_))


#print("now output state label")
#print(msm.state_labels_)
#np.savetxt('state_lables.txt', msm.state_labels_)


#print("now output TPM")
#print(msm.transmat_)
#np.savetxt('TPM.txt', msm.transmat_)

#print("now output populations")
#print(msm.populations_)
#np.savetxt("population.txt", msm.populations_)



from msmbuilder.msm import MarkovStateModel, implied_timescales

data=dataset('./cluster',mode='r',fmt='dir-npy',verbose=True)

lag_times=range(10,400,10)
msm_timescales = implied_timescales(data, lag_times, n_timescales=10,msm=MarkovStateModel(reversible_type='transpose'))
np.savetxt('msm_timescales_XXX.txt',msm_timescales)


exit()

from msmbuilder.lumping import PCCAPlus
pcca = PCCAPlus.from_msm(msm, n_macrostates=3)
macro_trajs = pcca.transform(clustered_trajs)

plt.hexbin(txx[:, 0], txx[:, 1], bins='log', mincnt=0.1, cmap="bone_r")
plt.scatter(clusterer.cluster_centers_[msm.state_labels_, 0],
            clusterer.cluster_centers_[msm.state_labels_, 1],
            s=5,
            c=pcca.microstate_mapping_,
)
plt.xlabel('tIC 1')
plt.ylabel('tIC 2')
plt.savefig('3macro.png')

np.savetxt('pcca_3.txt', pcca.microstate_mapping_)

