from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.cluster import KCenters
import numpy as np
from msmbuilder.decomposition import tICA
from msmbuilder.dataset import dataset
from matplotlib import pyplot as plt
import os
import sys
from matplotlib.colors import LogNorm

def fit_and_plot(transformed,opath, tIC_a, tIC_b):  #plot the tica projections
    transformed=np.concatenate(transformed)
    x = transformed[:, int(tIC_a)-1]
    y = transformed[:, int(tIC_b)-1]

    plt.axes(axisbg='w')
    plt.grid(False)
    plt.hist2d(x, y, bins=100, cmap='hot_r', norm=LogNorm())
    plt.xlabel('%s tIC'%(str(tIC_a)))
    plt.ylabel('%s tIC'%(str(tIC_b)))
    plt.title('tICA Heatmap (log color scale)')
    plt.colorbar()
    plt.savefig(opath)
    plt.close()


#####################begin to main program

#################3#######inputs
atom_pairs = np.loadtxt('pairlist.txt', dtype=int)   #indexes for the atom pairs you are interestd(index starts from 0): atom1 atom2
xtc_file_dir='trajectories/'  #folder to put xtc

featurizer = AtomPairsFeaturizer(pair_indices=atom_pairs)

traj_list_array=[]
for line in open("trajlist"):
    traj_list_array.append(line.strip())
print traj_list_array   #trajectory name


####################calculate the pairwise distances for tica
ticadist=[]
for trajfile in traj_list_array:
    xyz=dataset(xtc_file_dir+trajfile, topology='test.pdb')
    temp=featurizer.fit_transform(xyz)
    ticadist.append(temp[0])    #now we have the pairwise distance between the atoms of interest

###apart from the pairwise distances, other features you can try are "dihedral angle", "raw coordinates", "reciporacal distances", "contact map" and so on


###scan the tica parameters
'''   
lag_time_list=range(1,20,1)
for lag_time in lag_time_list:
    tica_model=tICA(lag_time=lag_time,n_components=10)
    tica_trajs=tica_model.fit(ticadist)  #projected tica coordinate
    print lag_time, tica_trajs.timescales_[0],tica_trajs.timescales_[1],tica_trajs.timescales_[2],tica_trajs.timescales_[3],tica_trajs.timescales_[4],tica_trajs.timescales_[5],tica_trajs.timescales_[6],tica_trajs.timescales_[7],tica_trajs.timescales_[8],tica_trajs.timescales_[9]
    np.savetxt('tic_eigenvector_lagtime_%d.txt'%(lag_time), tica_trajs.eigenvectors_)
    tica_trajs = tica_model.fit_transform(ticadist)
    fit_and_plot(tica_trajs, 'tICA_projection_'+str(lag_time)+'.png', 1,2)
'''
###


##########################calculate the tica projections, tica eigenvecor,tica eigenvectors and eigenvalues 
os.system("mkdir tica_projections/")

lag_time = 10  #tica lag time
tica_model=tICA(lag_time=lag_time,n_components=10)  #n_components: number of tica projections it will save into RAM

tica_trash = tica_model.fit(ticadist) #"fit" only calculate tica eigenvector (right eigenvector) and eigenvalue
np.savetxt('tica_eigenvector_lagtime_%d.txt'%(lag_time), tica_trash.eigenvectors_)
np.savetxt('tica_eigenvalue_lagtime_%d.txt'%(lag_time), tica_trash.eigenvalues_)
np.savetxt('tica_covariance_matrix.txt', tica_trash.covariance_)   #covariance matrix W, in tica we have v'*W*v = I; v is the right eigenvectors

tica_trajs=tica_model.fit_transform(ticadist)  #fit_transform: only return tica projections on the major eigenvectors

#tica_trajs is the data for your DBSCAN

for j in range(len(traj_list_array)):
    np.savetxt('tica_projections/%s_tica.txt'%(traj_list_array[j][:-4]), tica_trajs[j][:,0:4])   #we save the top tica projections

###after you have the tica projection (tica_trajs), you can do clustering like DBSCAN, KCENTERS


exit()



#########################below: kcenters to get the microstates
nMicro = 100
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



#below: perform pcca plus to get the macrostates
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

