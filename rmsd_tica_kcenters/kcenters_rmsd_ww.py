"""
Writer: Wei WANG (wwangat@gmail.com)
Function: kcenters, metric=aligned RMSD
"""

from msmbuilder.cluster import KCenters
from msmbuilder.featurizer import RMSDFeaturizer
import numpy
from msmbuilder.dataset import dataset
import matplotlib
import mdtraj as md
matplotlib.use('Agg')

import os
from msmbuilder.msm import MarkovStateModel

atomindex_file="AtomIndices.dat"   #input atom indices file, inputs, starting from 0

select_atoms = numpy.loadtxt(atomindex_file, dtype=int)

nMicro=100

print select_atoms
xtc_file_dir = '../trajectories/';   #inputs, where you put the xtc file

traj_list_array=[]
for line in open("trajlist"):
    traj_list_array.append(line.strip())
print traj_list_array

dataset=[]
for trajfile in traj_list_array:
    t = md.load(xtc_file_dir+trajfile, top='test.pdb', atom_indices=select_atoms)
    dataset.append(t)
print dataset
#ww: check whether they have aligned w.r.t reference

kcenters=KCenters(n_clusters=nMicro, metric='rmsd', random_state=0)

kcenters_sequences = kcenters.fit(dataset)

out_assignment_dir='Microassignment/'
out_kcenters_distances_dir='distances/'
os.system("mkdir %s"%(out_assignment_dir))
os.system("mkdir %s"%(out_kcenters_distances_dir))

tmp_counter=0
for ifile in traj_list_array:
    numpy.savetxt("%s/%s_assignment_.txt"%(out_assignment_dir, ifile[:-4]),kcenters.labels_[tmp_counter],fmt='%d')
    numpy.savetxt("%s/%s_distances_.txt"%(out_kcenters_distances_dir, ifile[:-4]),kcenters.distances_[tmp_counter],fmt='%18.5f')
    tmp_counter+=1

print kcenters.cluster_centers_.save_xyz('center_.xyz')
