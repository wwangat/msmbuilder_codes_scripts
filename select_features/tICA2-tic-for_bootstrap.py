"""
Writer: Wei WANG (wwangat@gmail.com)
Function: get the tICA eigenvalues and eigenvectors for the boostrapping step, to select significant atom pairs

"""
from sklearn.pipeline import Pipeline
from itertools import combinations
from msmbuilder.cluster import KCenters
from msmbuilder.featurizer import AtomPairsFeaturizer
import numpy
from msmbuilder.dataset import dataset
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import os
from msmbuilder.decomposition import tICA
from msmbuilder.msm import implied_timescales
from msmbuilder.msm import MarkovStateModel
import sys
from matplotlib.colors import LogNorm

def fit_and_plot(transformed,opath, tIC_a, tIC_b):
    transformed=numpy.concatenate(transformed)
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


pairdist = numpy.loadtxt('../pairwise_dist5',dtype=int)  #read in pairwise distance directly, inputs

#atomindex_file="../AtomIndices.dat"   #input atom indices file, inputs
#pairdist=create_pairwise_index(atomindex_file)
print pairdist
##format of this atom index file:
##number_of_atoms_included ..,..,..,(atom indices starting from 0)
##example:
##35 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34

#xtc_file_dir = '../../xtc-align/';   #inputs, where you put the xtc file
#xyz = dataset('%s/*.xtc'%(xtc_file_dir), topology = '../1-noWAT-1.pdb')  #inputs, where you put your pdb file
#os.system("for j in %s/*.xtc;do echo `basename $j .xtc`;done>trajlist"%(xtc_file_dir))


#traj_list_array=[]
#for line in open("trajlist"):
#    traj_list_array.append(line.strip())
#print traj_list_array

#featurizer = AtomPairsFeaturizer(pair_indices=pairdist)

print "begin to fit trajs to calculate the pairwise distances"

#pairdist4tica = xyz.fit_transform_with(featurizer, 'atompairsfeaturizer/', fmt='dir-npy')

print "begin to calculate the tica implied timescale, tica ITS can guide us how many tics to use"

#print out tica eigenvalues and eigenvectors
n_components = 50
lag_time = 20

pairdist4tica = dataset('atompairsfeaturizer/', mode='r', fmt='dir-npy')  #the atom features is bootstrapped

tica=tICA(lag_time=lag_time, n_components=n_components)
print "begin to calculate tica eiganvalues and eigenvectors"
tica.fit(pairdist4tica)


f=open('tica_eigenvalues_%d'%(lag_time),'w')
for j in range(n_components):
    f.write("%f\n  "%(tica.eigenvalues_[j]))
f.close()

numpy.savetxt('tica_eigenvectors_%d'%(lag_time), tica.eigenvectors_[:, :])

sys.exit()
