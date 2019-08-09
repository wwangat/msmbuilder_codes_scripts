
# coding: utf-8

# In[21]:

#get_ipython().magic(u'pylab inline')


# In[2]:

import numpy as numpy
from numpy import *
import scipy as scipy
from scipy import linalg
import mdtraj as mdtraj
import optparse
import os
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import KMeans
from msmbuilder.msm import MarkovStateModel
from msmbuilder.featurizer import AtomPairsFeaturizer
from sklearn.pipeline import Pipeline
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from msmbuilder.cluster import KCenters
from msmbuilder.utils.progressbar import ProgressBar, Percentage, Bar, ETA
from msmbuilder.cluster import MultiSequenceClusterMixin
from msmbuilder.utils import check_iter_of_sequences, array2d
from msmbuilder.msm import implied_timescales

coord = []
atomName = []

projMode = 2  #number of tics you want to interpret

for line in file('clean.pdb'): #clean pdb
    atomName.append(line[13])
    line=line.strip().split()
    coord.append([float(line[5]),float(line[6]),float(line[7])])

coord = numpy.array(coord)

perturb = []

for line in file('tica_components'):#tica components given by msmbuilder, first two columns are the atom id for each pairwise distance 
    line=line.strip().split()
    perturb.append([int(line[0]),int(line[1]), float(line[projMode+1])])

#print len(perturb)

atoms = []
for line in file('atom_indices'):#atoms used in the pairwise distances, index from 0
    atoms.append(int(line.strip()))

#atoms = numpy.concatenate([range(0,26), range(44,70)],0)
map = {}
for i in range(len(atoms)):
    map[atoms[i]] = i

refPos = numpy.array(coord[atoms])
#print refPos

distMtx = numpy.zeros((len(atoms),len(atoms)), float)
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        distMtx[i][j] = numpy.sqrt((coord[atoms[i]][0]-coord[atoms[j]][0])**2+(coord[atoms[i]][1]-coord[atoms[j]][1])**2+(coord[atoms[i]][2]-coord[atoms[j]][2])**2)
        distMtx[j][i] = distMtx[i][j]

#print distMtx


#distMtx = [[0,3,4],[3,0,5],[4,5,0]]
#distMtx = numpy.array(distMtx)

from sklearn import manifold
seed = numpy.random.RandomState(seed=3)

mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1) #adjust hte parameters

steps = 24 #number of snapshots

for i in range(steps):
    distMtxTMP = numpy.copy(distMtx)
    delta = 20*sin(i/float(steps)*2*numpy.pi) #modify 20 to other numbers, here 20 is kind of velocity
    for j in range(len(perturb)):
        a = map[perturb[j][0]]
        b = map[perturb[j][1]]
        distMtxTMP[a][b] += delta * perturb[j][2]
        distMtxTMP[b][a] += delta * perturb[j][2]

    pos = mds.fit(distMtxTMP).embedding_

    U, S, VT = numpy.linalg.svd(numpy.dot(refPos.T,pos))
    pos = numpy.dot(pos, numpy.dot(VT.T,U.T))

    print len(pos)
    print

    for j in range(len(pos)):
#        print atoms[j]
#        print map[atoms[j]]
        print '%s %8.3f %8.3f %8.3f' %(atomName[atoms[j]], pos[j][0], pos[j][1], pos[j][2])

