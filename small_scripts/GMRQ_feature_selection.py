'''

Written by Wang Wei on Mar 17, 2019.
Scan the parameters to get the best description of the slowest transition modes of the system
Evaluating the tICA-MSM models based on GMRQ, in the framework of cross-valiation
Parameters to evaluate: pairwise distance set, tICA correlation time, # of tICs used in clustering, # of microstates

dependency of the codes:python3, scikit-learn > v0.20, msmbuilder
please contact Wei Wang (wwangat@gmail.com) once you have encountered some problems

'''
#load all the modules needed in the MSM analysis of the alanine dipeptide system
import os
import sys
import warnings
warnings.filterwarnings('ignore')
import numpy as np

from msmbuilder.featurizer import AtomPairsFeaturizer, DihedralFeaturizer
from msmbuilder.cluster import KCenters
from msmbuilder.decomposition import tICA
from msmbuilder.dataset import dataset
from msmbuilder.msm import implied_timescales, MarkovStateModel
from msmbuilder.lumping import PCCAPlus
#import pyemma.msm as pyemma_msm

import matplotlib
#matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import mdtraj as md
from msmbuilder.io.sampling import sample_dimension
import random

from sklearn.model_selection import KFold


#functions that will be used in the analysis
#featurizing_the_conformations: transform the xyz coordinates into features, the output will be used for tICA
def featurizing_the_conformations(featurizer, xtc_traj_folder, traj_list_array, pdb_name):
    #we usually use dihedral, rmsd or pairwise distance featurizers in MSM
    output_features=[]
    for trajfile in traj_list_array:
        xyz=dataset(trajectory_dir+trajfile, topology=trajectory_dir+pdb_name) #xyz is the coordinates for the frames in the trajectory
        temp=featurizer.fit_transform(xyz) #we transform the coordinates into pairwise distaneces
        print("loaded %s into (%d,%d) dimensional file"%(trajfile, len(temp[0]), len(temp[0][0])))
        output_features.append(temp[0])
    return output_features

def draw_tica_projection_cross_validation(resultdir, opath, train_data_projection, test_data_projection, tIC_a, tIC_b):
    transformed1 = np.concatenate(train_data_projection)
    plt.scatter(transformed1[::20, tIC_a-1], transformed1[::20, tIC_b-1],s=2, color='black')
    transformed2 = np.concatenate(test_data_projection)
    plt.scatter(transformed2[::20, tIC_a-1], transformed2[::20, tIC_b-1],s=2, color='red')
    plt.xlabel('tIC %s'%(str(tIC_a)))
    plt.ylabel('tIC %s'%(str(tIC_b)))
    plt.legend(['train', 'test'])
    plt.title('tICA projections')
    plt.savefig(resultdir+'/'+opath)

def fit_predict_tica(new_frame_xyz, tica_model): #having tested to be correct
    tica_projection=[]
    for j in range(len(new_frame_xyz)):
        tica_projection.append(np.dot(new_frame_xyz[j]-tica_model.means_.T, tica_model.components_[:, :].T))
    return tica_projection

trajectory_dir='trajs/'  #folder to store the pdb and trajectories
pdb_name='ala2.pdb' #name for the topology file
trajname_list='trajname.list' #the name list for the trajectories in trajectory_dir
pairwise_distances_files_list='pairwise_distances_files_list.list'

resultdir = './results' #the output folder

if not os.path.exists(resultdir):
    os.makedirs(resultdir)

traj_list_array=[]
for line in open(trajname_list):
    traj_list_array.append(line.strip())
#shuffle the trajectories to make the cross-validation stage more unbiased
np.random.shuffle(traj_list_array)

print("the shutffled trajectory list reading into memory are as follows:")
print(traj_list_array)


#Parameter range we need to tune in the process
n_tics_range=range(2, 3, 1)
n_Micro_range=range(100, 200, 100)
tica_correlation_time_range=range(2, 3, 1)

n_splits = 5
temp_num=0
for features_file in open(pairwise_distances_files_list):
    temp_num+=1
    print('----------------------------------------------------------------------------------------')
    print("now we are handling the feature file:", features_file.strip())
    atom_pairs = np.loadtxt(features_file.strip(), dtype='int')
    print("the features we are handling are:\n", atom_pairs)

    sub_resultdir = resultdir+'/feature_list'+str(temp_num)+'/'
    if not os.path.exists(sub_resultdir):
        os.makedirs(sub_resultdir)

    featurizer = AtomPairsFeaturizer(pair_indices=atom_pairs)
    data = featurizing_the_conformations(featurizer, trajectory_dir, traj_list_array, pdb_name)

    cv = KFold(n_splits=n_splits, shuffle=False) #5-fold cross validation, exclusive
    fold = 0
    for (train_index, test_index) in cv.split(traj_list_array):
        fold += 1
        print("now we are handling fold %d"%(fold))
        print("training data:", [traj_list_array[i] for i in train_index])
        print("testing data", [traj_list_array[i] for i in test_index])

        train_data = [data[i] for i in train_index]
        test_data = [data[i] for i in test_index]

        for tica_correlation_time in tica_correlation_time_range:
            tICA_model = tICA(lag_time = tica_correlation_time, n_components = 10)
            tICA_model.fit(train_data)
            train_data_projection = fit_predict_tica(train_data, tICA_model)
            test_data_projection = fit_predict_tica(test_data, tICA_model) #just deduct the same mean as they are of the same ensemble
            plt.figure()
            draw_tica_projection_cross_validation(sub_resultdir, 'Fold_%d_tica_lagtime_%d_train_data_proj_tIC12.png'%(fold, tica_correlation_time), train_data_projection, test_data_projection, 1, 2)
            plt.figure()
            draw_tica_projection_cross_validation(sub_resultdir, 'Fold_%d_tica_lagtime_%d_train_data_proj_tIC13.png'%(fold, tica_correlation_time), train_data_projection, test_data_projection, 1, 3)

            for n_tics in n_tics_range:
                for n_Micro in n_Micro_range:
                    print("parameters: fold-", fold, ',tica_lagtime-', tica_correlation_time, ',n_tics-', n_tics, ',n_Micro-', n_Micro)
                    kcenters = KCenters(n_clusters=n_Micro, metric='euclidean', random_state=0)
                    kcenters.fit(train_data_projection)
                    train_data_sequence = kcenters.predict(train_data_projection)
                    test_data_sequence = kcenters.predict(test_data_projection)
                    msm = MarkovStateModel(n_timescales=3, lag_time=100, reversible_type='transpose', verbose=False, sliding_window=True, ergodic_cutoff='on') #the parameters may change
                    msm.fit(train_data_sequence)
                    train_score = msm.score(train_data_sequence)
                    test_score = msm.score(test_data_sequence)
                    f1 = open(sub_resultdir+'/Fold_%d_tica_lagtime_%d_ntics_%d_nMicro_%d_gmrq.summary'%(fold, tica_correlation_time, n_tics, n_Micro), 'w')
                    f1.write('train_score:%f'%(train_score))
                    f1.write('\n')
                    f1.write('test_score:%f'%(test_score))
                    f1.write('\n')
                    f1.close()
                    print('computing implied timescale for training data')   #the x-range to plot implied timescale should also change
                    train_msm_timescales = implied_timescales(train_data_sequence, range(1, 20, 1), n_timescales=10, msm=MarkovStateModel(reversible_type='transpose', ergodic_cutoff='on'))
                    np.savetxt(sub_resultdir+'/Fold_%d_tica_lagtime_%d_ntics_%d_nMicro_%d_traindata_its.dat'%(fold, tica_correlation_time, n_tics, n_Micro), train_msm_timescales)
                    print('computing implied timescale for testing data')
                    test_msm_timescales = implied_timescales(test_data_sequence, range(1, 20, 1), n_timescales=10, msm=MarkovStateModel(reversible_type='transpose', ergodic_cutoff='on'))
                    np.savetxt(sub_resultdir+'/Fold_%d_tica_lagtime_%d_ntics_%d_nMicro_%d_testdata_its.dat'%(fold, tica_correlation_time, n_tics, n_Micro), train_msm_timescales)
    print('We have finished the analysis for this feature set, now collecting the cross-validation data')
    train_score_collection = []
    test_score_collection = []
    for tica_correlation_time in tica_correlation_time_range:
        for n_tics in n_tics_range:
            for n_Micro in n_Micro_range:
                for fold in range(1, n_splits+1, 1):
                    temp_list=[]
                    filename=sub_resultdir+'/Fold_%d_tica_lagtime_%d_ntics_%d_nMicro_%d_gmrq.summary'%(fold, tica_correlation_time, n_tics, n_Micro)
                    for line in open(filename):
                        temp_list.append(line.strip().split(':')[1])
                    train_score_collection.append(float(temp_list[0]))
                    test_score_collection.append(float(temp_list[1]))
                #collect and compute the statistics
                f1 = open(sub_resultdir+'/tica_lagtime_%d_ntics_%d_nMicro_%d_gmrq_statistics.stat'%(tica_correlation_time, n_tics, n_Micro), 'w')
                f1.write('mean_train_score:%f, std_train_score:%f, mean_test_score:%f, std_test_score:%f\n'%(np.mean(train_score_collection), np.std(train_score_collection), np.mean(test_score_collection), np.std(test_score_collection)))
                f1.close()
