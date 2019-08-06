
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
#display images generated in line

'''
#Description about MSM

Markov State Model (MSM) is a powerful method to investigate the conformational dynamics of proteins. 
MSM achieves the analysis by building a master equation from the Molecular Dynamics (MD) dataset.
Most of the kinetic and thermodynamic properties can be computed from the master equation

The general protocol of MSM includes: 
1.Clustering the conformational space into states. This is usually achieved by the splitting-and-lumping algorithm.
   1.0. We perform tICA to select kinetically slow subspace for the clustering stage
   1.1. We  split the conformational space into hundreds of microstates based on the geometric proximity.
        KMeans or kcenters is commonly used in the splitting stage.
   1.2. We then lump the microstates into macrostates based on the kinetic proximity, Spectral-based clustering is commonly applied.

2. Build the master equation (P(\tau)=T(\tau) P(0), T(\tau) is the transition probability matrix)
   2.0. Select a proper lag time so that the state-model will be Markovian (lose memory). The lag time is often chosen based on the implied timescale plot
   2.1. Calculate the transition probability matrix (TPM). TPM is the key to MSM.
   (We can build MSM on microstate model or macrostate model, as long as the model becomes Markovian)
3. Extract kinetic properties from the MSM, for example: Mean first passage time (MFPT) and transition pathways (obtained by transition path theory (TPT))

References:[1]Constructing Markov State Models to elucidate the functional conformational changes of complex biomolecules
           [2]MSMBuilder: Statistical Models for Biomolecular Dynamics
           [3]PyEMMA 2: A Software Package for Estimation, Validation, and Analysis of Markov Models.

This notebook contains the whole pipline to perform MSM analysis. Specifically, we run tICA, splitting-and-lumping, and TPT.

This script is written by Wei Wang (wwangat@gmail.com)

Dependency of this script: msmbuilder (>3.6.0 version), pyemma 2

'''


# In[ ]:


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
import pyemma.msm as pyemma_msm

import matplotlib
#matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import mdtraj as md
from msmbuilder.io.sampling import sample_dimension
import random

from pyemma.plots import plot_free_energy
import msmtools.analysis


# In[ ]:


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


##plots:
#draw_tica_projection: after tica, we draw the MD data projected onto major tICs
#plot_states_on_tic_space: the datapoints in same color belong to the same states
#plot_impliedtimescale: ITS shows how many basins

def draw_tica_histogram_core(x, y, tIC_a, tIC_b):
#    plt.axes(axisbg='w')
    plt.grid(False)
    plt.hist2d(x, y, bins=100, cmap='hot_r', norm=LogNorm())
    plt.xlabel('tIC %s'%(str(tIC_a)))
    plt.ylabel('tIC %s'%(str(tIC_b)))
    plt.title('tICA Heatmap (log color scale)')

def draw_tica_projection(resultdir, tica_trajs, opath, tIC_a, tIC_b):
    transformed=np.concatenate(tica_trajs)
    x = transformed[:, int(tIC_a)-1]
    y = transformed[:, int(tIC_b)-1]
    draw_tica_histogram_core(x, y,tIC_a, tIC_b)
    plt.colorbar()
    plt.savefig(resultdir + '/'+opath)
#    plt.close()
    
def plot_states_on_tic_space(resultdir, opath, tica_trajs, msm_traj, tIC_a, tIC_b):
    transformed_tic=np.concatenate(tica_trajs)
    transformed_assignments=np.concatenate(msm_traj)
    for j in range(min(transformed_assignments), max(transformed_assignments)+1):
        index = np.where(transformed_assignments == j)[0]
        plt.scatter(transformed_tic[index[::20], tIC_a-1], transformed_tic[index[::20], tIC_b-1],s=2)
    if(max(transformed_assignments)-min(transformed_assignments)<=6):
        plt.legend(['S'+str(i) for i in range(min(transformed_assignments), max(transformed_assignments)+1)])
    plt.xlabel('tIC %s'%(str(tIC_a)))
    plt.ylabel('tIC %s'%(str(tIC_b)))
    plt.title('tICA projections')
    plt.savefig(resultdir + '/'+opath)

def plot_Ramachandran(resultdir, opath, phi_psi, msm_traj):
    transformed_phi_psi=np.concatenate(phi_psi)
    transformed_assignments=np.concatenate(msm_traj)
    for j in range(min(transformed_assignments), max(transformed_assignments)+1):
        index = np.where(transformed_assignments == j)[0]
        plt.scatter(transformed_phi_psi[index[::20], 0]*180/np.pi, transformed_phi_psi[index[::20], 1]*180/np.pi,s=2)
    if(max(transformed_assignments)-min(transformed_assignments)<=6):
        plt.legend(['S'+str(i) for i in range(min(transformed_assignments), max(transformed_assignments)+1)])
    plt.xlabel('phi')
    plt.ylabel('psi')
    plt.xlim([-180, 180])
    plt.ylim([-180, 180])
    plt.title('Ramachandran plots')
    plt.savefig(resultdir + '/'+opath)
    
def plot_impliedtimescale(resultdir, opath, lag_times, msm_timescales, axis_unit):
    for i in range(len(msm_timescales[0])):
        plt.plot(lag_times, msm_timescales[:, i], 'o-')
    plt.title('Discrete-time MSM Relaxation Timescales')
    plt.xlabel(axis_unit)
    plt.ylabel(axis_unit)
    plt.semilogy()
    plt.savefig(resultdir + '/'+ opath)

##sampling the conformations

#sampling_along_tIC: intepret the physical meaning of tIC by sampling conforamtions along the dominant tICs
#sampling_along_msm_eigenmode: intepret the physical meaning of the slowest transition by sampling conformation along the msm eigenvectors

def sampling_along_tIC(resultdir, opath, tica_trajs, xtc_traj_folder, traj_list_array, pdb_name, tIC_a):
    transformed=np.concatenate(tica_trajs)
    draw_tica_histogram_core(transformed[:, 0], transformed[:, 1], '1', '2')
    tica_trajs = {i:tica_trajs[i] for i in range(len(tica_trajs))}    #tica_trajs is now a dictionary
    inds = sample_dimension(tica_trajs, dimension=tIC_a-1, n_frames=200, scheme='random') #sample 200 conformations
    #make trajectory
    traj = md.join(md.load_frame(xtc_traj_folder+traj_list_array[i], index=frame_i, top=xtc_traj_folder+pdb_name)
                  for i, frame_i in inds)
    #save the trajectory
    traj.save("%s/tica-dimension-tIC%s.xtc"%(resultdir, tIC_a-1))
    #show the samples on tICA projections
    samples_coord = []
    for i, frame_i in inds:
        samples_coord.append([tica_trajs[i][frame_i][0], tica_trajs[i][frame_i][1]])
    samples_coord = np.array(samples_coord)
    print(samples_coord.shape)
    plt.plot(samples_coord[:, 0], samples_coord[:, 1], 'o-')
    plt.legend('sample')
    plt.savefig(resultdir + '/' + opath)

def sampling_along_msm_eigenmode(resultdir, msm_eigen_trajs, assignment_list, xtc_traj_folder, traj_list_array, pdb_name, mode_num):
    #mode_num=1: the slowest dynamic mode (2nd eigenvector of TPM)
    info_list=[] #elements in info_list: [traj_id, frame_id, class_id, eigenvector_value]
    for j in range(len(assignment_list)):
        for k in range(len(assignment_list[j])):
            info_list.append([j, k, assignment_list[j][k], msm_eigen_trajs[j][k][mode_num-1]]) #index starts from 0
    #randomly sample 200 conformations
    samples=random.sample(info_list, 200) #can play with the number 200
    sorted_samples=sorted(samples, key=lambda x: x[3])
    #get the xtc file containing all the sampled conformations
    traj = md.join(md.load_frame(xtc_traj_folder+traj_list_array[i], index=frame_i, top=xtc_traj_folder+pdb_name)
                  for i, frame_i, class_id, eigenvector_value in sorted_samples)        
    traj.save("%s/msm-%s-dynamic-mode.xtc"%(resultdir, mode_num))

def sampling_representative_structures_for_MSM(resultdir, msm_traj, xtc_traj_folder, traj_list_array, pdb_name):
    transformed_assignments=np.concatenate(msm_traj)
    info_list=[] #elements in info_list: [traj_id, frame_id, class_id, eigenvector_value]
    for j in range(len(msm_traj)):
        for k in range(len(msm_traj[j])):
            info_list.append([j, k, msm_traj[j][k]])
    info_list = np.array(info_list)
    for j in range(min(transformed_assignments), max(transformed_assignments)+1):
        samples = random.sample(list(info_list[np.where(info_list[:, 2] == j)[0]]), 50)
        traj = md.join(md.load_frame(xtc_traj_folder+traj_list_array[i], index=frame_i, top=xtc_traj_folder+pdb_name)
                  for i, frame_i, class_id in samples)
        traj.save("%s/representative_structure_in_state%d.xtc"%(resultdir, j))

#run TPT to get the kinetic properties, based on the macrostate MSM, analytically
def evaluate_dominant_paths(TPM, lagtime, source_state, sink_state, time_unit):
    M = pyemma_msm.markov_model(TPM)
    tpt = pyemma_msm.tpt(M, [source_state], [sink_state])
    (paths,pathfluxes) = tpt.pathways()
    cumflux = 0
    print("Dominant pathways from state %d to state %d:"%(source_state, sink_state))
    print("path\t percentage")
    for i in range(len(paths)):
        print(paths[i], '\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux))
    print('MFPT from state %d to state %d  = %f %s'% (source_state, sink_state, M.mfpt(source_state, sink_state)*lagtime, time_unit))
#msmbuilder's mpft



# In[ ]:


#Evaluate the MFPT and TPT based on the microstate model and the micro-to-macro mapping relationship, analytically
def calculate_macro_TPT_and_MFPT_basedon_micro_MSM(microTPM, mapping, source_state, sink_state, lagtime, time_unit, verbose):
    P = pyemma_msm.markov_model(microTPM) #TPM is row-normalized
    A = np.where(mapping==source_state)[0]
    B = np.where(mapping==sink_state)[0]
    print('MFPT from state %d to state %d = %f %s'%(source_state, sink_state, msmtools.analysis.mfpt(microTPM, B, A)*lagtime, time_unit))
    #get TPT,  #bug fixed on June 4, 2019
    tpt = pyemma_msm.tpt(P, A, B)
    (paths,pathfluxes) = tpt.pathways()
    cumflux = 0
    print("summarizing the pathways in macrostate form")
    temp_file=open('tpt_temp.log', 'w')
    for i in range(len(paths)):
        cumflux += pathfluxes[i]
        temp_file.write("%16.9f\t%16.9f\t%16.9f\t%s\n"%(pathfluxes[i], 100.0*pathfluxes[i]/tpt.total_flux, 100.0*cumflux/tpt.total_flux, mapping[paths[i]]))
    temp_file.close()
    os.system('bash lump.sh') #may do this in bash
    for line in open('tpt_pathlump.log'):
        print(line.strip())
    print('############################################################')
    cumflux = 0
    if(verbose == 'on'):
        print("more details about the paths in the microstate form")
        print("Path flux\t\t%path\t%of total\tpaths")
        for i in range(len(paths)):
            cumflux += pathfluxes[i]
            print(pathfluxes[i],'\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux),'%\t','%3.1f'%(100.0*cumflux/tpt.total_flux),'%\t',paths[i])


# In[ ]:


#specify the input files:
#1. folders to store the pdb file and the trajectories (in xtc or dcd format)
#2. the atom index or index pair for the pairwise distances used in tICA (index starts from 0 in python) 
#3. the name list for the trajectories

#In alanine dipeptide, the MD data is saved every 1 ps. Run at 300 K using Langevin

trajectory_dir='trajs/'  #folder to store the pdb and trajectories
pdb_name='ala2.pdb' #name for the topology file
trajname_list='trajname.list' #the name list for the trajectories in trajectory_dir
atom_pair_list='pairdist.list' #the pairwise distances list used in tICA
resultdir = './results'


# In[ ]:


if not os.path.exists(resultdir):
    os.makedirs(resultdir)
    
#load all the input files
print('now we are loading all the input files used in the tICA-MSM analysis')
atom_pairs = np.loadtxt(atom_pair_list, dtype=int) #import the pairwise distance index file as integer type
traj_list_array=[]
for line in open(trajname_list):
    traj_list_array.append(line.strip())


# In[ ]:


#step 1.0: tICA
#Select kinetic slow variables via tICA (time-lagged independent component analysis)
#tICA finds the linear combination of the input features that maximizing the normalized time-lagged correlation matrix
#In this example, we use pairwise distance of all heavy atoms as the input features for tICA.

#input: trajectories, output: tICA projections
#prepare data for tICA
featurizer = AtomPairsFeaturizer(pair_indices=atom_pairs) #In this example, we use pairwise distances
pairdist4tica = featurizing_the_conformations(featurizer, trajectory_dir, traj_list_array, pdb_name)
print("now we have prepared the data for tICA: the pairwise distances for all frames in all trajectories")

#run tICA
tica_model=tICA(lag_time=10,n_components=2) #tica lagged should be pre-specified, you can play with this number!
tica_trajs=tica_model.fit_transform(pairdist4tica)  #projected the MD data onto tica coordinates
#print("output of tica:", tica_trajs)
#plot the tica projections
draw_tica_projection(resultdir, tica_trajs,'tica_12.png', 1, 2)

#sample conformations along tIC1
print('now we are sampling representative conformations along tIC1')
plt.figure()
sampling_along_tIC(resultdir, 'samples_tic1.png', tica_trajs, trajectory_dir, traj_list_array, pdb_name, 1)
print("You can use vmd to visualize the tica-dimension-tIC1.xtc file")


# In[ ]:


#step 1.1: split the conformations into hundreds of microstates
#perform kCenters on the tIC subspace
#input:tICA projections, output:assignments indicating which microstate each conformation is assigned to
nMicro=100 #specified a priori
kcenters=KCenters(n_clusters=nMicro, metric='euclidean', random_state=0)
microstate_sequences = kcenters.fit(tica_trajs)
print("output of msm:", microstate_sequences.labels_)

plt.figure()
plot_states_on_tic_space(resultdir, 'micorstate.png', tica_trajs, microstate_sequences.labels_, 1, 2)


# In[ ]:


#plot the microstate implied timescale, which will show how many macrostates we need
plt.figure()
lag_times=range(2,50,2)
msm_timescales = implied_timescales(microstate_sequences.labels_, lag_times, n_timescales=10,msm=MarkovStateModel(reversible_type='transpose', ergodic_cutoff='off'))
plot_impliedtimescale(resultdir, 'microstate_its.png', lag_times, msm_timescales, 'ps')


# In[ ]:


#the first dynamic eigenvector is associated to the slowest transitions in the dataset
#we can understand the physical meaning of the first eigenmode through sampling the conformations
msm = MarkovStateModel(lag_time=4, reversible_type='transpose', n_timescales=3, ergodic_cutoff='off') #lag time should be chosen such that the model becomes Markovian
#n_timescale specify the number of dynamic mode (the 1st one is the 2nd eigenvector of TPM) that outputs
msm.fit(microstate_sequences.labels_)
msm_eigen_trajs = msm.eigtransform(microstate_sequences.labels_)
sampling_along_msm_eigenmode(resultdir, msm_eigen_trajs, microstate_sequences.labels_, trajectory_dir, traj_list_array, pdb_name, 1) #1:the slowest transition,play with the eigenmodes
print("the timescale associated with the slowest collective motion in the system is: %d ps"%(msm.timescales_[0]))
print("using vmd to open msm-1-dynamic-mode.xtc, to intepret the slowest dynamic mode of the system")


# In[ ]:


#judging from the above implied timescale, we need to lump the microstates into 3 macrostates
#lump kinetically close microstates into a few macrostates using PCCA+ algorithm, which will facilitate the visualization and intereptation of kinetics of the system
pcca = PCCAPlus.from_msm(msm, n_macrostates=3)
macro_trajs = pcca.transform(microstate_sequences.labels_)
#show the macrostates onto tICA space
plt.figure()
plot_states_on_tic_space(resultdir, 'macrostate.png', tica_trajs, macro_trajs, 1, 2)

#for alanine dipeptide, we can also validate the lumping using Ramanchandran plots
phi_psi_featurizer = DihedralFeaturizer(types=['phi', 'psi'], sincos=False) #if "True", then output sin, cos of the angles
phi_psi = featurizing_the_conformations(phi_psi_featurizer, trajectory_dir, traj_list_array, pdb_name)
plt.figure()
plot_Ramachandran(resultdir, 'macrostate_Ramachandran.png', phi_psi, macro_trajs)


# In[ ]:


#sample representative structures for each macrostate, used for visualization

sampling_representative_structures_for_MSM(resultdir, macro_trajs, trajectory_dir, traj_list_array, pdb_name)
print("you can now visualize the representative structures using vmd")


# In[ ]:


#plot the macrostate implied timescale, we can then select proper lag time to build the macrostate MSM
plt.figure()
lag_times=range(2,50,2)
msm_timescales = implied_timescales(macro_trajs, lag_times, n_timescales=3,msm=MarkovStateModel(reversible_type='transpose', ergodic_cutoff='off'))
plot_impliedtimescale(resultdir, 'macrostate_its.png', lag_times, msm_timescales, 'ps')


# In[ ]:


#evaluate the dynamics based on the microstate MSM
lagtime = 10 #microstate Markovian lag time
microstate_msm = MarkovStateModel(lag_time=lagtime, reversible_type='transpose', n_timescales=3, ergodic_cutoff='off') #lag time should be chosen such that the model becomes Markovian
#n_timescale specify the number of dynamic mode (the 1st one is the 2nd eigenvector of TPM) that outputs
microstate_msm.fit(microstate_sequences.labels_)
source_state = 1 #specify the source state you want to investigate
sink_state = 2 #specify the sink macrostate you want to investigate
calculate_macro_TPT_and_MFPT_basedon_micro_MSM(microstate_msm.transmat_, pcca.microstate_mapping_, source_state, sink_state, lagtime, '0.1 ps', 'off')


# In[ ]:


#Build the macrostate MSM and extract kinetics & thermodynamics based on the MSM
lag_time = 12 #the lag time where implied timescale reaches plateau 
msm = MarkovStateModel(lag_time, reversible_type='transpose', n_timescales=2) #lag time should be chosen such that the model becomes Markovian
#n_timescale specify the number of dynamic mode (the 1st one is the 2nd eigenvector of TPM) that outputs
msm.fit(macro_trajs)

#get the first dynamic eigenmode
print("the first dynamic eigenmode:", msm.right_eigenvectors_[:, 1])
print("the slowest transition happens between the macrostates with +ev & -ev eigenvector components")

#get the stationary population
print("Stationary population for the macrostates:", msm.populations_)

#get the transition pathways and mean first passage time from one state to another, using transition path theory (tpt)
source_state = 1
sink_state = 2
evaluate_dominant_paths(msm.transmat_, lag_time, source_state, sink_state, 'ps') #play with source state and sink state


# In[ ]:


print("end of the example of alanine dipeptide")

