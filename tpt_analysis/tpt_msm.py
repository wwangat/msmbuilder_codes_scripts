# In[1]:


import numpy as numpy
from numpy import *
import optparse 
import os
from optparse import OptionParser
from msmbuilder.msm import MarkovStateModel
import pyemma.msm as msm
# In[3]:

def create_pairwise_index(atomindex_file):
    tmp_array=[]
    counter=0
    for line in file(atomindex_file):
        line=line.strip().split()
        for i in range(len(line)):
            if ((counter==0)&(i==0)):  # the first number in the index file of msmbuilder1 is the total number
                continue
            else:
                tmp_array.append(int(line[i])) # starting from 0 and shift the atom index in pdb by 1; 
        counter+=1
    output_array=numpy.unique(tmp_array)
    return list(combinations(output_array, 2))

def read_traj(path="", topo=None):
    tmp_traj=mdtraj.load(path, top=topo,discard_overlapping_frames=True)
    return tmp_traj

def read_trajs_from_msm_workflow(rootDIR, traj_files_array=[],topo=None):
    #Read the trajs based on the trajlist 
    trajs_array=[]
    print "Beware, I'll try to read %d trajlists. I am memory innefficient so don't be surprised if your computer breaks"%len(traj_files_array) 
    counter=0
    for ifile in traj_files_array:
        tmp_traj_array=[]
        counter+=1
        print "Reading traj #%03d from trajfile %s"%(counter, ifile)
        for ipath in file("%s/trajectories/%s"%(rootDIR, ifile)):
            ipath=ipath.strip()
            if (len(ipath.split()) > 1 ):
                print "Something is wrong with this name (don't use spaces), I'll beter stop: "
            ipath_ext=rootDIR+"/"+ipath+".xtc"
            print ipath_ext
            print os.path.exists(ipath_ext)
            tmp_traj_array.append(read_traj(ipath_ext,topo))
        trajs_array.append(tmp_traj_array)
    return trajs_array

def fit_and_plot(pipeline, trajectories,opath, tIC_a, tIC_b):
    transformed = pipeline.fit_transform(trajectories)
    transformed = numpy.concatenate(transformed)
    
    print('Eiegenvalue sum', pipeline.named_steps['tica'].eigenvalues_.sum())

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

parser = OptionParser()
parser.add_option('-l',"--trajlist_file",default = "trajlist", help="the trajlist indicating all the microstate assignment")
parser.add_option('-m',"--mapping_file",default = "mapping.txt", help="the microstate to macrostate mapping relationship, start from 0")
#parser.add_option('-o',"--resultdir", default='.', help="where to put the results")
parser.add_option('-t','--lag_time_msm',help='lagtime is needed',default=1,type='int')
parser.add_option('-i','--source',help='indicate the source, index start from 0',type='int')
parser.add_option('-e','--sink',help='indicate the sink, index start from 0',type='int')

(options, args) = parser.parse_args()
if options.trajlist_file == "noInput" or options.mapping_file == "noInput" or options.lag_time_msm == "noInput":
    print "Please input all the necessary file according to python tpt_msm.py -h"
if options.source == "noInput" or options.sink == "noInput":
    print "Both source and sink state should be indicated"

#reading the mapping relationship
mapping = numpy.loadtxt(options.mapping_file, dtype='int')
sink=options.sink
source=options.source
#Read the trajlist
traj_files_array=[]
for iline in file(options.trajlist_file):
    iline=iline.strip().rstrip('\n').split()
    if (len(iline) > 1 ):
        print "Something is wrong with this line, I'll beter stop: "
        print iline
    traj_files_array.append(iline[0])
print "I found %d trajfiles"% len(traj_files_array)
assignment_array=[]
for i in xrange(len(traj_files_array)):
    assignment_array.append(numpy.loadtxt("%s"%( traj_files_array[i]),dtype='int'))
#####having read microstate assignment
#lag_time_msm=200 ## lag time to build MSM, in the unit of frame interval
## set parameters to build MSM

msm_model = MarkovStateModel(lag_time = options.lag_time_msm, reversible_type='mle', sliding_window=1, ergodic_cutoff='on', prior_counts=0, verbose=True)
msm_model.fit(assignment_array)

print msm_model.transmat_
TPM = msm.markov_model(msm_model.transmat_)

########################using pyemma to get transition paths
set1 = numpy.where(mapping == source)[0]
set2 = numpy.where(mapping == sink)[0]
#print set1
#print set2
tpt=msm.tpt(TPM, set1, set2)

## get tpt gross flux
#print "grpss flux matrix"
#print tpt.gross_flux
#print "forward committor matrix"
#print tpt.committor
#print "backward committor matrix"
#print tpt.backward_committor_matrix
#
## get tpt net flux
#Fp = tpt.net_flux
#print "net flux matrix"
#print Fp
#
#print "###########################################################"
#
#print 'Total TPT flux = ', tpt.total_flux
##print 'Rate from TPT flux = ', tpt.rate
#print 'Rate from TPM flux = %f microseconds' % (1.0*0.1*1e-3/tpt.rate) #since time unit is 200ps, and lagtime is 500 steps
#
##mplt.plot_flux(tpt, pos=tptpos, flux_scale=100.0/tpt.total_flux, arrow_label_format="%3.1f")
##ylabel("committor")
#
#print "#############################################################"
##print "all the pathways from pre-translocation state to backtracked state are:"
#tpt.pathways()
(paths,pathfluxes) = tpt.pathways(fraction=1.0) #maximum iteratively give 1000 different paths. You can refer to modified Dijkstra's algorithm and flux decomposition to learn more about it
cumflux = 0
print "in macrostate form"
print "Path flux\t\t%path\t%of total\tpath"
for i in range(len(paths)):
    cumflux += pathfluxes[i]
#    print pathfluxes[i],'\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux),'%\t','%3.1f'%(100.0*cumflux/tpt.total_flux),'%\t\t',paths[i]
    print pathfluxes[i],'\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux),'%\t','%3.1f'%(100.0*cumflux/tpt.total_flux),'%\t\t',mapping[paths[i]]
##print "state 0: S4-little, state 1: fray state S3, state 2: fray state S2, state 3: backtracked state, state 4: S4-middle, state 5: S4-3rd, state 6: S4-2nd, state 7: pre-translational state??"
#
#
##for i in range(len(transformed)):
##    numpy.savetxt('pcca_macroAssign%d' %(i),transformed[i], fmt='%d')
#
#
##net_flux=net_fluxes(lumping.metastable_sets[0],lumping.metastable_sets[1],msm)
## top 10 paths from macro3 to macro2
##paths(lumping.metastable_sets[0],lumping.metastable_sets[1],net_flux,num_paths=10)
