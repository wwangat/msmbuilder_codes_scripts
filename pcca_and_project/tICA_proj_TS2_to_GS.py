#project GS to TS major tICs or the other way around

from msmbuilder.cluster import KCenters
import numpy
from msmbuilder.dataset import dataset
import matplotlib
matplotlib.use('Agg')
from msmbuilder.msm import MarkovStateModel
import matplotlib.pyplot as plt
import os
from msmbuilder.decomposition import tICA
from msmbuilder.msm import implied_timescales

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

outputdir='project_GS_AND_TS_to_GS_tICs/';  #adjust variables
os.system("mkdir %s"%(outputdir));

GS_pairdist_dir='GS/atompairsfeaturizer/'    #where you put GS pairwise distance, adjust variables
TS_pairdist_dir='TS2/atompairsfeaturizer/'    #where you put TS pairwise distance, adjust variables

#read in GS and TS pairwise distance directory
os.system("for j in %s/*.npy;do echo `basename $j .npy`;done>pairdist_GS"%(GS_pairdist_dir))
os.system("for j in %s/*.npy;do echo `basename $j .npy`;done>pairdist_TS"%(TS_pairdist_dir))

GS_pair_list_array=[]
for line in open("pairdist_GS"):
    GS_pair_list_array.append(line.strip())

TS_pair_list_array=[]
for line in open("pairdist_TS"):
    TS_pair_list_array.append(line.strip())

print "You are dealing with pairwise distance with dimensions:", len(TS_pair_list_array)

n_components = 10  #adjust variable
tic_lag_time = 10   #adjust variable
nMicro=800   #adjust variable

GS_pairdist4tica = dataset(GS_pairdist_dir, mode='r', fmt='dir-npy')   #just use it to calculate major tic eigenvectors and mean of each dimension, no need to care about the order of read trajs

GS_tica=tICA(lag_time=tic_lag_time, n_components=n_components)
print "begin to calculate GS tica eiganvalues and eigenvectors under the params: tic lagtime=", tic_lag_time, "#components=", n_components
tica_sequences=GS_tica.fit_transform(GS_pairdist4tica)

f=open('%s/GS_tica_eigenvalues_lagtime%d_tIC%d'%(outputdir, tic_lag_time, n_components),'w')
for j in range(n_components):
    f.write("%f\n  "%(GS_tica.eigenvalues_[j]))
f.close()

numpy.savetxt('%s/GS_tica_eigenvectors_lagtime%d_tIC%d'%(outputdir, tic_lag_time, n_components), GS_tica.eigenvectors_[:, :])   #we thus save all the top tic eigenvectors

numpy.savetxt('%s/GS_pairwise_means'%(outputdir), GS_tica.means_)  #GS_tica.means_ is a vector

print "just for reference of GS projs,tested to be correct by ww:"
print tica_sequences

print "We are now re-read the pairwise distance and calculate the GS tic projections"

os.system("mkdir %s/GS_project_onto_GS_tics"%(outputdir))
for line in range(len(GS_pair_list_array)):
    temp = numpy.load("%s/%s.npy"%(GS_pairdist_dir, GS_pair_list_array[line]))
    #we begin to project
    numpy.savetxt("%s/GS_project_onto_GS_tics/%s_ticproj.txt"%(outputdir, GS_pair_list_array[line]), numpy.dot((temp-GS_tica.means_.T),GS_tica.eigenvectors_[:, :]))

TS2_pairdist4tica = dataset(TS_pairdist_dir, mode='r', fmt='dir-npy')   #just use it to calculate major tic eigenvectors and mean of each dimension, no need to care about the order of read trajs

TS2_tica=tICA(lag_time=tic_lag_time, n_components=n_components)

print "We are now read the pairwise distance and calculate the TS2 tic projections (must first do tica to get the tica.means because i am lazy)"

TS2_tica.fit(TS2_pairdist4tica)

numpy.savetxt('%s/TS2_pairwise_means'%(outputdir), TS2_tica.means_)  #GS_tica.means_ is a vector

os.system("mkdir %s/TS2_project_onto_GS_tics"%(outputdir))
for line in range(len(TS_pair_list_array)):
    temp = numpy.load("%s/%s.npy"%(TS_pairdist_dir, TS_pair_list_array[line]))
    #we begin to project
    numpy.savetxt("%s/TS2_project_onto_GS_tics/%s_ticproj.txt"%(outputdir, TS_pair_list_array[line]), numpy.dot((temp-TS2_tica.means_.T), GS_tica.eigenvectors_[:, :]))


#we now begin to do kcenters on reduced tic space of GS
print "begin to do kcenters on reduced tic space of GS"
#here we use the previous defined tic correlation time

num_tics_for_clustering=3   #adjust variables,get insight from tic implied timescale

#we read coordinates from outside to do kcenters in case trajectories are not of the same length

os.system("for j in %s/GS_project_onto_GS_tics/*ticproj.txt;do echo `basename $j _ticproj.txt`;done>ticproj_GS"%(outputdir))
os.system("for j in %s/TS2_project_onto_GS_tics/*ticproj.txt;do echo `basename $j _ticproj.txt`;done>ticproj_TS2"%(outputdir))

tica_sequences=[]
GS_ticproj_list_array=[]
for line in open("ticproj_GS"):
    line=line.strip()
    GS_ticproj_list_array.append(line)
    temp=numpy.loadtxt("%s/GS_project_onto_GS_tics/%s_ticproj.txt"%(outputdir, line))
    temp=temp[:, 0:num_tics_for_clustering]
    tica_sequences.append(temp)
    
TS2_ticproj_list_array=[]
tica_TS2_sequences=[]
for line in open("ticproj_TS2"):
    TS2_ticproj_list_array.append(line.strip())
    temp1=numpy.loadtxt("%s/TS2_project_onto_GS_tics/%s_ticproj.txt"%(outputdir, line.strip()))
    temp1=temp1[:, 0:num_tics_for_clustering]
    tica_TS2_sequences.append(temp1)

tmp_counter=0

kcenters = KCenters(n_clusters=nMicro)
#kcenters = KCenters(n_clusters=num_tics_for_clustering)        # Fr :)

kcenters_sequences = kcenters.fit_predict(tica_sequences)  #here it is ground state tica sequences

print "begin to plot the microstate implied timescale into the objective dir"
            #plot implied timescale

lag_times=range(10,100,10);  #adjust variables
n_timescales=5  #adjust variables

msm_timescales = implied_timescales(kcenters_sequences, lag_times, n_timescales=n_timescales, msm=MarkovStateModel(verbose=True, reversible_type='transpose'))

outfile_name="%s/GS_ITS_tic%d_lagtime%d_clustersize%d.dat"%(outputdir, num_tics_for_clustering, tic_lag_time, nMicro)
print msm_timescales
print msm_timescales.shape

for k in range(n_timescales):
    plt.plot(lag_times, msm_timescales[:, k], 'o-')
f2=open(outfile_name, 'w')
for i in range(len(lag_times)):
    f2.write("%d    "%(lag_times[i]))
    for j in range(n_timescales):
        f2.write("%f    "%(msm_timescales[i, j]))
    f2.write('\n')
f2.close()
plt.title('Discrete-time MSM Relaxation Timescales')
plt.semilogy()
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,y1,100000))  #need to change this number
outfile_name="%s/GS_ITS_tic%d_lagtime%d_clustersize%d.png"%(outputdir, num_tics_for_clustering, tic_lag_time, nMicro)
plt.savefig(outfile_name)
plt.close()

print "begin to write the GS microstate assignments into the objective dir"
tmp_counter=0
GS_assignment=outputdir+"/GS_microAssignment/"
os.system("mkdir %s"%(GS_assignment))
for ifile in range(len(GS_pair_list_array)):
    outfile_name="%s/%s_%d_%d_%d.txt"%(GS_assignment, GS_pair_list_array[ifile], num_tics_for_clustering, tic_lag_time, nMicro)
    numpy.savetxt(outfile_name,kcenters.labels_[tmp_counter],fmt='%s')
    tmp_counter+=1

print "We project TS2 onto GS cluster centers"
print "the cluster centers are:", kcenters.cluster_centers_

numpy.savetxt("%s/GS_tic%d_lagtime%d_clustersize%d_cluster_center_coords.dat" % (outputdir, num_tics_for_clustering, tic_lag_time, nMicro), kcenters.cluster_centers_)

print "predict TS2 microstate assignments based on GS cluster centers"
print tica_TS2_sequences


TS2_ticproj_list_array=[]
kcenters_TS2_sequences=[]
for line in open("ticproj_TS2"):
    TS2_ticproj_list_array.append(line.strip())
    temp1=numpy.loadtxt("%s/TS2_project_onto_GS_tics/%s_ticproj.txt"%(outputdir, line.strip()))
    temp1=temp1[:, 0:num_tics_for_clustering]
#    print 'HELLO', temp1
    temp1=kcenters.partial_predict(numpy.array(temp1))
    kcenters_TS2_sequences.append(temp1)

print "begin to plot the microstate implied timescale into the objective dir"
            #plot implied timescale

lag_times=range(10,100,10);  #adjust variables
n_timescales=5  #adjust variables

msm_timescales = implied_timescales(kcenters_TS2_sequences, lag_times, n_timescales=n_timescales, msm=MarkovStateModel(verbose=True, reversible_type='transpose'))

outfile_name="%s/TS2_ITS_tic%d_lagtime%d_clustersize%d.dat"%(outputdir, num_tics_for_clustering, tic_lag_time, nMicro)
print msm_timescales
print msm_timescales.shape

for k in range(n_timescales):
    plt.plot(lag_times, msm_timescales[:, k], 'o-')
f2=open(outfile_name, 'w')
for i in range(len(lag_times)):
    f2.write("%d    "%(lag_times[i]))
    for j in range(n_timescales):
        f2.write("%f    "%(msm_timescales[i, j]))
    f2.write('\n')
f2.close()
plt.title('Discrete-time MSM Relaxation Timescales')
plt.semilogy()
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,y1,100000))  #need to change this number
outfile_name="%s/TS2_ITS_tic%d_lagtime%d_clustersize%d.png"%(outputdir, num_tics_for_clustering, tic_lag_time, nMicro)
plt.savefig(outfile_name)
plt.close()

tmp_counter=0
TS2_assignment=outputdir+"/TS2_microAssignment/"
os.system("mkdir %s"%(TS2_assignment))
for ifile in range(len(TS_pair_list_array)):
    outfile_name="%s/%s_%d_%d_%d.txt"%(TS2_assignment, TS_pair_list_array[ifile], num_tics_for_clustering, tic_lag_time, nMicro)
    numpy.savetxt(outfile_name,kcenters_TS2_sequences[tmp_counter],fmt='%s')
    tmp_counter+=1


sys.exit()
