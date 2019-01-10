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

def readtrajs_from_folder(traj_folder):
    os.system("for j in %s/*.npy;do echo `basename $j .npy`;done>temp_trajlist"%(traj_folder))
    traj_list_array = []
    for line in open('temp_trajlist'):
        traj_list_array.append(line.strip())
    pairwise_distance = []
    for line in traj_list_array:
        pairwise_distance.append(numpy.load(traj_folder+line+'.npy'))  #npy files
    return traj_list_array, pairwise_distance

def fit_predict_tica_embeddings(traj_folder, reference_tica_information, outputdir, names, n_components, tica_lagtime):
    traj_list_array, pairwise_distance = readtrajs_from_folder(traj_folder)
    #calculate the tica means for the current system
    test_tica=tICA(lag_time=tica_lagtime, n_components=n_components)
    test_tica.fit(pairwise_distance)
    numpy.savetxt('%s/%s_pairwise_means'%(outputdir, names), test_tica.means_)

    for line in range(len(traj_list_array)):
        temp = numpy.load("%s/%s.npy"%(traj_folder, traj_list_array[line]))
        #we begin to project
        results_to_store = numpy.dot((temp-test_tica.means_.T),reference_tica_information.eigenvectors_[:, :])
        numpy.savetxt("%s/%s_ticproj.txt"%(outputdir, traj_list_array[line]), results_to_store[:, 0:n_components])

def read_sequences_txt(filename):
    microstate = []
    for line in open(filename):
        microstate.append(numpy.loadtxt(line.strip()))
    return microstate

def calculate_its(kcenters_sequences, lag_times, n_timescales, outfile_name, ergodic_cutoff_option):
	msm_timescales = implied_timescales(kcenters_sequences, lag_times, n_timescales=n_timescales, msm=MarkovStateModel(verbose=True, reversible_type='transpose', ergodic_cutoff=ergodic_cutoff_option))
	for k in range(n_timescales):
	    plt.plot(lag_times, msm_timescales[:, k], 'o-')
	f2=open(outfile_name+'.dat', 'w')
	for i in range(len(lag_times)):
	    f2.write("%d    "%(lag_times[i]))
	    for j in range(n_timescales):
	        f2.write("%f    "%(msm_timescales[i, j]))
	    f2.write('\n')
	f2.close()
	plt.title('Discrete-time MSM Relaxation Timescales')
	plt.semilogy()
	x1,x2,y1,y2 = plt.axis()
	plt.savefig(outfile_name+'.png')
	plt.close()
