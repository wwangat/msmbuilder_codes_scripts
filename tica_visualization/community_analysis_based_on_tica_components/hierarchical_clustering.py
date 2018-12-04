############clustering the atoms into communities, generate the good initial clustering for the optimization code
import numpy as np
import os
from lof import LocalOutlierFactor

X= np.loadtxt('tic1_importance_matrix.mtx')
nState=4
option='ward'

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, single, fcluster

pDist = []
Y=X
for i in range(len(Y)):
        for j in range(i+1, len(Y)):
                pDist.append(Y[i][j])


pDist = np.array(pDist)
print(pDist.shape)
Z = linkage(pDist, option)

plt.figure(figsize = (25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')

dendrogram(
	Z,
	leaf_rotation = 0. , #rotate the x axis labels
	leaf_font_size = 10. , #font size for the x axis labels
)

#os.mkdir(option)
outputdir=option
plt.savefig(outputdir+'/'+option+'.png')
#plt.show()

cluster_result=fcluster(Z, nState, criterion='maxclust')
print cluster_result
final_result=cluster_result-1
print final_result
np.savetxt(outputdir+"/lumping_into_"+str(nState)+"states.txt", final_result, fmt='%d')
#numpy.savetxt(outputdir+"/lumping_into_"+str(nState)+"states.txt", fcluster(Z, nState, criterion='maxclust')-1, fmt='%d')

np.savetxt(outputdir+"/linkage.txt", Z)
