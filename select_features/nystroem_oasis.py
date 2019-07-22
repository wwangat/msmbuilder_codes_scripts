import numpy as np
from pyemma.coordinates.transform.tica import *
from pyemma.coordinates.transform._tica_base import *
from pyemma.coordinates.transform.nystroem_tica import *

input_feature_data=[] #input features
for line in open('feature_trajlist.list'):
  print(line.strip())
  temp=np.loadtxt(line.strip())
  input_feature_data.append(temp)

lag=80 #tica correlation lag time
num_features = 18 #revise
for i in range(2,19,1): #modify
  max_columns=i #maximum # of features to use
  t = NystroemTICA(lag, max_columns, initial_columns=np.random.choice(num_features,1,replace=False), nsel=1) #add nsel features each time
  print("descriptions for this run:")
  print(t.describe())
  t.estimate(input_feature_data)   #######running oasis_tica
  print("initial columns are")
  print(t.initial_columns)
  print("the final selected feature indexes are:")
  print(t.column_indices)
  print("the estimated generalized eigenvalues using the selected columns")
  print(t.eigenvalues)
#  t = NystroemTICA(txx, lag, max_columns, initial_columns=np.random.choice(5000,1,replace=False), nsel=20)
  np.savetxt("t_timescales" + str(max_columns) +".txt", t.timescales) 
  ######saving the useful informations:eigenvalues and column indexes#######
  np.savetxt("selected_feature_indexes_maxcolumn"+str(max_columns)+".txt", t.column_indices, fmt='%d')
  np.savetxt("eigenvalues_based_on_selected_feature_indexes_maxcolumn"+str(max_columns)+".txt", t.eigenvalues)
