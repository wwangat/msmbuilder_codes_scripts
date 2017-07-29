"""
Author: Wei WANG (wwangat@gmail.com)
Function: Plot contour map
Input: A file that has the following format
    Microstate_trajectory_name    frameid(start fro 0)  microstate_id(start from 0)  macrostate_id(start from 0)  tIC1_coords    tIC2_coords    tIC3_coords

"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from optparse import OptionParser
def plot_sub(xcoor,ycoor,addition):
	scale_x = 50
	scale_y = 50

	upper_x = max(xcoor)
	lower_x = min(xcoor)
	upper_y = max(ycoor)
	lower_y = min(ycoor)

	print('upper_x: ', upper_x)
	print('lower_x: ', lower_x)
	print('upper_y: ', upper_y)
	print('lower_y: ', lower_y)

	step_x = (upper_x - lower_x) / (scale_x-1)
	step_y = (upper_y - lower_y) / (scale_y-1)
	X = []
	Y = []
	for i in range(scale_x):
		X.append(step_x*i+lower_x)
	for i in range(scale_y):
		Y.append(step_y*i+lower_y)

	density_matrix = []
	for i in range(scale_x):
		density_matrix.append([])
		for j in range(scale_y):
			density_matrix[i].append(0)

	point_list = []
	for i in range(len(xcoor)):
		point_list.append( (xcoor[i], ycoor[i]) )
	for p in point_list:
		i = int((p[0]-lower_x) / step_x)
		j = int((p[1]-lower_y) / step_y)
		density_matrix[i][j] += 1

	density_matrix = np.asarray(density_matrix)
	max_den = np.max(np.max(density_matrix))*1.0
	for i in range(len(density_matrix)):
		for j in range(len(density_matrix[i])):
			if (density_matrix[i][j] != 0):
				density_matrix[i][j] = np.log(density_matrix[i][j]/max_den)
			else:
				density_matrix[i][j]=-9
	print(density_matrix)
	print max_den
	plt.title('2D density picture')
	Z = np.asarray(density_matrix)
	plt.contour(X,Y,Z.T, 20, cmap=addition)

parser = OptionParser()
parser.add_option('-f',"--filename",default = "filename", help="the filename")
parser.add_option('-s','--tIC1',help='tIC1',type='int')
parser.add_option('-t','--tIC2',help='tIC2',type='int')

(options, args) = parser.parse_args()
if options.tIC1 == "noInput" or options.tIC2 == "noInput":
    print "Both tIC1 and tIC2 should be indicated"

#reading the mapping relationship
tIC1=options.tIC1
tIC2=options.tIC2
filename = options.filename

assign=[]
xcoor=[]
ycoor=[]
for line in file(filename):
	line=line.split()
	assign.append(int(line[3]))
	xcoor.append(float(line[4+tIC1]))
	ycoor.append(float(line[4+tIC2]))
assign = np.asarray(assign)
xcoor = np.asarray(xcoor)
ycoor = np.asarray(ycoor)
#now starting to make statistics of the frequency
nState=np.max(assign)

nbin=100

colorid=['Reds','Greens','Oranges','Blues','Purples','Greys']   #just add more colors here
plot_sub(xcoor,ycoor, colorid[5])
plt.xlabel('tIC%d'%(tIC1))
plt.ylabel('tIC%d'%(tIC2))

plt.title('Free energy landscape')
plt.savefig("%s_overall_tIC%d%d.png"%(filename,tIC1,tIC2))
plt.colorbar()
plt.close()



plt.figure()
for stateid in range(nState+1):
	print stateid
	index=np.where(assign==stateid)
	plot_sub(xcoor[index],ycoor[index],colorid[stateid])

plt.xlabel('tIC%d'%(tIC1))
plt.ylabel('tIC%d'%(tIC2))
plt.title('Free energy landscape')
plt.xlim([min(xcoor),max(xcoor)])
plt.ylim([min(ycoor),max(ycoor)])
plt.colorbar()
plt.savefig("%s_overall_split_tIC%d%d.png"%(filename,tIC1,tIC2))
plt.close()
