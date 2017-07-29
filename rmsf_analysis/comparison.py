import numpy as np
import matplotlib.pyplot as plt 
import os

def hist_plot(id,num_bin):
    os.system("for j in GS-rmsf/rmsf-*.xvg;do sed -n '13,$p' $j | grep -w %s;done | awk '{print $2}'>temp1;for j in INT2-rmsf/rmsf-*.xvg;do sed -n '13,$p' $j | grep -w %s;done | awk '{print $2}'>temp2;for j in TS22-rmsf/rmsf-*.xvg;do sed -n '13,$p' $j | grep -w %s;done | awk '{print $2}'>temp3;"%(id,id,id))
    temp1=np.loadtxt('temp1');
    temp2=np.loadtxt('temp2');
    temp3=np.loadtxt('temp3');
    maxvalue=np.max([np.max(temp1), np.max(temp2), np.max(temp3)])
    minvalue=np.min([np.min(temp1), np.min(temp2), np.min(temp3)])
    interval=(maxvalue-minvalue)/num_bin
    nbins=np.arange(minvalue, maxvalue, interval)
    a=np.histogram(temp1, nbins)
    b=np.histogram(temp2, nbins)
    c=np.histogram(temp3, nbins)

    plt.plot(a[1][:-1], a[0]/sum(a[0]*1.0), b[1][:-1], b[0]/sum(b[0]*1.0), c[1][:-1],c[0]/sum(c[0]*1.0))
    plt.scatter(a[1][:-1], a[0]/sum(a[0]*1.0))
    plt.scatter(b[1][:-1], b[0]/sum(b[0]*1.0))
    plt.scatter(c[1][:-1],c[0]/sum(c[0]*1.0))
    plt.legend(['GS','INT', 'TS2'])
    plt.xlabel('RMSF(nm)')
    plt.ylabel('percentage of trajectories')
    plt.savefig('Res%s.png'%(id))
    plt.close()
num_bin=20
for j in range(375,428,1):
    hist_plot(j, num_bin)
