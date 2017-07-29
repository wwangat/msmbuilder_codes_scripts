#! /usr/bin/env python
###############################################################################
#      Filename:  test_pyemma_tpt.py
#       Created:  2016-08-19 01:09
#        Author:  Wei WANG        (wwangat@gmail.com)
#        Input : a row normalized microstate transition proabability matrix and the micro to macro mapping
###############################################################################

###############################################################################

"""
test the pyemma code on tpt calculation
"""

###############################################################################

# -*- coding: utf-8 -*-
# vim:fenc=utf-8

###############################################################################
import pyemma.plots as mplt
import pyemma.msm as msm
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-x',"--TPM",default = "TPM.txt", help="the microstate TPM, mle is preferred")
parser.add_option('-m',"--mapping_file",default = "mapping.txt", help="the microstate to macrostate mapping relationship, start from 0")
#parser.add_option('-o',"--resultdir", default='.', help="where to put the results")
parser.add_option('-i','--source',help='indicate the source, index start from 0',type='int')
parser.add_option('-e','--sink',help='indicate the sink, index start from 0',type='int')

(options, args) = parser.parse_args()
if options.source == "noInput" or options.sink == "noInput" or options.mapping_file=="noInput":
    print "Both source and sink state and micro2macro mapping should be indicated"

#reading the mapping relationship
mapping = np.loadtxt(options.mapping_file, dtype='int')
sink=options.sink
source=options.source
#P = np.array([[0.8,  0.15, 0.05,  0.0,  0.0],
#            [0.1,  0.75, 0.05, 0.05, 0.05],
#            [0.05,  0.1,  0.8,  0.0,  0.05],
#            [0.0,  0.2, 0.0,  0.8,  0.0],
#            [0.0,  0.02, 0.02, 0.0,  0.96]])

#ours
P=np.loadtxt(options.TPM)
print "row-normalized macro TPM: "
print P

print "###########################################################"
M = msm.markov_model(P)
#pos = np.array([[2.0,-1.5],[1,0],[2.0,1.5],[0.0,-1.5],[0.0,1.5]])
#mplt.plot_markov_model(M, pos=pos);

#state 1: fold state, stat 2: misfold state, state 3: unfold state, state 4: hairpin 1 state
A = [source]
B = [sink] 
tpt = msm.tpt(M, A, B)

# get tpt gross flux
F = tpt.gross_flux
print '**Flux matrix**: '
print F
print '**forward committor**: '
print tpt.committor
print '**backward committor**: '
print tpt.backward_committor
# we position states along the y-axis according to the commitor
#tptpos = np.array([tpt.committor, [0,0,0.5,-0.5,0]]).transpose()
#print '\n**Gross flux illustration**: '
#mplt.plot_flux(tpt, pos=tptpos, arrow_label_format="%10.4f", attribute_to_plot='gross_flux')

# get tpt net flux
Fp = tpt.net_flux
# or: tpt.flux (it's the same!)
print '**Net-Flux matrix**: '
print Fp
# visualize
#mplt.plot_flux(tpt, pos=tptpos, arrow_label_format="%10.4f", attribute_to_plot='net_flux')

print "###########################################################"

print 'Total TPT flux = ', tpt.total_flux
#print 'Rate from TPT flux = ', tpt.rate
#print 'A->B transition time = %f ps' % (1.0*0.1*200*1e-3/tpt.rate)

#print 'mfpt(unfold,fold) = %f '% (M.mfpt(2, 0)*150*0.1*1e-3) # same as above about the unit

#mplt.plot_flux(tpt, pos=tptpos, flux_scale=100.0/tpt.total_flux, arrow_label_format="%3.1f")
#ylabel("committor")

print "#############################################################"
#print "all the pathways from unfold state to fold state are:"
tpt.pathways()
(paths,pathfluxes) = tpt.pathways(fraction=0.95)
cumflux = 0
print "Path flux\t\t%path\t%of total\tpath"
for i in range(len(paths)):
    cumflux += pathfluxes[i]
    print pathfluxes[i],'\t','%3.1f'%(100.0*pathfluxes[i]/tpt.total_flux),'%\t','%3.1f'%(100.0*cumflux/tpt.total_flux),'%\t\t',paths[i] 

print "##############################################################"
#print "description: state 0 -> fold state, state 1 -> hairpin 1 state, state 2 -> unfold state, state 3 -> misfold state"
#mplt.plot_network(Fsubpercent, pos=tptpos, state_sizes=tpt.stationary_distribution, arrow_label_format="%3.1f")


