#! /usr/bin/env python
###############################################################################
#      Filename:  get_atompair_info.py
#       Created:  2019-01-23 18:54
#        Author:  Wei WANG        (wwangat@gmail.com)
###############################################################################

###############################################################################

"""

"""

###############################################################################

# -*- coding: utf-8 -*-
# vim:fenc=utf-8

###############################################################################


#input: the atoms we used to generate the pairwise distance pairs, atom number starts from 1
import numpy
import optparse

p = optparse.OptionParser()
p.add_option('--input', '-f', help="please input the file containing the atom information, starts from 1")
p.add_option('--output', '-o', help="the filename for pairwise distance index")
options, arguments = p.parse_args()

atoms=[]
for line in open(options.input):
    line=line.strip()
    #line[0]: atom number from 1, line[1]:atom name, line[2]:residue number
    atoms.append(line)

fp = open(options.output, 'w')
for j in range(len(atoms)):
    for k in range(j+1, len(atoms)):
        fp.write("%s %s"%(atoms[j], atoms[k]))
        fp.write('\n')
fp.close()
