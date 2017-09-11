import os
import os.path
from numpy import *
import optparse
import sys
import linecache

def main():

        p = optparse.OptionParser()
        p.add_option('--lagstep', '-l', default="1")
        p.add_option('--headDir', '-d', default="./")
        p.add_option('--trajlist', '-t', default="-1")
        p.add_option('--maxLagstep', '-m', action="store", type="int", help="(required) Maximum lag time to go up to.")
        p.add_option('--dt','-i',default="2")
        p.add_option('--mapping','-p',help="mapping file from micro to macro")   
	p.add_option('--micronum','-k',help="microstate num")
	p.add_option('--bigstates','-s', help="file of bigstates to get the residence time", default="")
	p.add_option('--stepsahead', '-a', help = "the steps to correct recrossing", default="0")
        options, arguments = p.parse_args()
        dt=int(options.dt)
        tLag=int(options.lagstep)
        headDir=(options.headDir)
        trajlistfiles=(options.trajlist)
        filename=(options.mapping)
        tMaxLag=int(options.maxLagstep)
	micronum=int(options.micronum)
	stepsahead=int(options.stepsahead)
	bigstates=options.bigstates
	states=[]
	if (len(bigstates)!=0):
		fp=open(bigstates,"r")
		for line in fp:
			states.append(int(float(line)))
		fp.close
		print states

	mapMicroToMacro = zeros(micronum)
        fq=open(str(filename),"r")
        j=0
        macronum=0
        for line in fq:
                line = line.strip()
                mapMicroToMacro[j]=int(line)
                if int(line) > macronum:
                        macronum=int(line)
                j+=1
        fq.close
        #print "macronum", macronum

        for i in range(1,tMaxLag/tLag +1):
                tCount = getMicroTransitionsFromAssignements(tLag*i, headDir, trajlistfiles, stepsahead, micronum)
                tCount=(tCount+tCount.transpose())/2.0
		tCountMac=zeros((macronum+1, macronum+1))
		indices=array(tCount.nonzero()).transpose()
		for k in range(len(indices)):
			micro1 = indices[k,0]
			micro2 = indices[k,1]
			if (mapMicroToMacro[micro1]!= -1 and mapMicroToMacro[micro2] != -1) :
				tCountMac[mapMicroToMacro[micro1], mapMicroToMacro[micro2]] += tCount[micro1, micro2]
		del tCount
#		savetxt("tCountMac_"+str(tLag*i*dt)+".dat", tCountMac)
		savetxt("POP.DAT", tCountMac.sum(axis=1))
		if (len(bigstates)==0):
			for k in range(macronum+1):
				p=tCountMac[k][k]/tCountMac.sum(axis=1)[k]
				error= sqrt(i*(p-p*p)/tCountMac.sum(axis=1)[k])
				fr=open("residence_time_md_"+str(k)+".dat", "a")		
				fr.write("%d %f %f\n" %(i*dt*tLag, p, error))
				fr.close()
		else :
			for k in range(len(states)):
				p=tCountMac[states[k]][states[k]]/tCountMac.sum(axis=1)[states[k]]
				error= sqrt(i*(p-p*p)/tCountMac.sum(axis=1)[states[k]])
				fr=open("residence_time_md_"+str(states[k])+".dat", "a")
				fr.write("%d %f %f\n" %(i*dt*tLag, p, error))
				fr.close()
## add the error calculation part


		del tCountMac
		del micro1
		del micro2

def getMicroTransitionsFromAssignements(tLag, headDir, trajlistfiles, stepsahead, micronum):
    """Get the transition matrix.

    ARGUMENTS:
      subsample = If false use sliding window to count transitions, if true use every lagTime'th step (bool).
      useNoe = Determine whether or not to use Noe method to get transition probability matrix (bool).  If useNoe is true then subsample should also be true in order to get accurate statistics.
      nIter = Number iteration Noe sampling to do, only matters if useNoe is true (int).
      freqSample = Frequency to store samples from Noe Markov chain, only matters if useNoe is true (int).  Stores every freqSample'th step.
    """
    tCount=zeros((micronum, micronum))
    # read assignment files to determine transition count matrix
    for trajFileName in file(trajlistfiles):
      # skip if file doesn't exist
      # read file
      trajFileName = trajFileName.strip()
      trajFileName = "%s/assignments/%s" % (headDir,trajFileName)
      f = open(trajFileName, 'r')
      assignments = f.readlines()
      f.close()

      # get transition counts for given lag time
      trajLen = len(assignments)
      i = 0
#========================================================= original
#      while i < trajLen-tLag:
#        origState = self.getState(assignments[i].strip().split())
#        newState = self.getState(assignments[i+tLag].strip().split())
#========================================================= end of original


#========================================================= add correction of recrossing 
      trajLen = len(assignments)
      i = stepsahead  #Start ahead
      while i < (trajLen - (tLag + stepsahead)):
        origState = int((assignments[i].strip().split())[0])
        newState = int((assignments[i+tLag].strip().split())[0])
        avoidCountIt=0
        for jj in range((i+tLag), (i+tLag+stepsahead)):
           tmpstate= int((assignments[jj].strip().split())[0])
           if (tmpstate != newState):
             avoidCountIt=1
           if (avoidCountIt == 1):
             break
#========================================================= end of recrossing

        # state will be -1 if using ST flags in hierarchical clustering and temperature not in desired range.
        # this has no effect if not using ST

        # add count
        tCount[origState,newState] += 1 - avoidCountIt

        # move to next window.
        # sliding window
        #i += tLag
        #i += 100
        #i += 400
        WINDOW
    return tCount


if __name__ == '__main__':
 main()

