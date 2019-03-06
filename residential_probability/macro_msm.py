#code contributed by qiao qin
from numpy import *
import optparse
import sys
import linecache

def main():

        p = optparse.OptionParser()
        p.add_option('--headDir', '-d', default="./")
        p.add_option('--trajlist', '-t', default="-1")
        p.add_option('--dt','-i',default="2")
	p.add_option('--lagstep', '-l', default="1")
	p.add_option('--propagatenum', '-g', default="1")
        p.add_option('--mapping','-p',help="mapping file from micro to macro")
	p.add_option('--micronum','-k',help="microstate num")
	p.add_option('--bigstates','-s', help="file of bigstates to get the residence time", default="")
	p.add_option('--stepsahead', '-a', help="stepsahead for correction of recrossing")
        options, arguments = p.parse_args()
	propnum=int(options.propagatenum)
        dt=int(options.dt)
        tLag=int(options.lagstep)
        headDir=(options.headDir)
        trajlistfiles=(options.trajlist)
        filename=(options.mapping)
	stepsahead=int(options.stepsahead)
	micronum=int(options.micronum)
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

	tCount= getMicroTransitionsFromAssignements(tLag, headDir, trajlistfiles, stepsahead, micronum)
	tCount=(tCount+tCount.transpose())/2.0
	tCountMac=zeros((macronum+1, macronum+1))	
	indices=array(tCount.nonzero()).transpose()
	for k in range(len(indices)):
		micro1 = indices[k,0]
		micro2 = indices[k,1]
                if (mapMicroToMacro[micro1]!= -1 and mapMicroToMacro[micro2] != -1) :
                        tCountMac[mapMicroToMacro[micro1],mapMicroToMacro[micro2]] += tCount[micro1, micro2]

	del tCount
	del micro1
	del micro2
# normalize to get Prob
	savetxt("POPMACRO.DAT", tCountMac.sum(axis=1))
        savetxt("tCountMac.DAT", tCountMac)

	weight=array(tCountMac.sum(axis=1).flatten())
	weight[weight==0] = 1
	tProbMac = tCountMac / weight
	del tCountMac
	tProbMacProp=tProbMac.copy()
        for i in range(propnum):
		if i==0:
			tProbMacProp=tProbMac.copy()
		else:
			tProbMacProp=dot(tProbMacProp, tProbMac)
#		savetxt("tPROB_"+str(tLag*i*dt)+".DAT", tProbMacProp)
		if (len(bigstates)==0):
			for k in range(macronum+1):
				fr=open("residence_time_msm"+str(k)+".dat", "a")		
				fr.write("%d %f\n" %((i+1)*dt*tLag, tProbMacProp[k][k]))
				fr.close()
		else :
			for k in range(len(states)):
				fr=open("residence_time_msm"+str(states[k])+".dat", "a")
				fr.write("%d %f\n" %((i+1)*dt*tLag, tProbMacProp[states[k]][states[k]]))
				fr.close()
	del tProbMac
	del tProbMacProp

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
        if origState<0 or newState<0:
          # move to next window.
          # if subsample is true then move forward lagTime steps.
          # otherwise move forward one step (sliding window).
          if subsample:
            i += tLag
          else:
            i += 1
          continue

        # add count
        tCount[origState,newState] += 1

        # move to next window.
        # sliding window
        #i += tLag
	#i += 1
        WINDOW
    return tCount


if __name__ == '__main__':
 main()

