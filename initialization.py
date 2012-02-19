#Perform some initial setup
import numpy
import scipy
import pylab
import time

import config as _conf
from config import *
from geneSimClasses import *

def checkCollision(newNode, existingNodes, nodeSize):
    sizeSquared = nodeSize**2 #avoid performing the same calculation repeatedly!
    for node in existingNodes:
        if (((newNode[0]-node[0])**2 + (newNode[1]-node[1])**2 + (newNode[2]-node[2])**2) < sizeSquared):
            print "collision averted"
            return True
    return False

def placeBindSites(bindSites, clusters=False, cSize=0, cProbability=0):
    i=0
    if clusters==False:
        while (i<len(bindSites)):
            bindSites[i][0] = numpy.random.uniform(-maxX + TFsize/2.0,maxX-TFsize/2.0)
            bindSites[i][1] = numpy.random.uniform(-maxY+TFsize/2.0,maxY-TFsize/2.0)
            bindSites[i][2] = numpy.random.uniform(-maxZ+TFsize/2.0,maxZ-TFsize/2.0)
            bindSites[i][3] = -1 #all initially unbound
            bindSites[i][4] = False #none are in clusters here
            if not checkCollision(bindSites[i],bindSites[0:i],bindSize):
                i+=1
    else:
        print "Running With with clustering and a clusterprobability of " + str(clusterProbability)
        clusterX = numpy.random.uniform(-maxX + clusterSize/2.0,maxX-clusterSize/2.0)
        clusterY = numpy.random.uniform(-maxY + clusterSize/2.0,maxY-clusterSize/2.0)
        clusterZ = numpy.random.uniform(-maxZ + clusterSize/2.0,maxZ-clusterSize/2.0)
        while (i<len(bindSites)):
            rand = numpy.random.uniform(0,1)
            if (rand <= cProbability): #we place this one in a cluster!
                placed = 0
                while (placed==0):
                    bindSites[i][0] = clusterX + numpy.random.uniform(-clusterSize/2.0,+clusterSize/2.0)
                    bindSites[i][1] = clusterY + numpy.random.uniform(-clusterSize/2.0,+clusterSize/2.0)
                    bindSites[i][2] = clusterZ + numpy.random.uniform(-clusterSize/2.0,+clusterSize/2.0)
                    bindSites[i][3] = -1 #all initially unbound
                    bindSites[i][4] = True #in a cluster
                    if not checkCollision(bindSites[i],bindSites[0:i],bindSize):
                        placed=1
            else:
                placed = 0
                while (placed==0):
                    bindSites[i][0] = numpy.random.uniform(-maxX + TFsize/2.0,maxX-TFsize/2.0)
                    bindSites[i][1] = numpy.random.uniform(-maxY+TFsize/2.0,maxY-TFsize/2.0)
                    bindSites[i][2] = numpy.random.uniform(-maxZ+TFsize/2.0,maxZ-TFsize/2.0)
                    bindSites[i][3] = -1 #all initially unbound
                    bindSites[i][4] = False #NOT in a cluster
                    if not checkCollision(bindSites[i],bindSites[0:i],bindSize):
                        placed=1
            i+=1
        return (clusterX,clusterY,clusterZ)
    return

def placeTFs(tfs, size):
    i=0
    while (i<len(tfs)):
        tfs[i][0] = numpy.random.uniform(-maxX + TFsize/2.0,maxX-TFsize/2.0)
        tfs[i][1] = numpy.random.uniform(-maxY+TFsize/2.0,maxY-TFsize/2.0)
        tfs[i][2] = numpy.random.uniform(-maxZ+TFsize/2.0,maxZ-TFsize/2.0)
        tfs[i][3] = -1 #set them all unbound initially
        i+=1
    return

def checkDiffusionCollision(node, allNodes, allNodesPrev):
    #do nothing here now, yet... maybe never at all, since collisions between TFs are kind of irrevelant
    return

def diffuseTFs(TFs,bindSites,dt,diffC):
    #treating each dimension as a gaussian probability of movement, with sigma = (2Dt)^.5
    sigma = (2*dt*(10**-6)*diffC)**.5 #The 10^-6 is due to dt being in _microseconds_

    #avoid a lot of menial, repetitive calculations just for comparison purposes
    maxXminusHalfSize = maxX - TFsize/2.0
    minXplusHalfSize = -maxX + TFsize/2.0
    maxYminusHalfSize = maxY - TFsize/2.0
    minYplusHalfSize = -maxY + TFsize/2.0
    maxZminusHalfSize = maxZ - TFsize/2.0
    minZplusHalfSize = -maxZ + TFsize/2.0

    degreesOfFreedom = 3*len(TFs)

    randArr = numpy.random.normal(0,sigma,3*len(TFs))

    i=0
    tf=0
    while (i<degreesOfFreedom):
        if (TFs[tf][3] == -1): #only currently unbound elements may move!
            TFs[tf][0] += randArr[i]
            if (TFs[tf][0] > maxXminusHalfSize): #prevent leaving the box!
                TFs[tf][0] = 2.0*maxXminusHalfSize - TFs[tf][0]
            if (TFs[tf][0] < minXplusHalfSize): #prevent leaving the box!
                TFs[tf][0] = 2.0*minXplusHalfSize - TFs[tf][0]
            TFs[tf][1] += randArr[i+1]
            if (TFs[tf][1] > maxYminusHalfSize): #prevent leaving the box!
                TFs[tf][1] = 2.0*maxYminusHalfSize - TFs[tf][1]
            if (TFs[tf][1] < minYplusHalfSize): #prevent leaving the box!
                TFs[tf][1] = 2.0*minYplusHalfSize - TFs[tf][1]
            TFs[tf][2] += randArr[i+2]
            if (TFs[tf][2] > maxZminusHalfSize): #prevent leaving the box!
                TFs[tf][2] = 2.0*maxZminusHalfSize - TFs[tf][2]
            if (TFs[tf][2] < minZplusHalfSize): #prevent leaving the box!
                TFs[tf][2] = 2.0*minZplusHalfSize - TFs[tf][2]
        i+=3
        tf+=1
    return

def iterate():
    global TFarray
    #TFarrayPrev = TFarray.copy() #dont think we need to know this at all, currently
    TFunbind()
    unBoundTFs = filter(lambda x: x[3] != -1, TFarray)
    diffuseTFs(unBoundTFs,bindSitesPosn,delta_t,TFdiffusionC)
    TFbind(unBoundTFs)
    return

def TFunbind(): #allow TFs to become unbound based on their probability of doing so
    i=0

    numTFs = len(TFarray)
    c = numpy.random.uniform(0.0,1.0, numTFs)

    while (i<numTFs):
        if (TFarray[i][3] !=-1 ):
            if (c[i] <= pUnbind):
                bindSitesPosn[TFarray[i][3]][3] = -1 #unbind site
                siteNum = TFarray[i][3]
                TFarray[i][3] = -1 #unbind TF
                print "UNBIND!"
                if (i == tfToTrack):
                    writeTFStatus(t + delta_t , i)
                if (clustering and clusterProbability<1):
                    if (siteNum == nonClusterSiteToTrack):
                        writeSiteStatus(t + delta_t , siteNum, False) # False for the non-cluster site
                    elif (siteNum == clusterSiteToTrack):
                        writeSiteStatus(t + delta_t , siteNum)
                else:
                    if (siteNum == siteToTrack):
                        writeSiteStatus(t + delta_t , siteNum)
        i+=1
    return

def TFbind(tfs):
    i=0
    numTFs = len(tfs)
    while (i<numTFs):
        siteNum = checkTFsiteCollision(tfs[i],bindSitesPosn)
        if (siteNum != -1):
            c = numpy.random.uniform(0.0,1.0)
            if (c <= pBind):
                bindSitesPosn[siteNum][3] = i #bind site
                tfs[i][3] = siteNum #bind TF
                print "BIND!"
                if (i == tfToTrack):
                    writeTFStatus(t + delta_t , i)
                if (clustering and uSitesNum):
                    if (siteNum == nonClusterSiteToTrack):
                        writeSiteStatus(t + delta_t , siteNum, False) # False for the non-cluster site
                    elif (siteNum == clusterSiteToTrack):
                        writeSiteStatus(t + delta_t , siteNum)
                else:
                    if (siteNum == siteToTrack):
                        writeSiteStatus(t + delta_t , siteNum)

        i+=1
    return

def checkTFsiteCollision(tf,bindSites):
    i=0
    while (i<len(bindSites)):
        if (bindSitesPosn[siteNum][3] != -1): # ensure only unbound sites are considered for binding
            distSquared = (bindSites[i][0] - tf[0])**2 + (bindSites[i][1] - tf[1])**2 + (bindSites[i][2] - tf[2])**2
            if (distSquared <= bindDistanceSquared):
                return i
        i+=1
    return -1 #no collision

def calcGlobalStats():
    stats = {} #an associative array
    #Percentage of TFs Bound
    tfNum = len(TFarray) #in case we allow them to ever be created or destroyed
    tfBound = 0
    tfx = 0
    tfy = 0
    tfz = 0
    for tf in TFarray:
        tfx += tf[0] #for center of mass of the TFs calculation
        tfy += tf[1]
        tfz += tf[2]
        if (tf[3] != -1):
            tfBound +=1
    stats['tfPercentBound'] = float(tfBound)/tfNum
    #Percentage of Binding Sites Occupied
    bNum = len(bindSitesPosn) #in case we allow them to ever be created or destroyed (less likely than with TFs)
    bBound = 0
    for b in bindSitesPosn:
        if (b[3] != -1):
            bBound +=1
    stats['bindPercentOccupied'] = float(bBound) / bNum
    #TF Center of Mass
    cm = numpy.zeros((3))
    cm[0] = tfx / tfNum
    cm[1] = tfy / tfNum
    cm[2] = tfz / tfNum
    stats['tfmu-x'] = cm[0]
    stats['tfmu-y'] = cm[1]
    stats['tfmu-z'] = cm[2]
    return stats
#to output results, x = calcGlobalStats(); for i in x {print str(i)+':' + str(x[i])}

def writeSysStatus(time,verbose=False): #output current systen status to a file -- need to think about output format and what data is useful
    #verbose is intended to mean write out the status of the entire system
    # (as in all positions of binding sites and tfs), but I'm not sure how to
    # make that gnuplot-friendly
    sysStats = calcGlobalStats()
    if (time==0):
        statsComment = '#timestamp(mu-s)'
    statsValues = str(time)
    statNames = sysStats.keys()
    statNames.sort()
    for statName in statNames:
        if (time==0):
            statsComment += ','+statName
        statsValues += '\t'+str(sysStats[statName])
    if (time==0):
        systemFileOut.write(statsComment+'\n')
    systemFileOut.write(statsValues+'\n')
    return

def writeTFStatus(time,tf):
    if (time==0):
        tfComment = '#TF being tracked: #'+str(tf)+' (of '+str(len(TFarray))+')\n#timestamp(mu-s),x,y,z,status,boundTo(-1 means unbound!)'
    tfVals = str(time)+'\t'+str(TFarray[tf][0])+'\t'+str(TFarray[tf][1])+'\t'+str(TFarray[tf][2])+'\t'+str(int(TFarray[tf][3]!=-1))+'\t'+str(TFarray[tf][3])
    if (time==0):
        tfFileOut.write(tfComment+'\n')
    tfFileOut.write(tfVals+'\n')
    return

def writeSiteStatus(time,site, c = True): # c is False for the non-cluster site, when applicable
    if (c):
        fOut = bindSiteFileOut
    else:
        fOut = bindSiteFileOut2
    if (time==0):
        siteComment = '#Site being tracked: #'+str(site)+' (of '+str(len(bindSitesPosn))+')\n'
        siteComment += '#located at position: '+str(bindSitesPosn[site][0])+'\t'+str(bindSitesPosn[site][1])+'\t'+str(bindSitesPosn[site][2])+'\n'
        if bindSitesPosn[site][4]:
            siteComment += 'site is within a cluster'
        siteComment += '#timestamp(mu-s),status,boundTo(-1 means unbound!)'
        fOut.write(siteComment+'\n')

    siteVals = str(time)+'\t'+str(int(bindSitesPosn[site][3]!=-1))+'\t'+str(bindSitesPosn[site][3])
    fOut.write(siteVals+'\n')
    return

def printStatus(time, verbose=False):
    print "System simulation %.2f%% complete (t=%d) " % (((100.*time)/duration), time)
    if (verbose):
        print "\n TF Info:"
        for tf in TFarray:
            print tf[0:4]

def writeConfig(f):
    f.write('### Configuration used to generate this run ###\n')
    var = dir(_conf)
    var.sort(key=str.upper)
    for v in var:
        if (v[0] != '_'):
            f.write('# ' + v + ': ' + str(eval(v)) + '\n') #output the name of v, then the value of v
    f.write('\n\n')
    return

def runSimulation():
    global t
    iteration = 0

    while (t < duration):
        if (iteration % printOutAt == 0):
            printStatus(t, False)
        if (iteration % writeSysOutAt == 0):
            writeSysStatus(t)
        iterate();
        t += delta_t
        iteration += 1
    print "Simulation Complete"
    return

#initial variable settings
numpy.random.seed(seed)

maxX = systemWidth/2
maxY = systemLength/2
maxZ = systemHeight/2

#the array containing all the positions of the binding sites
bindSitesPosn = numpy.zeros((bindSitesNum,5)) # values 0, 1, and 2 are posn values. value 3 is the boundStatus (-1 is not bound!)
TFarray = numpy.zeros((TFnumberInit,4)) # values 0, 1, and 2 are posn values. value 3 is the boundStatus (-1 is not bound!)

#create initial files and write the config to the top of each, to ensure the config and data are never separated
start_t = time.localtime() #all files will share the same timestamp as part of their name
timestamp = str(start_t.tm_mon)+'-'+str(start_t.tm_mday)+'-'+str(start_t.tm_year)+'-h'+str(start_t.tm_hour)+'m'+str(start_t.tm_min)+'s'+str(start_t.tm_sec)
systemFileOut = open(fileOut+'-system-'+timestamp,'w',0)
tfFileOut = open(fileOut+'-TF-'+timestamp,'w',0)
bindSiteFileOut = open(fileOut+'-bindingSite-'+timestamp,'w',0)
if (clustering and clusterProbability<1):
    bindSiteFileOut2 = open(fileOut+'-bindingSite(non-cluster)-'+timestamp,'w',0)
    writeConfig(bindSiteFileOut2)
writeConfig(systemFileOut)
writeConfig(tfFileOut)
writeConfig(bindSiteFileOut)

t=0
try:
    #The actual flow of the program occurs below.
    clusterPosn = placeBindSites(bindSitesPosn,clustering,clusterSize,clusterProbability)
    if (clusterPosn):
        systemFileOut.write("# Cluster centered at: " + str(clusterPosn) + "\n\n")
    placeTFs(TFarray,TFsize)

    unclusteredSites = filter(lambda x: x[4] == False, bindSitesPosn)
    uSitesNum = len(unclusteredSites)

    print "TF's and Binding Sites now in place."

    if (clustering):
        print "Cluster centered at: " + str(clusterPosn) + "\n"

    print "Binding Sites are fixed at the following locations:"
    systemFileOut.write("# Binding Sites are fixed at the following locations:\n")
    for site in bindSitesPosn:
        print site[0:3]
        systemFileOut.write("# " + str(site) + "\n")
    systemFileOut.write("\n\n")
    print "\n TFs are initially located here:"
    for tf in TFarray:
        print tf[0:3]

    print "Beginning time evolution..."

    tfToTrack = numpy.random.randint(len(TFarray))
    writeTFStatus(0,tfToTrack)
    if (clustering and clusterProbability<1):
        clusterSiteToTrack = numpy.random.randint(len(bindSitesPosn))
        nonClusterSiteToTrack = numpy.random.randint(len(bindSitesPosn))
        while (not bindSitesPosn[clusterSiteToTrack][4]):
            clusterSiteToTrack = numpy.random.randint(len(bindSitesPosn))
        writeSiteStatus(0, clusterSiteToTrack)
        while (bindSitesPosn[nonClusterSiteToTrack][4]):
            nonClusterSiteToTrack = numpy.random.randint(len(bindSitesPosn))
        writeSiteStatus(0, nonClusterSiteToTrack, False)
    else:
        siteToTrack = numpy.random.randint(len(bindSitesPosn))
        writeSiteStatus(0, siteToTrack)


    runSimulation()
finally: #this way the files are closed even if an error is encountered
    systemFileOut.close()
    tfFileOut.close()
    bindSiteFileOut.close()