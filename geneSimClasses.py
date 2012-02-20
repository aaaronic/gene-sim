# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 23:26:17 2012

@author: Ahaarnos
"""
import numpy as np
import time



class BindingSite:
    def __init__(self, siteNum, size, posn=(0,0,0), boundTo=None, inCluster=False):
        self.siteNum=siteNum
        self.size=size
        self.x=posn[0]
        self.y=posn[1]
        self.z=posn[2]
        self.boundTo=boundTo
        self.inCluster=inCluster
        self.file=None
    
    def __str__(self):
        return str(self.x)+'\t'+str(self.y)+'\t'+str(self.z)+'\t'+'\t'+str(self.boundTo)

    def siteOverlap(self, otherSites):
        for site in otherSites:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < (self.size+site.size)**2):
                print "Avoiding site collision in placement"
                return True
        return False
    
    def writeStatus(self, time, numSites):
        if (self.file is not None): # make sure there's something to write to before continuing
            if (time==0):
                comment = '#Site being tracked: #'+str(self.siteNum)+' (of '+str(numSites)+')\n'
                comment += '#located at position: '+str(self.x)+'\t'+str(self.y)+'\t'+str(self.z)+'\n'
                if (self.inCluster):
                    comment += 'site is within a cluster'
                comment += '#timestamp(mu-s),status(0/1 for unbound/bound),boundTo(None means unbound!)'
                self.file.write(comment+'\n')
            tfNum = -1
            if (self.boundTo is not None):
                tfNum = self.boundTo.tfNum
            siteVals = str(time)+'\t'+str(self.boundTo is not None)+'\t'+str(tfNum)
            self.file.write(siteVals+'\n')



class TF:
    def __init__(self, tfNum, size, posn=(0,0,0), boundTo=None):
        self.tfNum=tfNum
        self.size=size
        self.x=posn[0]
        self.y=posn[1]
        self.z=posn[2]
        self.prevX=posn[0]
        self.prevY=posn[1]
        self.prevZ=posn[2]
        self.boundTo=boundTo
        self.file=None

    def __str__(self):
        return str(self.x)+'\t'+str(self.y)+'\t'+str(self.z)+'\t'+'\t'+str(self.boundTo)

    def tfCollision(self, otherTFs, tfSize): #not used for now... maybe never
        sizeSquared = (2*tfSize)**2
        for site in otherTFs:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < sizeSquared):
                print "Avoiding site collision in placement"
                return True
        return False
    
    def checkTFsiteCollision(self, bindSites, bindDistanceSquared):
        for s in bindSites:
            if (s.boundTo is None): # ensure only unbound sites are considered for binding
                distSquared = (s.x - self.x)**2
                if (distSquared < bindDistanceSquared): #check each coordinate individually before bothering to continue to the full comparison
                    distSquared += (s.y - self.y)**2 #dy**2
                    if (distSquared < bindDistanceSquared):
                        distSquared += (s.z - self.z)**2
                        if (distSquared <= bindDistanceSquared):
                            return s
        return None #no collision
    
    def writeStatus(self, time, numTFs):
        if (self.file is not None): # make sure there's something to write to before continuing
            if (time==0):
                comment = '#TF being tracked: #'+str(self.tfNum)+' (of '+str(numTFs)+')\n#timestamp(mu-s),x,y,z,status(1 is bound),boundTo(-1 is unbound)'
            siteNum = -1
            if (self.boundTo is not None):
                siteNum = self.boundTo.siteNum
            tfVals = str(time)+'\t'+str(self.x)+'\t'+str(self.y)+'\t'+str(self.z)+'\t'+str(int(self.boundTo is not None))+'\t'+str(siteNum)
            if (time==0):
                self.file.write(comment+'\n')
            self.file.write(tfVals+'\n')



class Simulation:
    def __init__(self, config):
        np.random.seed(config.seed)

        self.time=0
        self.iteration=0
        self.config=config
        self.tfs=[]
        self.sites=[]
        self.cluster=None

        self.start_t = time.localtime() #all files will share the same timestamp as part of their name
        self.timestamp = str(self.start_t.tm_mon)+'-'+str(self.start_t.tm_mday)+'-'+str(self.start_t.tm_year)+'-h'+str(self.start_t.tm_hour)+'m'+str(self.start_t.tm_min)+'s'+str(self.start_t.tm_sec)
        self.fileOut = open(config.fileOut+'-system-'+self.timestamp,'w',0)
        self.writeConfig(self.fileOut)

        maxX = config.systemWidth/2
        self.maxX = maxX
        maxY = config.systemLength/2
        self.maxY = maxY
        maxZ = config.systemHeight/2
        self.maxZ = maxZ

        for i in range(config.TFnumberInit):
            xPosn = np.random.uniform(-maxX + config.TFsize/2.0,maxX-config.TFsize/2.0)
            yPosn = np.random.uniform(-maxY + config.TFsize/2.0,maxY-config.TFsize/2.0)
            zPosn = np.random.uniform(-maxZ + config.TFsize/2.0,maxZ-config.TFsize/2.0)
            self.tfs.append(TF(i, config.TFsize, (xPosn,yPosn,zPosn)))

        tfToTrack = np.random.randint(len(self.tfs))
        self.tfs[tfToTrack].file = open(self.config.fileOut+'-TF-'+self.timestamp,'w',0)
        self.writeConfig(self.tfs[tfToTrack].file)
        self.tfs[tfToTrack].writeStatus(0,len(self.tfs))

        ## Place the binding sites either in a cluster or not, depending on the config
        i=0
        if (not config.clusterProbability):
            while(i<config.bindSitesNum):
                xPosn = np.random.uniform(-maxX + config.bindSize/2.0,maxX-config.bindSize/2.0)
                yPosn = np.random.uniform(-maxY + config.bindSize/2.0,maxY-config.bindSize/2.0)
                zPosn = np.random.uniform(-maxZ + config.bindSize/2.0,maxZ-config.bindSize/2.0)
                site = BindingSite(i, config.bindSize, (xPosn,yPosn,zPosn))
                if (not site.siteOverlap(self.sites)):
                    self.sites.append(site)
                    i+=1 #now we've placed it and can move on to the next one
            siteToTrack = np.random.randint(len(self.sites))
            self.sites[siteToTrack].file = open(self.config.fileOut+'-bindingSite-'+str(siteToTrack)+'-'+self.timestamp,'w',0)
            self.writeConfig(self.sites[siteToTrack].file)
            self.sites[siteToTrack].writeStatus(0, len(self.sites))
        else:
            clusterX = np.random.uniform(-maxX + config.clusterSize/2.0,maxX-config.clusterSize/2.0)
            clusterY = np.random.uniform(-maxY + config.clusterSize/2.0,maxY-config.clusterSize/2.0)
            clusterZ = np.random.uniform(-maxZ + config.clusterSize/2.0,maxZ-config.clusterSize/2.0)
            self.cluster = SiteCluster(config.clusterSize, (clusterX, clusterY, clusterZ))
            print "Running With with a clusterprobability of " + str(config.clusterProbability)
            print "Cluster centered at:\n" + str(self.cluster.x)+','+str(self.cluster.y)+','+str(self.cluster.z)
            while (i<config.bindSitesNum):
                rand = np.random.uniform(0,1)
                if (rand <= config.clusterProbability): #we place this one in a cluster!
                    placed = False
                    while (not placed): #we keep the clustering decision while looking for a place to put it
                        xPosn = clusterX + np.random.uniform(-config.clusterSize/2.0,+config.clusterSize/2.0)
                        yPosn = clusterY + np.random.uniform(-config.clusterSize/2.0,+config.clusterSize/2.0)
                        zPosn = clusterZ + np.random.uniform(-config.clusterSize/2.0,+config.clusterSize/2.0)
                        site = BindingSite(i, config.bindSize, (xPosn,yPosn,zPosn), None, True)
                        if (not site.siteOverlap(self.sites)):
                            self.sites.append(site)
                            self.cluster.addSite(site)
                            placed=True
                else: #not in a cluster
                    placed = False
                    while (not placed): #we keep the clustering decision while looking for a place to put it
                        xPosn = np.random.uniform(-maxX + config.bindSize/2.0,maxX-config.bindSize/2.0)
                        yPosn = np.random.uniform(-maxY + config.bindSize/2.0,maxY-config.bindSize/2.0)
                        zPosn = np.random.uniform(-maxZ + config.bindSize/2.0,maxZ-config.bindSize/2.0)
                        site = BindingSite(i, config.bindSize, (xPosn,yPosn,zPosn))
                        if (not site.siteOverlap(self.sites)):
                            self.sites.append(site)
                            placed=True
                i+=1
            c = True
            while (c):
                siteToTrack = np.random.randint(len(self.sites))
                c = self.sites[siteToTrack].inCluster #we found one not it a cluster that we'll be tracking
            self.sites[siteToTrack].file = open(self.config.fileOut+'-bindingSite-'+str(siteToTrack)+'-nonCluster-'+self.timestamp,'w',0)
            self.writeConfig(self.sites[siteToTrack].file)
            self.sites[siteToTrack].writeStatus(0, len(self.sites))

            clusterSiteToTrack = self.cluster.sites[np.random.randint(len(self.cluster.sites))].siteNum
            self.sites[clusterSiteToTrack].file = open(self.config.fileOut+'-bindingSite-'+str(siteToTrack)+'-inCluster-'+self.timestamp,'w',0)
            self.writeConfig(self.sites[clusterSiteToTrack].file)
            self.sites[clusterSiteToTrack].writeStatus(0, len(self.sites))
        print "TF's and Binding Sites now in place."

        print "Binding Sites are fixed at the following locations:"
        self.fileOut.write("# Binding Sites are fixed at the following locations:\n")
        for site in self.sites:
            print str(site.x)+","+str(site.y)+","+str(site.z)
            self.fileOut.write("# " + str(site) + "\n")
        self.fileOut.write("\n\n")
        print "\n TFs are initially located here:"
        for tf in self.tfs:
            print str(tf.x)+","+str(tf.y)+","+str(tf.z)

    def runSimulation(self):
        while (self.time < self.config.duration):
            if (self.iteration % self.config.printOutAt == 0):
                self.printStatus(False)
            if (self.iteration % self.config.writeSysOutAt == 0):
                self.writeStatus(False)
            self.iterate();
            self.time += self.config.delta_t
            self.iteration += 1
        print "Simulation Complete"

    def writeConfig(self, f):
        f.write('### Configuration used to generate this run ###\n')
        var = dir(self.config)
        var.sort(key=str.upper)
        for v in var:
            if (v[0] != '_'):
                f.write('# ' + v + ': ' + str(eval('self.config.'+str(v))) + '\n') #output the name of v, then the value of v
        f.write('\n\n')

    def iterate(self):
        self.unbindTFs()
        unBoundTFs = filter(lambda x: x.boundTo is None, self.tfs)
        self.diffuseTFs(unBoundTFs)
        self.bindTFs(unBoundTFs)

    def bindTFs(self, tfs):
        i=0
        c = np.random.uniform(0.0,1.0, len(tfs))
        
        for tf in tfs:
            site = tf.checkTFsiteCollision(self.sites,self.config.bindDistanceSquared)
            if (site is not None):
                if (c[i] <= self.config.pBind):
                    site.boundTo = tf #bind site
                    tf = site #bind TF
                    print "BIND!"
                    if (tf.file is not None):
                        self.writeTFStatus(self.time + self.config.delta_t, tf)
                    if (site.file is not None):
                        self.writeSiteStatus(self.time + self.config.delta_t, site)
            i+=1

    def unbindTFs(self):
        i=0
        c = np.random.uniform(0.0,1.0, len(self.tfs))

        for tf in self.tfs:
            if (tf.boundTo is not None):
                if (c[i] <= self.config.pUnbind):
                    site = tf.boundTo
                    tf.boundTo.boundTo = None
                    tf.boundTo = None
                    print "UNBIND!"
                    if (tf.file is not None):
                        self.writeTFStatus(self.time + self.config.delta_t, tf)
                    if (site.file is not None):
                            self.writeSiteStatus(self.time + self.config.delta_t, site)
            i+=1
    
    def diffuseTFs(self, tfs):
        #treating each dimension as a gaussian probability of movement, with sigma = (2Dt)^.5
        sigma = (2*self.config.delta_t*(10**-6)*self.config.TFdiffusionC)**.5 #The 10^-6 is due to dt being in _microseconds_
    
        #avoid a lot of menial, repetitive calculations just for comparison purposes
        maxXminusHalfSize = self.maxX - self.config.TFsize/2.0
        minXplusHalfSize = -self.maxX + self.config.TFsize/2.0
        maxYminusHalfSize = self.maxY - self.config.TFsize/2.0
        minYplusHalfSize = -self.maxY + self.config.TFsize/2.0
        maxZminusHalfSize = self.maxZ - self.config.TFsize/2.0
        minZplusHalfSize = -self.maxZ + self.config.TFsize/2.0
    
        degreesOfFreedom = 3*len(tfs)
        randArr = np.random.normal(0,sigma,3*len(tfs))
    
        i=0
        tf=0
        while (i<degreesOfFreedom):
            if (tfs[tf].boundTo is None): #only currently unbound elements may move!
                tfs[tf].x += randArr[i]
                if (tfs[tf].x > maxXminusHalfSize): #prevent leaving the box!
                    tfs[tf].x = 2.0*maxXminusHalfSize - tfs[tf].x
                if (tfs[tf].x < minXplusHalfSize): #prevent leaving the box!
                    tfs[tf].x = 2.0*minXplusHalfSize - tfs[tf].x
                tfs[tf].y += randArr[i+1]
                if (tfs[tf].y > maxYminusHalfSize): #prevent leaving the box!
                    tfs[tf].y = 2.0*maxYminusHalfSize - tfs[tf].y
                if (tfs[tf].y < minYplusHalfSize): #prevent leaving the box!
                    tfs[tf].y = 2.0*minYplusHalfSize - tfs[tf].y
                tfs[tf].z += randArr[i+2]
                if (tfs[tf].z > maxZminusHalfSize): #prevent leaving the box!
                    tfs[tf].z = 2.0*maxZminusHalfSize - tfs[tf].z
                if (tfs[tf].z < minZplusHalfSize): #prevent leaving the box!
                    tfs[tf].z = 2.0*minZplusHalfSize - tfs[tf].z
            i+=3
            tf+=1
        return
    
    def printStatus(self, verbose=False):
        print "System simulation %.2f%% complete (t=%d) " % (((100.*self.time)/self.config.duration), self.time)
        if (verbose):
            print "\n TF Info:"
            for tf in self.tfs:
                print str(tf.x)+','+str(tf.y)+','+str(tf.z)+','+str(int(tf.boundTo is not None))
    
    def writeStatus(self, verbose=False):
        # output current systen status to a file -- need to think about output format and what data is useful
        # verbose is intended to mean write out the status of the entire system (as in all positions of
        # binding sites and tfs), but I'm not sure how to make that gnuplot-friendly...
        sysStats = self.calcGlobalStats()
        if (self.time==0):
            statsComment = '#timestamp(mu-s)'
        statsValues = str(self.time)
        statNames = sysStats.keys()
        statNames.sort()
        for statName in statNames:
            if (self.time==0):
                statsComment += ','+statName
            statsValues += '\t'+str(sysStats[statName])
        if (self.time==0):
            self.fileOut.write(statsComment+'\n')
        self.fileOut.write(statsValues+'\n')
    
    def calcGlobalStats(self):
        stats = {} #an associative array
        #Percentage of TFs Bound
        tfNum = len(self.tfs) #in case we allow them to ever be created or destroyed
        tfBound = 0
        tfx = 0
        tfy = 0
        tfz = 0
        for tf in self.tfs:
            tfx += tf.x #for center of mass of the TFs calculation
            tfy += tf.y
            tfz += tf.z
            if (tf.boundTo is not None):
                tfBound +=1
        stats['tfPercentBound'] = float(tfBound)/tfNum
        #Percentage of Binding Sites Occupied
        bNum = len(self.sites) #in case we allow them to ever be created or destroyed (less likely than with TFs)
        bBound = 0
        for site in self.sites:
            if (site.boundTo is not None):
                bBound +=1
        stats['bindPercentOccupied'] = float(bBound) / bNum
        #TF Center of Mass
        stats['tfmu-x'] = tfx / tfNum
        stats['tfmu-y'] = tfy / tfNum
        stats['tfmu-z'] = tfz / tfNum
        return stats
        #to output results, x = calcGlobalStats(); for i in x {print str(i)+':' + str(x[i])}
    
    def closeFiles(self):
        for tf in self.tfs:
            if (tf.file is not None):
                tf.file.close()
        for site in self.sites:
            if (site.file is not None):
                site.file.close()
        if (self.fileOut is not None):
            self.fileOut.close()



class SiteCluster:
    def __init__(self, size, posn=(0,0,0)):
        self.x=posn[0]
        self.y=posn[1]
        self.z=posn[2]
        self.size=size #as a radius
        self.sites=[] #the sites in the cluster

    def addSite(self, site):
        if (site is not None):
            self.sites.append(site)


class Config:
    def __init__(self):
        # This is the seed for the random number generator.  Each run's sequence of
        # random numbers is completely determined by this variable.  This enables us
        # to run multiple runs with the same sequence of random numbers.
        self.seed = 0 # any integer
        #numpy.random.randint(100000000) is a good way to get a new one (higher than 4 billion is unadvisable).
        
        # System-as-a-whole parameters
        ## Delta_t of each step through the simulation
        self.delta_t       = 1000 # (micro-seconds)
        ## Duration of the simulation
        self.duration      = 2000 * self.delta_t # (either multiply delta_t or set a raw time - in microseconds)
        #how often should we print out the current system status?
        self.printOutAt    = self.duration / (self.delta_t * 100) #(a raw number of timesteps into the simulation)
        #how often should we wite the system status to a file?
        self.writeSysOutAt = 100 #(a raw number of timesteps into the simulation)
        #Prefix of the names of the files to output to (the timestamp of the run will also be included as a suffix)
        self.fileOut       = 'Spring2012_Data/tf10b100-20min_clus0.8-r1'
        ## Dimensions of the system as a whole.  These values are in nm!
        self.systemWidth   = 500 # x-axis (nm)
        self.systemLength  = 500 # y-axis (nm)
        self.systemHeight  = 500 # z-axis (nm)
        
        
        # Various parameters relating to the Transcription Factor (TF) molecules.
        ## This value is the size of the TF molecules; if they are cubes, this is the
        ## side length, and if they are spheres, this is their diameter.
        self.TFsize       = 10/3.0 # (nm)
        ## The TF Diffusion Coefficient describes the scale of the random, Brownian
        ## movement exhibited by TF molecules.
        self.TFdiffusionC = 100 # (nm^2/s) this number and delta_t need to be set reasonably,
                             # so that collisions aren't very unlikely due to large displacements per timestep
                             # the mean displacement per timestep is (6*diff_C*delta_t)**0.5
        self.TFnumberInit = 100 # This is the number of TF molecules at system start
        
        
        # Binding Site Parameters
        ## This value is the size of binding sites; if they are cubes, this is the
        ## side length, and if they are spheres, this is their diameter.
        self.bindSize     = 10 # (nm)
        self.bindDistance = self.bindSize + self.TFsize
        self.bindDistanceSquared = self.bindDistance**2 #for expediency of the simulation
        ## Probability of binding for a TF which within the volume of the binding site
        ## during a single time-step of the system.
        self.pBind        = 1  # absolute probability [0,1]
        ## Probability of unbinding for a TF which within the volume of the binding
        ## site during a single time-step of the system.
        self.pUnbind      = (1.0 / (25.0 * (1.0 / (self.delta_t * 10**-6))))  # absolute probability [0,1] expressed this way, 25 seconds is the average bound time
        self.bindSitesNum = 100 # Number of binding sites in the system
        ## probability of each placed binding site to be within the cluster (if enabled)
        self.clusterProbability = 0.8 # a value of zero here turns off clustering behaviour completely
        ## side-length/diameter of the cluster region
        self.clusterSize = 3 * self.bindSize * self.bindSitesNum * self.clusterProbability # (nm)