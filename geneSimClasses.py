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

    def siteOverlap(self, otherSites):
        for site in otherSites:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < (self.size+site.size)**2):
                print "Avoiding site collision in placement"
                return True
        return False

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
            if (s.boundTo is not None): # ensure only unbound sites are considered for binding
                distSquared = (s.x - self.x)**2
                if (distSquared < bindDistanceSquared): #check each coordinate individually before bothering to continue to the full comparison
                    distSquared += (s.y - self.y)**2 #dy**2
                    if (distSquared < bindDistanceSquared):
                        distSquared += (s.z - self.z)**2
                        if (distSquared <= bindDistanceSquared):
                            return s
        return None #no collision

class Simulation:
    def __init__(self, config):
        np.random.seed(config.seed)

        self.t=0
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
        maxY = config.systemLength/2
        maxZ = config.systemHeight/2

        for i in range(config.TFnumberInit):
            xPosn = np.random.uniform(-maxX + config.TFsize/2.0,maxX-config.TFsize/2.0)
            yPosn = np.random.uniform(-maxY + config.TFsize/2.0,maxY-config.TFsize/2.0)
            zPosn = np.random.uniform(-maxZ + config.TFsize/2.0,maxZ-config.TFsize/2.0)
            self.tfs.append(TF(i, config.TFsize, (xPosn,yPosn,zPosn)))

        tfToTrack = np.random.randint(len(self.tfs))
        self.tfs[tfToTrack].file = open(self.config.fileOut+'-TF-'+self.timestamp,'w',0)
        self.writeConfig(self.tfs[tfToTrack].file)
        self.writeTFStatus(0,self.tfs[tfToTrack])

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
            self.writeSiteStatus(0, self.sites[siteToTrack])
        else:
            clusterX = np.random.uniform(-maxX + config.clusterSize/2.0,maxX-config.clusterSize/2.0)
            clusterY = np.random.uniform(-maxY + config.clusterSize/2.0,maxY-config.clusterSize/2.0)
            clusterZ = np.random.uniform(-maxZ + config.clusterSize/2.0,maxZ-config.clusterSize/2.0)
            self.cluster = SiteCluster(config.clusterSize, (clusterX, clusterY, clusterZ))
            print "Running With with a clusterprobability of " + str(config.clusterProbability)
            print "Cluster centered at:\n" + str(self.clusterPosn)
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
            self.writeSiteStatus(0, self.sites[siteToTrack])

            clusterSiteToTrack = self.cluster.sites[np.random.randint(len(self.cluster.sites))].siteNum
            self.sites[clusterSiteToTrack].file = open(self.config.fileOut+'-bindingSite-'+str(siteToTrack)+'-inCluster-'+self.timestamp,'w',0)
            self.writeConfig(self.sites[clusterSiteToTrack].file)
            self.writeSiteStatus(0, self.sites[clusterSiteToTrack])
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
        diffuseTFs(unBoundTFs,delta_t,TFdiffusionC)
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

    def writeTFStatus(self, time,  tf):
        if (tf.file is not None): # make sure there's something to write to before continuing
            if (time==0):
                tfComment = '#TF being tracked: #'+str(tf.tfNum)+' (of '+str(len(self.tfs))+')\n#timestamp(mu-s),x,y,z,status(1 is bound),boundTo(None is unbound)'
            tfVals = str(time)+'\t'+str(tf.x)+'\t'+str(tf.y)+'\t'+str(tf.z)+'\t'+str(int(tf.boundTo is not None))+'\t'+str(tf.boundTo.siteNum)
            if (time==0):
                tf.file.write(tfComment+'\n')
            tf.file.write(tfVals+'\n')

    def writeSiteStatus(self, time, site):
        if (site.file is not None): # make sure there's something to write to before continuing
            if (time==0):
                siteComment = '#Site being tracked: #'+str(site.siteNum)+' (of '+str(len(self.sites))+')\n'
                siteComment += '#located at position: '+str(site.x)+'\t'+str(site.y)+'\t'+str(site.z)+'\n'
                if site.inCluster:
                    siteComment += 'site is within a cluster'
                siteComment += '#timestamp(mu-s),status(0/1 for unbound/bound),boundTo(None means unbound!)'
                site.file.write(siteComment+'\n')

            siteVals = str(time)+'\t'+str(site.boundTo is not None)+'\t'+str(site.boundTo.tfNum)
            site.file.write(siteVals+'\n')










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