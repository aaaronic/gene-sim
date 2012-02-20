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
        
    def __str__(self):
        return str(self.x)+'\t'+str(self.y)+'\t'+str(self.z)+'\t'+'\t'+str(self.boundTo)
    
    def tfCollision(self, otherTFs, tfSize): #not used for now... maybe never
        sizeSquared = (2*tfSize)**2
        for site in otherTFs:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < sizeSquared):
                print "Avoiding site collision in placement"
                return True
        return False

class Simulation:
    def __init__(self, config):
        np.random.seed(config.seed)
        
        self.t=0
        self.config=config
        self.tfs=[]
        self.sites=[]
        self.cluster=None
        
        self.start_t = time.localtime() #all files will share the same timestamp as part of their name
        self.timestamp = str(self.start_t.tm_mon)+'-'+str(self.start_t.tm_mday)+'-'+str(self.start_t.tm_year)+'-h'+str(self.start_t.tm_hour)+'m'+str(self.start_t.tm_min)+'s'+str(self.start_t.tm_sec)
        self.fileOut = open(config.fileOut+'-system-'+self.timestamp,'w',0)

        maxX = config.systemWidth/2
        maxY = config.systemLength/2
        maxZ = config.systemHeight/2
        
        for i in range(config.TFnumberInit):
            xPosn = np.random.uniform(-maxX + config.TFsize/2.0,maxX-config.TFsize/2.0)
            yPosn = np.random.uniform(-maxY + config.TFsize/2.0,maxY-config.TFsize/2.0)
            zPosn = np.random.uniform(-maxZ + config.TFsize/2.0,maxZ-config.TFsize/2.0)
            self.tfs.append(TF(i, config.TFsize, (xPosn,yPosn,zPosn)))

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
    
    def writeConfig(self, f):
        f.write('### Configuration used to generate this run ###\n')
        var = dir(self.config)
        var.sort(key=str.upper)
        for v in var:
            if (v[0] != '_'):
                f.write('# ' + v + ': ' + str(eval('self.config.'+str(v))) + '\n') #output the name of v, then the value of v
        f.write('\n\n')










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