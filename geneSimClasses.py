# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 23:26:17 2012

@author: Ahaarnos
"""

class BindSite:
    def _init_(self, posn=(0,0,0), boundTo=None, inCluster=False):
        self.x=posn[0]
        self.y=posn[1]
        self.z=posn[2]
        self.boundTo=boundTo
        self.inCluster=inCluster
    
    def siteOverlap(self, otherSites, siteSize):
        sizeSquared = (2*siteSize)**2
        for site in otherSites:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < sizeSquared):
                print "Avoiding site collision in placement"
                return True
        return False

class TF:
    def _init_(self, posn=(0,0,0), boundTo=None):
        self.x=posn[0]
        self.y=posn[1]
        self.z=posn[2]
        self.prevX=posn[0]
        self.prevY=posn[1]
        self.prevZ=posn[2]
        self.boundTo=boundTo
        
    def _str_(self):
        return str(self.x)+'\t'+str(self.y)+'\t'+str(self.z)+'\t'+'\t'+str(self.boundTo)
    
    def tfCollision(self, otherTFs, tfSize): #not used for now... maybe never
        sizeSquared = (2*tfSize)**2
        for site in otherTFs:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < sizeSquared):
                print "Avoiding site collision in placement"
                return True
        return False

class Simulation:
    def _init(self):
        self.time=0
        self.dt=0 #may not need to be here
        self.maxX=0
        self.maxY=0
        self.maxZ=0

class SiteCluster:
    def _init(self):
        self.x=0
        self.y=0
        self.z=0
        self.halfWidth=0
        self.sites=[] #the sites in the cluster