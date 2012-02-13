# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 23:26:17 2012

@author: Ahaarnos
"""

class bindSite:
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
    
    def tfCollision(self, otherTFs, tfSize): #not used for now... maybe never
        sizeSquared = (2*tfSize)**2
        for site in otherTFs:
            if (((self.x-site.x)**2 + (self.y-site.y)**2 + (self.z-site.z)**2) < sizeSquared):
                print "Avoiding site collision in placement"
                return True
        return False