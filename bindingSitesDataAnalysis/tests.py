# -*- coding: utf-8 -*-
"""
Created on Wed Feb 29 21:22:58 2012

@author: Aaron
"""

import os
import sys
import fnmatch

x = os.listdir("../Spring2012_Data/tf10b100-30min_clus0.8_r1")

systemFile = fnmatch.filter(x, '*system*')[0]

class Run:
    def __init__(self, folder):
        self.sites = []
        self.folder = folder
        self.filenames = os.listdir(folder)
        self.systemfile = os.path.join(folder, fnmatch.filter(self.filenames, '*system*')[0])
        self.buildSitesArray()
        
        self.duration = 0
        self.seed = 0
        self.parseSysFileForConfig()
        
        
        self.numSites = len(self.sites)
    
    def buildSitesArray(self):
        for fname in fnmatch.filter(self.filenames, "*ClusterSite*"):
            self.sites.append(Site(os.path.join(self.folder, fname), True))
            
        for fname in fnmatch.filter(self.filenames, "*FreeSite*"):
            self.sites.append(Site(os.path.join(self.folder, fname), False))
            
    def parseSysFileForConfig(self):
        f = open(self.systemFile)

        line = f.readline()
        while (len(line) > 3):  #stop reading after we've read in all of the config
            if (line.find("duration: ") != -1):
                self.duration = float(line[line.find("duration: ") + len("duration: "):])
            if (line.find("seed: ") != -1):
                self.seed = float(line[line.find("seed: ") + len("seed: "):])
            line = f.readline()
        f.close()
        
        
class Site:
    def __init__(self, filePath, inCluster):
        self.file = filePath
        self.inCluster = inCluster
        self.bindDurations = []
