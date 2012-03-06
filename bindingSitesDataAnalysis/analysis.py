# -*- coding: utf-8 -*-
"""
Created on Wed Feb 29 20:24:28 2012

@author: Ahaarnos
"""

import os
import fnmatch
import scipy
import matplotlib.pyplot as plt

class RunSet:
    def __init__(self, folders):
        self.folders = folders # an array of folder names
        self.runs = [] # an array of run objects
        
        for folder in folders:
            self.runs.append(Run(folder))

class Run:
    def __init__(self, folder):
        self.sites = []
        self.folder = folder
        self.filenames = os.listdir(folder)
        self.systemFile = os.path.join(folder, fnmatch.filter(self.filenames, '*system*')[0])
        self.buildSitesArray()
        self.numSites = len(self.sites)
        
        self.duration = 0
        self.seed = 0
        self.parseSysFileForConfig()
        
    
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
        self.siteNum = -1
        self.totalBoundTime = 0
        
        self.boundStatus = {"t":[], "state":[]}
        
        self.previousState = [0, None, -1] # a line-to line state comparison
        self.parseSiteFile()
        
    def parseSiteFile(self):
        f = open(self.file)

        line = f.readline()
        while (line):
            if (line.find("#Site being tracked: #") != -1):
                xOfY = line[line.find("#Site being tracked: #")+len("#Site being tracked: #"):]
                self.siteNum = int(xOfY.split()[0])
            if (line[0] != "#" and line[0] != "s" and line[0] != '\n'): #the s case is a legacy case based on a data output error
                self._readSiteStateLine(line)
            line = f.readline()
        f.close()

    def _readSiteStateLine(self, line):
            params = line.split()
            if (params[1] == self.previousState[1]):
                # two True/Falses in a row... multibinding!   ...BAD!
                raise Exception("There is something malformed about "+str(self.siteNum))
            else:
                if (params[1] == "True" and self.previousState[1] == "False"):
                    #add a final 0 point to the array to show the ending of the unbound state
                    self.boundStatus["t"].append(int(params[0])-1)
                    self.boundStatus["state"].append(0)

                    #add the bound statu to the array
                    self.boundStatus["t"].append(int(params[0]))
                    self.boundStatus["state"].append(1)
                
                if (params[1] == "False" and self.previousState[1] == "True"):
                    delta_t = int(params[0]) - int(self.previousState[0])
                    self.bindDurations.append(delta_t)
                    self.totalBoundTime += delta_t
                    
                    #add a final 1 point to the array to show the ending of the bound state
                    self.boundStatus["t"].append(int(params[0])-1)
                    self.boundStatus["state"].append(1)

                    #add the unbound status to the array
                    self.boundStatus["t"].append(int(params[0]))
                    self.boundStatus["state"].append(0)
                    
            self.previousState = params
            


#Things to calculate/plot:
    # % of simulation time spent bound for each site
        # graph this distribution, identifying clustered and non-clustered sites
        # calculate the average as a whole, and the average for each group separately
    # plot a distribution of tf-site binding durations
        # calculate this average time and compare to you desired average used from the paper
    # avg number of bindings per site
        # total
        # in cluster
        # outside cluster
        
#Do all of the above per-ensemble, and for the sum of common ensembles

folders = ["../Spring2012_Data/tf100b100-30min_clus0.8_r1",
           "../Spring2012_Data/tf100b100-30min_clus0.8_r2",
           "../Spring2012_Data/tf100b100-30min_clus0.8_r3"]

print "Setting up the system"
system = RunSet(folders)

print "Beginning statistics calculations per-run"

print "Beginning stas calc for the whole system"

#draw the histogram
fig = plt.figure()
ax = fig.add_subplot(111)

data = []

for run in system.runs:
    for site in run.sites:
        data += site.bindDurations

ax.hist(data, 75, facecolor='green', alpha=0.75)

ax.set_xlabel("Binding Duration")

#plt.show()

total = 0.
for d in data:
    total += d
    
print total / len(data)


clusterPercentages = []
freePercentages = []

for run in system.runs:
    for site in run.sites:
        if site.inCluster:
            clusterPercentages.append((site.totalBoundTime + 0.0) / run.duration)
        else:
            freePercentages.append((site.totalBoundTime + 0.0) / run.duration)
            
totalPercentages = freePercentages + clusterPercentages

#draw the histogram
fig = plt.figure()
ax = fig.add_subplot(111)

data = []

for run in system.runs:
    for site in run.sites:
        data += site.bindDurations

ax.hist(clusterPercentages, 75, facecolor='red', alpha=0.75)

ax.set_xlabel("% Bound (of total Simulation)")


"""
for a list of folders corresponding to the same config:

/    get the system file for each run

/    parse its config data to know how many binding sites and the run duration, then and identify the run by the seed number

/    then load each free site file
    
/        name the site based on the original site number #Site being tracked: #blah (of 100) line

/        store it as an object in the sites array for the run
        
/        create an array of bind durations
        
/        determine total duration bound
        
            divide by the simulation duration to determine % of simulation time spent bound
        
/        create the array needed to graph the bound status
        
    then calculate per-run stats using all sites for the run
    
        get numbers for total, in-cluster, and outside-cluster for each stat
        
        be able to graph the above distributions for the run
        
    then calculate multi-run stats for all runs combined
    
        be able to graph the above distributions for the multi-run case
        
"""