# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 23:05:01 2012

@author: Ahaarnos
"""

import geneSimClasses as gs

"""
Set various configuration variables for each run here.  All values need to be 
rigourously documented!
"""

config = gs.Config();

# This is the seed for the random number generator.  Each run's sequence of
# random numbers is completely determined by this variable.  This enables us
# to run multiple runs with the same sequence of random numbers.
config.seed = 6483 # any integer
#numpy.random.randint(100000000) is a good way to get a new one (higher than 4 billion is unadvisable).

# System-as-a-whole parameters
## Delta_t of each step through the simulation
config.delta_t       = 1000 # (micro-seconds)
## Duration of the simulation
config.duration      = 1800000 * config.delta_t # (either multiply delta_t or set a raw time - in microseconds)
#how often should we print out the current system status?
config.printOutAt    = config.duration / (config.delta_t * 100) #(a raw number of timesteps into the simulation)
#how often should we wite the system status to a file?
config.writeSysOutAt = 500 #(a raw number of timesteps into the simulation)
## Dimensions of the system as a whole.  These values are in nm!
config.systemWidth   = 500 # x-axis (nm)
config.systemLength  = 500 # y-axis (nm)
config.systemHeight  = 500 # z-axis (nm)


# Various parameters relating to the Transcription Factor (TF) molecules.
## This value is the size of the TF molecules; if they are cubes, this is the
## side length, and if they are spheres, this is their diameter.
config.TFsize       = 10/3.0 # (nm)
## The TF Diffusion Coefficient describes the scale of the random, Brownian
## movement exhibited by TF molecules.
config.TFdiffusionC = 100 # (nm^2/s) this number and delta_t need to be set reasonably,
                     # so that collisions aren't very unlikely due to large displacements per timestep
                     # the mean displacement per timestep is (6*diff_C*delta_t)**0.5
config.TFnumberInit = 10 # This is the number of TF molecules at system start


# Binding Site Parameters
## This value is the size of binding sites; if they are cubes, this is the
## side length, and if they are spheres, this is their diameter.
config.bindSize     = 10 # (nm)
config.bindDistance = config.bindSize + config.TFsize
config.bindDistanceSquared = config.bindDistance**2 #for expediency of the simulation
config.trackAllSites = True
## Probability of binding for a TF which within the volume of the binding site
## during a single time-step of the system.
config.pBind        = 1  # absolute probability [0,1]
## Probability of unbinding for a TF which within the volume of the binding
## site during a single time-step of the system.
config.pUnbind      = (1.0 / (25.0 * (1.0 / (config.delta_t * 10**-6))))  # absolute probability [0,1] expressed this way, 25 seconds is the average bound time
config.bindSitesNum = 100 # Number of binding sites in the system
## probability of each placed binding site to be within the cluster (if enabled)
config.clusterProbability = 0.8 # a value of zero here turns off clustering behaviour completely
## side-length/diameter of the cluster region
config.clusterSize = 3 * config.bindSize * (config.bindSitesNum)**(1/3.) * config.clusterProbability # (nm)

#Prefix of the names of the files to output to (the timestamp of the run will also be included as a suffix)
config.fileOut       = 'Spring2012_Data/tf10b100-30min_clus0.8-r1'