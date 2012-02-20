#Perform some initial setup
#import numpy
#import scipy
#import pylab
#import time

import geneSimConfig as c
import geneSimClasses as gs

if (__name__ == "__main__"):
    print "Initializing the system"
    system = gs.Simulation(c.config)
    print "Beginning time evolution..."
    system.runSimulation()
    system.closeFiles()