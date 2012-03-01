# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 21:04:55 2012

@author: Ahaarnos
"""

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