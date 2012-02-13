"""
Set various configuration variables for each run here.  All values need to be 
rigourously documented!
"""

# This is the seed for the random number generator.  Each run's sequence of
# random numbers is completely determined by this variable.  This enables us
# to run multiple runs with the same sequence of random numbers.
seed = 66143 # any integer
#numpy.random.randint(100000000) is a good way to get a new one (higher than 4 billion is unadvisable).

# System-as-a-whole parameters
## Delta_t of each step through the simulation
delta_t       = 1000 # (micro-seconds)
## Duration of the simulation
duration      = 600000 * delta_t # (either multiply delta_t or set a raw time - in microseconds)
#how often should we print out the current system status?
printOutAt    = duration / (delta_t * 100) #(a raw number of timesteps into the simulation)
#how often should we wite the system status to a file?
writeSysOutAt = 100 #(a raw number of timesteps into the simulation)
#Prefix of the names of the files to output to (the timestamp of the run will also be included as a suffix)
fileOut       = 'Spring2012_Data/tf10b100-20min_clus0.8-r1'
## Dimensions of the system as a whole.  These values are in nm!
systemWidth   = 1000 # x-axis (nm)
systemLength  = 1000 # y-axis (nm)
systemHeight  = 1000 # z-axis (nm)


# Various parameters relating to the Transcription Factor (TF) molecules.
## This value is the size of the TF molecules; if they are cubes, this is the
## side length, and if they are spheres, this is their diameter.
TFsize       = 10/3.0 # (nm)
## The TF Diffusion Coefficient describes the scale of the random, Brownian
## movement exhibited by TF molecules.
TFdiffusionC = 100 # (nm^2/s) this number and delta_t need to be set reasonably,
                     # so that collisions aren't very unlikely due to large displacements per timestep
                     # the mean displacement per timestep is (6*diff_C*delta_t)**0.5
TFnumberInit = 1000 # This is the number of TF molecules at system start


# Binding Site Parameters
## This value is the size of binding sites; if they are cubes, this is the
## side length, and if they are spheres, this is their diameter.
bindSize     = 10 # (nm)
bindDistanceSquared = (bindSize+TFsize)**2 #for expediency of the simulation
## Probability of binding for a TF which within the volume of the binding site
## during a single time-step of the system.
pBind        = 1  # absolute probability [0,1]
## Probability of unbinding for a TF which within the volume of the binding
## site during a single time-step of the system.
pUnbind      = (1.0 / (25.0 * (1.0 / (delta_t * 10**-6))))  # absolute probability [0,1] expressed this way, 25 seconds is the average bound time
bindSitesNum = 100 # Number of binding sites in the system
clustering   = True # whether or not binding sites should apppear in clusters
## probability of each placed binding site to be within the cluster (if enabled)
clusterProbability = 0.8
## side-length/diameter of the cluster region
clusterSize = 3 * bindSize * bindSitesNum * clusterProbability # (nm)