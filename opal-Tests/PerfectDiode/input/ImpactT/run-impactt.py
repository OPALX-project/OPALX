#!/usr/bin/python
# Filename: run-impactt.py
#
# Purpose: Run an ImpactT simulation of an ideal pulsed diode gun for comparison with OPAL.
#
# Running ImpactT by just editing its input file(s) can be tedious, especially when you
# want to do many runs with different parameter settings. This Python script was written
# in an attempt to simplify the ordeal.
#
# To start a simulation, such as an electron gun, is a little bit complicated because
# to start the simulation properly we must break the beam into several energy bins
# (typically 32 or more). Each bin requires its own input file. These files are created using
# a script file written by Andreas Adelmann, a template file and a parameter file written by this script.
#
# Routine can also take arguments. If no arguments present they have no effect.
#
# Argument Number     Value: Description
# ===============     ===========
#        1            run/write:
#                                    run (or nothing) ==> Run simulation normally.
#                                    write            ==> Just write script files to run the problem, but don't start run.
import math, os, string, sys
import subroutines


# Parse command line arguments.
runEnabled = True
modInputFiles = False
writeLastBatchFile = False
if len(sys.argv) >= 2:
    if sys.argv[1] == 'run':
        runEnabled = True
    else:
        runEnabled = False
    if sys.argv[1] == 'write':
        writeLastBatchFile = True


# Main program.

#-------------------------------------------------------------------------#
# Define directories where simulations take place and results are stored. #
#-------------------------------------------------------------------------#
#
# Directories must already exist.

runPath = '/home2/russel/scratch/sims/my-OPAL-tests/ImpactT/'                                                           # Directory for running simulation.
storePath = '/home2/russel/my-OPAL-tests/test-sims/ideal-diode/uniform-dist/ImpactT/10-ebins/field-100MVperm/0.1-amps/' # Directory for storing results.


#-----------------------------------#
# Define parameters for simulation. #
#-----------------------------------#

filePath = './input-files/'            # Location of gun creation files.
fieldPath = '../../fields/'            # Directory containing field maps.
    
templateFile = 'input-file.tmpl'       # Template file for ImpactT input file.
rfFieldFile = 'ideal-diode-impactt.T7' # Electric field file.

numXProc = 1                    # Number of columns in processor layout. (Best to use power of 2.)
if numXProc > 2:
    # Make sure we don't use to many processors
    # on Merlin.
    numXProc = 2
numYProc = numXProc             # Number of rows in processor layout.
                                # (numXProx = numYProc as far as I know, or program crashes.)
                                # on the node are used. Be careful that you use an even number of nodes in this case.

dt = 1.0e-13                    # Integration time step (s).
numSteps = 1000                 # Number of time steps to take.
energyBins = 10                 # Number of energy bins to use in space charge calculation.

numPartPerBin = 5000            # Number of particles per energy bin (approximate, ImpactT may change this slightly).
imageOn = 1                     # Flag to turn image charge calculation on or off. (1 is on, 0 is off.)
imageCutoff = 0.02              # Distance from cathode where image charge calculation will be stopped (m).

numXMesh = 32                   # Number of space charge mesh points in x direction. (Best to use power of 2.)
numYMesh = 32                   # Number of space charge mesh points in y direction. (Best to use power of 2.)
numZMesh = 32                   # Number of space charge mesh points in z direction.
gridRx = 0.003                  # Radius of space charge grid domain in x direction (m).
gridRy = 0.003                  # Radius of space charge grid domain in y direction (m).
gridLz = 0.04                   # Length of space charge domain (m). This should be at least long enough to include all beam
                                # line elements you wish to simulate.
beamStop = 0.02                 # Location of beam stop (m). If greater than 0.0 simulation will stop simulation/beam
                                # at this location. Otherwise set equal to gridLz (above). (Note that particles that have
                                # gone beyond gridLz will not be written to the .h5 output file. This can be a problem
                                # if you allow beamStop to be equal to gridLz and are interested in the beam's behavior
                                # at the beamStop position.)

distributionType = 17           # Type of initial distribution to use.
emitSlices = 400                # Number of time steps to take during beam emission time (beam emission = pulseLength).
                                # The time step will be adjusted in ImpactT to make sure this happens.
pulseLength = 30.0e-12          # Total lengh of initial beam pulse (s).

rxBeam = 0.0005                 # Radius of beam in x direction (m).
ryBeam = 0.0005                 # Radius of beam in y direction (m).

current = 0.1                   # Total beam current (A). Equals bunch charge * reference frequency.
initialEnergy = 1.0             # Initial beam energy (eV).
referenceFreq = 1498.95e6       # Reference energy (Hz).

binMergeLocation = 0.0145       # Longitudinal location where ImpactT is to merge the separate beam bins into one bunch (m).

fieldMag = -100.0               # Field magnitude multiplier.
fieldLength = 0.02              # Length of diode field (m).
fieldFreq = 0.0                 # Frequency of field oscillation (Hz).
fieldRadius = 0.001             # Radius of field map (m).

rmsCalcStep = 5                 # Number of time steps between RMS beam property and slice emittance calculations.
phaseDumpStep = 5               # Number of time steps between phase dumps.


# Define dictionaries.
paths = dict(runPath = runPath,\
             storePath = storePath,\
             filePath = filePath,\
             fieldPath = fieldPath)
files = dict(template = templateFile,\
             rfField = rfFieldFile)
parameters = dict(_PROCX_ = numXProc,\
                  _PROCY_ = numYProc,\
                  _DT_ = dt,\
                  _NUMSTEPS_ = numSteps,\
                  _EBINS_ = energyBins,\
                  _NP_ = numPartPerBin,\
                  _IMAG_ON_ = imageOn,\
                  _IMAG_CUTOFF_ = imageCutoff,
                  _NX_ = numXMesh,\
                  _NY_ = numYMesh,\
                  _NZ_ = numZMesh,\
                  _GRIDRX_ = gridRx,\
                  _GRIDRY_ = gridRy,\
                  _GRIDLZ_ = gridLz,\
                  _BEAMSTOP_ = beamStop,\
                  _DIST_TYPE_ = distributionType,\
                  _SLICES_ = emitSlices,\
                  _TEMISSION_ = pulseLength,\
                  _RX_BEAM_ = rxBeam,\
                  _RY_BEAM_ = ryBeam,\
                  _CURRENT_ = current,\
                  _INIT_ENERGY_ = initialEnergy,\
                  _REF_FREQ_ = referenceFreq,\
                  _MERGE_LOCATION_ = binMergeLocation,\
                  _FIELD_MAG_ = fieldMag,\
                  _FIELD_LENGTH_ = fieldLength,\
                  _FIELD_FREQ_ = fieldFreq,\
                  _FIELD_RADIUS_ = fieldRadius,\
                  rmsCalcStep = rmsCalcStep,\
                  phaseDumpStep = phaseDumpStep)


#-----------------------------------------------------#
# Run simulation.                                     #
#-----------------------------------------------------#

# Create input files and links.

# Copy files for creating simulation to working directory.
subroutines.copyFiles(paths, files)

#Create input files.
subroutines.createSimFiles(paths, files, parameters, True)

# Run simulation. (Or just write scripts.)
subroutines.runMerlin(paths, files, parameters, runEnabled)



