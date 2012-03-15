#!/usr/bin/python
# Filename: run-opal.py
#
# Purpose: Run an OPAL simulation.
#
# Argument Number     Value: Description
# ===============     ===========
#        1            run/write:
#                                    run (or nothing) ==> Run simulation normally.
#                                    write            ==> Just write script files to run the problem, but don't start run.
#
import math, os, string, sys
import subroutines

# Parse command line arguments.
runEnabled = True
if len(sys.argv) >= 2:
    if sys.argv[1] == 'run':
        runEnabled = True
    else:
        runEnabled = False
        

# Main program.

#-------------------------------------------------------------------------#
# Define directories where simulations take place and results are stored. #
#-------------------------------------------------------------------------#
#
# Directories must already exist.

runPath = '/home2/russel/scratch/sims/my-OPAL-tests/OPAL/' # Run directory.
storePath = '/home2/russel/my-OPAL-tests/test-sims/ideal-diode/uniform-dist/OPAL/10-ebins/field-100MVperm/0.1-amps/' # Directory for storing results.

#-----------------------------------#
# Define parameters for simulation. #
#-----------------------------------#

filePath = './'                    # Location of OPAL input file.
fieldPath = '../../fields/'        # Directory containing field maps.
    
inputFile = 'input-OPAL.in'        # Input file for OPAL.
cavityFieldFile = 'ideal-diode.T7' # Electric field file.

numberOfProcessors = 4             # Number of processors to be used.

# Define dictionaries.
paths = dict(runPath = runPath,\
             storePath = storePath,\
             filePath = filePath,\
             fieldPath = fieldPath)
files = dict(inputFile = inputFile,\
             cavityFieldFile = cavityFieldFile)
parameters = dict(numberOfProcessors = numberOfProcessors)

#-----------------#
# Run Simulation. #
#-----------------#

# Create input files and links.

# Copy files to working directory.
subroutines.copyFiles(paths, files)

# Run simulation. (Or just write scripts.)
subroutines.runMerlin(paths, files, parameters, runEnabled)
