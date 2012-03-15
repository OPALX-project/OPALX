# Filename: subroutines.py
#
# Purpose: Provide routines to modify input parameters and run simulation.
import datetime, math, os, string
from datetime import datetime

def runMerlin(paths, files, parameters, runEnabled):
    # Run simulation.

    # Write run script to file.
    outputFile = open(paths['runPath'] + '/dorun.sh', 'w')
    outputFile.write('#!/bin/bash\n')
    outputFile.write('source $HOME/.bashrc\n')
    outputFile.write('mpd --pid=$HOME/.mpd.pid&\n')
    outputFile.write('sleep 1\n')
    outputFile.write('mpiexec -np ' + str(parameters['numberOfProcessors']) \
                     + ' opal opal.in --commlib mpi --info 0 --warn 0 | tee opal.out\n')
    outputFile.write('kill `cat ~/.mpd.pid`\n')
    outputFile.write('cp ' + paths['runPath'] + '/opal* ' + paths['storePath'])
    outputFile.close()
    
    # Run script/simulation.
    os.system('cd ' + paths['runPath'] + '; chmod +x dorun.sh')
    if runEnabled:
        os.system('cd ' + paths['runPath'] + '; ./dorun.sh')


def copyFiles(paths, files):    
    # Routine to copy and link necessary files for simulation.

    # Create run directory to store all of the input files and results. Use current time and date.
    currentTime = datetime.now() # Current time.
    runDirectory = currentTime.strftime("%d-%m-%Y-%H:%M:%S")
    os.system('cd ' + paths['runPath'] + '; mkdir ./' + runDirectory)
    paths['runPath'] += runDirectory

    # Copy input file to local directory.
    os.system('cp ' + paths['filePath'] + files['inputFile'] + ' ' + paths['runPath'] + '/opal.in')

    # Create links to field files.
    os.system('cd ' + paths['runPath'] + '; ln -sf ' + paths['fieldPath'] + files['cavityFieldFile'])
