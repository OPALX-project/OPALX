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
    outputFile.write('mpirun -np ' + str(parameters['_PROCX_'] * parameters['_PROCY_']) \
                     + ' impactt-1.5.1 ImpactT.in ' + str(parameters['rmsCalcStep']) + ' ' \
                     + str(parameters['phaseDumpStep']) + '\n')
    outputFile.write('kill `cat ~/.mpd.pid`\n')
    outputFile.write('cp ' + paths['runPath'] + '/ImpactT* ' + paths['storePath'])
    outputFile.close()
    
    # Run script/simulation.
    os.system('cd ' + paths['runPath'] + '; chmod +x dorun.sh')
    if runEnabled:
        os.system('cd ' + paths['runPath'] + '; ./dorun.sh')


def createSimFiles(paths, files, parameters, doEverything):
    # Routine writes simulation parameters to template file and writes gun input files.

    # Initial cycle through parameters. Make replacements in template file.
    for key, value in parameters.iteritems():
        if key == '_INIT_ENERGY_':
            # Calculate initial beta gamma of beam.
            kRestMass = 0.51099906 # Electron rest mass in MeV/c^2.
            gammaInitial = value / (1.0e6 * kRestMass) + 1
            betaGammaInitial = math.sqrt(math.pow(value / (1.0e6 * kRestMass) + 1.0, 2.0) - 1.0)
            print "Initial gamma: " + str(gammaInitial)
            print "Initial beta:  " + str(betaGammaInitial / gammaInitial)
            os.system('cd ' + paths['runPath'] + '; sed -i \'s/_INIT_BETAGAMMA_/' + str(betaGammaInitial) + '/\' ' + files['template'])
            os.system('cd ' + paths['runPath'] + '; sed -i \'s/' + key + '/' + str(value) + '/\' ' + files['template'])
        elif key == '_CURRENT_':
            os.system('cd ' + paths['runPath'] + '; sed -i \'s/' + key + '/' + str(value / parameters['_EBINS_']) + '/\' ' + files['template'])
        else:
            os.system('cd ' + paths['runPath'] + '; sed -i \'s/' + key + '/' + str(value) + '/\' ' + files['template'])

    # Calculate beam slice center offsets for energy binning.
    kCLight = 299792458.0 # Speed of light in vacuum (m/s).
    beamVelocity = kCLight * betaGammaInitial / gammaInitial
    beamLength = parameters['_TEMISSION_'] * beamVelocity
    deltaZ = beamLength / (2.0 * parameters['_EBINS_'])
    beamCenter = -deltaZ

    # Create input files.
    for index in range(parameters['_EBINS_']):
        beamCenter = -deltaZ - 2.0 * index * deltaZ
        if index == 0:
            fileName = 'ImpactT.in'
        else:
            fileName = 'ImpactT' + str(index + 1) + '.in'
        os.system('cd ' + paths['runPath'] + '; cp ' + files['template'] + ' ' + fileName)
        os.system('cd ' + paths['runPath'] + '; sed -i \'s/_SIGMA_/' + str(deltaZ) + '/\' ' + fileName)
        os.system('cd ' + paths['runPath'] + '; sed -i \'s/_BCENTER_/' + str(beamCenter) + '/\' ' + fileName)

    # Erase template file.
    os.system('cd ' + paths['runPath'] + '; rm -rf ' + files['template'])


def copyFiles(paths, files):    
    # Routine to copy and link necessary files for simulation.

    # Create run directory to store all of the input files and results. Use current time and date.
    currentTime = datetime.now() # Current time.
    runDirectory = currentTime.strftime("%d-%m-%Y-%H:%M:%S")
    os.system('cd ' + paths['runPath'] + '; mkdir ./' + runDirectory)
    paths['runPath'] += runDirectory

    # Copy template file to local directory.
    os.system('cp ' + paths['filePath'] + files['template'] + ' ' + paths['runPath'])

    # Create links to field files.
    os.system('cd ' + paths['runPath'] + '; ln -sf ' + paths['fieldPath'] + files['rfField'] + ' 1T1.T7')

        



        
    
