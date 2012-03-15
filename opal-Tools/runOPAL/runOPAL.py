#!/gpfs/homefelsim/adelmann/extlib/bin/python

"""
Script that launches simulation of the SwissFEL project 

@author: Andreas Adelmann <andreas.adelmann@psi.ch>
@author: Yves Ineichen <yves.ineichen@psi.ch>
@version: 0.1

export FIELDMAPS=/home5/schietinger/Fieldmaps/FIN/
export TEMPLATES=/home5/schietinger/Templates/FIN/
"""

import sys,re,os,string,fileinput,shutil,glob
import commands

from simulation import Simulation
from opaldict import OpalDict

def getNearestRestartStep(restart_pos,fn):
    res = commands.getoutput('H5getStep ' + str(restart_pos) + "  " + fn);
    return int(res.split("=")[2])

def getBaseName():
    path ='.'
    ext1 ='*.tmpl'
    ext2 ='*.data'
    str1 = glob.glob(os.path.join(path,ext1))[0]
    str2 = glob.glob(os.path.join(path,ext2))[0]
    str1spl = str1.split('.')
    str2spl = str2.split('.')

    if str1spl[1] == str2spl[1]:
        name = str2spl[1]
        name = name.split('/')[1]
    else:
        print 'template and data filename do not match'
        sys.exit()
    return name

def printUsage():
    print "./runOPAL.py [--restart-file=FILE [--restart-step=STEPNR | --restart-pos=POS]] [--test] [ATTR=SCANVALUE] {[ATTR=VALUE]}"
    print "Unit of POS is meter, SCANVALUE=start:end:step, example TFWHM=0.85:0.90:0.01 "
    print "--test does everything but submitting the job"

"""
Traverse all possible combinations of range variable values. Start simulation
once all range variables are fixed to a value. A list entry has the following
structure:
    ['name of var', start_value, end_value, step_value]
"""
def traverseRanges(list, opaldict, args):
    head = list[0]
    tail = list[1:]
    curval = head[1][0]
    while curval <= head[1][1]:
            opaldict[head[0]] = curval
            if len(tail) == 0:
                #run simultaion
                sim = Simulation(opaldict)
                sim.run(*args)
            else:
                traverseRanges(tail, opaldict, args)
            curval = curval + head[1][2]
 
"""
main method
"""
def main(argv):

    if (os.environ.get('TEMPLATES')):
        inputfilePath = os.environ.get('TEMPLATES')
	os.system('cp ' + inputfilePath + '/*.tmpl .')	
    else:
	ext1 ='*.tmpl'
        if (glob.glob(os.path.join('.',ext1))):
            inputfilePath = '../'
        else:
            print 'Template file unknown -> exiting ...'
            sys.exit()
            
    baseFileName = getBaseName()

    if os.environ.get('TEMPLATES'):
        if os.path.isfile(baseFileName + ".tmpl"):
	        os.system("mv " + baseFileName + ".tmpl tmplbak17.tmpl")	

    N = -1               # a running number; if given use it to label directory!
    doTest = -99
    restart_step = -99
    restart_file = ""

    dataFile = baseFileName + '.data'
    tmplFile = baseFileName + '.tmpl'
    oinpFile = baseFileName + '.in' # the resulting OPAL input file

    #create the dictionary
    opaldict = OpalDict(dataFile)
    # check if template values must be changed
    # if so add update the dictionary with the default values
    opaldict.addUserValues(argv)

    for arg in argv:
        if arg.startswith("--test"):
		    doTest=1
        elif arg.startswith("--restart-file"):
            restart_file = arg.partition("=")[2]
        elif arg.startswith("--restart-step"):
            restart_step = arg.split("=")[1]
        elif arg.startswith("--restart-pos"):
            restart_pos = arg.split("=")[1]
            restart_step = str(getNearestRestartStep(restart_pos,restart_file))
        elif arg.startswith("--help"):
            printUsage()
            exit()
        
    opaldict.scale()
    if not opaldict.hasRanges():
        sim = Simulation(opaldict)
        sim.run(N, baseFileName, restart_step, inputfilePath, tmplFile, oinpFile, restart_file, doTest)
    else:
        ranges = opaldict.Range()

        #create range toplevel dir
        dirname = baseFileName
        for p in opaldict.uservars:
            dirname += "_" + str(p[0]) + "=" + str(p[1])
        for (k, v) in ranges.items():
            dirname += "_" + k + "=" + str(v[0]) + ":" + str(v[1]) + ":" + str(v[2])
        # If there's already an file remove it...
        if os.path.isdir(dirname):
            print '\n\n'
            print 'REMOVE existing directory ', dirname
            shutil.rmtree(dirname)

        # create directory and change to the directory
        print dirname
        os.mkdir(dirname)
        os.chdir(dirname)
        
        #run simulations of all possible combinations
        args = [N, baseFileName, restart_step, inputfilePath, tmplFile, oinpFile, restart_file, doTest]
        traverseRanges(ranges.items(), opaldict, args)

        os.chdir("..")
        
    if os.path.isfile("tmplbak17.tmpl"):
	    os.system("mv tmplbak17.tmpl " + baseFileName + ".tmpl")	
    
#call main
if __name__ == "__main__":
    main(sys.argv[1:])
