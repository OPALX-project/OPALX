"""
Simulation class handles sge job related things

@author: Andreas Adelmann <andreas.adelmann@psi.ch>
@author: Yves Ineichen <yves.ineichen@psi.ch>
@version: 0.1
"""

import sys,re,os,string,fileinput,shutil,glob,commands

class Simulation:

    def __init__(self, opaldict):
        self.opaldict = opaldict
        self.dirname = ""

    def isOldMpiVersion(self):
        mpistr = os.environ.get('OPENMPI')
        mpiver = mpistr.partition("-")[2].partition('-')[0]
        mpiv = mpiver.split(".")

        if mpiv[0] == "1" and mpiv[1] >= "4":
            return False
        else:
            return True

    def createDirectory(self, dirname):
        # If there's already an file remove it...
        if os.path.isdir(self.dirname):
            print '\n\n'
            print 'REMOVE existing directory ', self.dirname
            shutil.rmtree(self.dirname)

        # create directory and change to the directory
        print self.dirname
        os.mkdir(self.dirname)

    def run(self,N, baseFileName, restart_step, inputfilePath, tmplFile, oinpFile, restart_file, doTest):
        # make directory name indicating chanched values
        self.dirname = baseFileName
        if N >= 0:
            self.dirname += str(N)

        self.dirname += self.opaldict.generateDirectoryName()

        if restart_step != -99:
            self.dirname += "_RESTART-STEP=" + restart_step

        CORES = int(self.opaldict['CORES'])

        self.createDirectory(self.dirname)
        os.chdir(self.dirname)

        if (os.environ.get('FIELDMAPS')):
            fieldmapPath = os.environ.get('FIELDMAPS')
        else:
            fieldmapPath = '../fieldmaps'
            if not (os.path.isdir(fieldmapPath)):
                print 'Fieldmap directory unknown exiting ...'
                sys.exit()
            
        # linking magnet and rf files
        os.system('lndir ' + fieldmapPath + ' ')

        # os.system('ln -fs /gpfs/homefelsim/amas/Simulations/Misc/ReferenceParticles.dat .')
        os.system('cp ' + inputfilePath + tmplFile + ' ' +  oinpFile)

        # do the replacements in the templatefile
        for s in self.opaldict:
            os.system('sed -i \'s/_'+ s +'_/'+str(self.opaldict[s])+'/\' '+oinpFile)

        #opalexe = os.environ.get('OPAL_EXE_PATH')
        opalexe = '$OPAL_EXE_PATH/opal'

        print 'Simulation directory is ' + self.dirname + ' using OPAL at ', os.environ.get('OPAL_EXE_PATH')
        print 'Using templatefile at ' + inputfilePath + ' using fieldmaps at ' + fieldmapPath + ' \n'
        print 'Parameter set in ' + oinpFile + ' are:'

        for s in self.opaldict:
            print ' :::: ' + s + '= ' + str(self.opaldict[s])

        self.WriteSGE(opalexe, oinpFile, restart_step, restart_file)
       
        if restart_step != -99 and restart_file != "":
            # this should be fixed in H5merge
            step = int(restart_step) - 1
            commands.getoutput('/gpfs/homefelsim/amas/bin/H5merge ' + restart_file + "[0:" + str(step) + "]" + " " + baseFileName + ".h5")
        
        #check if queue is specified in environment variable
        queue = os.getenv("SGE_QUEUE", "all.q")
        if doTest==-99:
            os.system('qsub -q ' + queue + ' -pe orte '+str(CORES)+' run.sge')
            print 'Done with setup of the OPAL simulation and submitting the job to ',queue, ' with ',CORES,'cores \n\n\n'
        else:
            print 'Done with setup of the OPAL simulation but not submitting the job (--test) \n\n\n'
        os.chdir('..')

    def WriteSGE(self, opalexe, oinpFile, restart_step, restart_fn):
        title=oinpFile.partition(".")[0]
        myfile = open('run.sge','w')
        s1 = "#!/bin/bash\n"
        s1 += "#$ -cwd\n"
        s1 += "#$ -j y\n"
        s1 += "#$ -pe orte 16\n"
        s1 += "#$ -N " + title + "\n"
        s1 += "#$ -v LD_LIBRARY_PATH,OPAL_EXE_PATH,OPENMPI\n\n"
        s1 += "MACHINE_FILE=$TMPDIR/machinefile\n"
        s1 += "awk '/^felsim/ {print $1\" slots=\"$2}' $PE_HOSTFILE > $MACHINE_FILE\n"
        s1 += "cp $MACHINE_FILE .machinefile.last\n\n"
        if restart_step != -99:
            s1 += "OPAL=\"" + opalexe + " " + oinpFile + " -restartfn " + restart_fn +   " -restart -1  --commlib mpi --info 0 --warn 0 | tee .out \"\n"
        else:
            s1 += "OPAL=\"" + opalexe + " " + oinpFile + " --commlib mpi --info 0 --warn 0 \"\n"
        s1 += "CMD=\"$OPENMPI/bin/mpirun -x LD_LIBRARY_PATH -machinefile $MACHINE_FILE -np $NSLOTS  $OPAL \"\n"
        s1 += "echo \"Running $CMD\"\n"
        s1 += "echo\n"
        s1 += "$CMD"

        myfile.write(s1)
        myfile.close()
