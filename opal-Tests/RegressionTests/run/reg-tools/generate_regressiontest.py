#!/usr/bin/python 
import datetime
import time
import commands
import os
import sys

"""
This method checks if all files in the reference directory are present
and if their md5 checksums still concure with the ones stored after
the simulation run
"""
def validateReferenceFiles(simname):
    allok = False

    if os.path.isfile("reference/" + simname + ".stat") and \
       os.path.isfile("reference/" + simname + ".stat.md5") and \
       os.path.isfile("reference/" + simname + ".out") and \
       os.path.isfile("reference/" + simname + ".out.md5") and \
       os.path.isfile("reference/" + simname + ".lbal") and \
       os.path.isfile("reference/" + simname + ".lbal.md5"):

        statout = commands.getoutput("md5sum --check reference/" + simname + ".stat.md5")
        outout = commands.getoutput("md5sum --check reference/" + simname + ".out.md5")
        lbalout = commands.getoutput("md5sum --check reference/" + simname + ".lbal.md5")

        print statout
        print outout
        print lbalout

        allok = statout == "reference/" + simname + ".stat: OK" and \
                outout == "reference/" + simname + ".out: OK" and \
                lbalout == "reference/" + simname + ".lbal: OK"

    return allok

"""
method to generate default rt file
"""
def writeDefaultRTFile(regtestname):
    rtfile = "\"" + regtestname + "\"\n"
    rtfile += "stat \"rms_x\" last 1E-15    #this is a comment\n"
    rtfile += "stat \"rms_y\" last 1E-15    #this is a comment\n"
    rtfile += "stat \"rms_s\" last 1E-15    #this is a comment"

    rtoutfile = open(regtestname + ".rt", "w")
    rtoutfile.writelines(rtfile)
    rtoutfile.close()

"""
method to generate sge file
"""
def writeSgeFile(regtestname):

    sgefile = "#!/bin/bash \n"
    sgefile += "#$ -cwd\n"
    sgefile += "#$ -j y\n"
    sgefile += "#$ -pe mpi 4\n"
    sgefile += "#$ -N " + regtestname + "-RT\n"
    sgefile += "#$ -v MPIHOME=/opt/mpi/openmpi-1.2.6-intel-10.0,LD_LIBRARY_PATH=/opt/intel-mkl/mkl-10.0/lib/em64t:/opt/mpi/openmpi-1.2.6-intel-10.0/lib:/opt/intel/intel-10.0/fce-10.0/lib:/opt/intel/intel-10.0/cce-10.0/lib,OPAL_EXE_PATH,REG_TEST_DIR\n"
    sgefile += "\n"
    sgefile += "MACHINE_FILE=$TMPDIR/machinefile\n"
    sgefile += "awk '/^felsim/ {print $1\" slots=\"$2}' $PE_HOSTFILE > $MACHINE_FILE\n"
    sgefile += "cp $MACHINE_FILE machinefile.last\n"
    sgefile += "\n"
    sgefile += "echo \"PE_HOSTFILE:\"\n"
    sgefile += "cat  $PE_HOSTFILE\n"
    sgefile += "echo \"MACHINE_FILE:\"\n"
    sgefile += "cat $MACHINE_FILE\n"
    sgefile += "echo \"SLOTS=$NSLOTS\"\n"
    sgefile += "\n"
    sgefile += "cd $REG_TEST_DIR\n"
    sgefile += "OPAL=\"$OPAL_EXE_PATH/opal " + regtestname + ".in --commlib mpi --info 0 --warn 0 2>&1\"\n"
    sgefile += "CMD=\"$MPIHOME/bin/mpirun -machinefile $MACHINE_FILE -np $NSLOTS  --mca ras localhost --mca pls rsh  $OPAL \"\n"
    sgefile += "$CMD"
    
    sgeoutfile = open(regtestname + ".sge", "w")
    sgeoutfile.writelines(sgefile)
    sgeoutfile.close()

"""
main method
assuming we have a running simulation in the current directory (T7's and input file)
FIXME:
    o currently cannot handle restarted simulations
    o currently just copies all T7, dat und phases files to reg test dir
"""
def main(argv):
    d = datetime.date.today()
    #chdir to path of script
    #os.chdir(sys.path[0])

    regtestname = ""
    if ".in" in argv[0]:
        #extract input file
        regtestname = (str.split(argv[0], "."))[0]
    else:
        print "./generate-regressiontest.py RegressionTestName.in"
        return

    #check if user has already set an OPAL executable
    #if not use the one from the last build test
    env = os.getenv("OPAL_EXE_PATH")
    if env is None:
        print "OPAL_EXE_PATH not set! Please set it to the directory containing the opal binary and rerun the script."
        return

    print "Starting to generate Regression Test: " + regtestname

    #only generate regression tests if opal executable is valid
    if os.path.isfile(os.getenv("OPAL_EXE_PATH") + "/opal"):
        commands.getoutput("mkdir " + regtestname)
        commands.getoutput("cp " + regtestname + ".in " + regtestname + "/")
        commands.getoutput("cp *.T7 " + regtestname + "/")
        commands.getoutput("cp *.dat " + regtestname + "/")
        commands.getoutput("cp *.phases " + regtestname + "/")
        os.chdir(regtestname)
        writeSgeFile(regtestname)
        os.environ["REG_TEST_DIR"] = os.getcwd()
        runout = commands.getoutput("qsub " + regtestname + ".sge -v REG_TEST_DIR=" + os.getcwd() + ",OPAL_EXE_PATH=" + os.getenv("OPAL_EXE_PATH"))
        username = commands.getoutput("whoami")
        #wait for job to finish for now do this serially
        qstatout = commands.getoutput("qstat -u " + username)
        while len(qstatout) > 0:
            #we only check every minute if job has finished
            time.sleep(10)
            qstatout = commands.getoutput("qstat -u " + username)

        #copy o to out file
        commands.getoutput("mv " + regtestname + "-RT.o* " + regtestname + ".out")
        commands.getoutput("rm " + regtestname + "-RT.po*")
        commands.getoutput("mkdir reference")
        commands.getoutput("mv " + regtestname + ".out reference/")
        commands.getoutput("mv " + regtestname + ".stat reference/")
        commands.getoutput("mv " + regtestname + ".lbal reference/")
        commands.getoutput("md5sum reference/" + regtestname + ".out > reference/" + regtestname + ".out.md5")
        commands.getoutput("md5sum reference/" + regtestname + ".lbal > reference/" + regtestname + ".lbal.md5")
        commands.getoutput("md5sum reference/" + regtestname + ".stat > reference/" + regtestname + ".stat.md5")
        if validateReferenceFiles(regtestname):
            print "Reference data verified successfully"
        else:
            print "Reference data verification failed!"
            return
        writeDefaultRTFile(regtestname)

    else:
        print "OPAL_EXE_PATH not set! Please set it to the directory containing the opal binary and rerun the script."
        return
    
    print "Finished generating Regression Test " + regtestname
    commands.getoutput("rm machinefile.last")
    commands.getoutput("rm " + regtestname + ".h5")
    print "Please edit the " + regtestname + ".rt file to suit your needs and cleanup the regression test directory. Once this is done move it to opal-Tests/RegressionTests/, svn add and commit."

#call main
if __name__ == "__main__":
    main(sys.argv[1:])
