#!/usr/bin/python 
import datetime
import commands
import os
import sys

from reporter import Reporter
from reporter import TempXMLElement

from builder import Builder

from regressiontest import RegressionTest

from documentation import OpalDocumentation
from documentation import OpalDoxygen

from tools import readfile
from tools import module_load
from tools import sendmails

#FIXME: ugly global variables
totalNrTests = 0
totalNrPassed = 0

"""
This method traverses the directory tree. It will check and execute regression tests for the following directory-layouts:

    DIR Structure:
    name/name.in
         name.rt
         run
         run-parallel
         *.T7
         reference/name.lbal
         reference/name.out
         reference/name.stat
         reference/name.lbal.md5
         reference/name.out.md5
         reference/name.stat.md5

Please make sure you use this naming scheme!
"""
def callback(arg, dirname, fnames):

    global totalNrPassed
    global totalNrTests

    if not ".svn" in dirname and not "reference" in dirname: #exclude svn and reference dirs

        dir = str.split(dirname, "/")
        simname = dir[len(dir)-1]

        # check if all files required are available
        if os.path.isfile(dirname + "/" + simname + ".in") and \
           os.path.isfile(dirname + "/" + simname + ".rt") and \
           os.path.isfile(dirname + "/" + simname + ".sge") and \
           os.path.isdir(dirname + "/" + "reference"):

                rep = Reporter()
                rep.appendReport("Found valid test in %s \n" % dirname)

                # check if we really want to run this test
                # empty list = run all
                runtests = arg[0]
                if runtests:
                    if not simname in runtests:
                        rep.appendReport("User decided to skip regression test %s \n" % simname)
                        return

                d = datetime.date.today()
                resultdir = "results/" + d.isoformat() + "/" + simname
                if not os.path.isdir(resultdir):
                    commands.getoutput("mkdir -p " + resultdir)
                
                simulation_report = TempXMLElement("Simulation")
                simulation_report.addAttribute("name", simname)
                simulation_report.addAttribute("date", "%s" % d)

                rt = RegressionTest(dirname, simname, resultdir)
                rt.run(simulation_report)
                totalNrTests += rt.totalNrTests
                totalNrPassed += rt.totalNrPassed

def bailout(runAsUser):
    rep = Reporter()
    d = datetime.date.today()
    rep.appendReport("==========================================================\n")
    rep.appendReport("Finished Regression Test on %s \n" % datetime.datetime.today())
    rep.appendReport("\n")
    rep.appendReport("http://amas.web.psi.ch/regressiontests/results_%s_%s_%s.xml \n\n" % (d.day, d.month, d.year))

    #send/print report
    if not runAsUser:
        os.chdir(sys.path[0])
        emails = readfile("email-list")
        sendmails(emails, rep.getReport(), totalNrTests)
    else:
        print rep.getReport()
    

def main(argv):
    rep = Reporter()
    #various paths needed
    www_folder = "/afs/psi.ch/project/amas/www/regressiontests"
    os.chdir(sys.path[0]) #chdir to path of script
    rundir = os.getcwd()
    regdir = (str.split(rundir, "/run"))[0]
    srcdir = (str.split(rundir, "/opal-Tests"))[0]
    builddir = rundir + "/build"
    d = datetime.date.today()
    global totalNrPassed
    global totalNrTests
    totalNrTests = 0
    totalNrPassed = 0
    runAsUser = False
    runtests = list()

    if "--user" in argv:
        runAsUser = True
        #build a list of tests user wants to run
        for arg in argv:
            if arg.startswith("--tests"):
                tests = str.split(arg, "=")[1]
                runtests = str.split(tests, ",")
    else:
        #"load" modules need to compile and run regression tests
        modules = readfile("modules")
        for module in modules:
            module_load(module)
            os.environ["SGE_CELL"]="sgefelsim"
            os.environ["SGE_EXECD_PORT"]="6445"
            os.environ["SGE_QMASTER_PORT"]="6444"
            os.environ["SGE_ROOT"]="/gpfs/homefelsim/export/sge"
            os.environ["SGE_CLUSTER_NAME"]="sgeclusterfelsim"
            os.environ["PATH"]= os.getenv("PATH") + ":/gpfs/homefelsim/export/sge/bin/lx24-amd64:/usr/kerberos/bin"
            print commands.getoutput("/bin/env")
        #rep.appendReport(commands.getoutput("/bin/env"))
    
        #klog to be able to do svn stuff
        os.environ["KRB5_CONFIG"] = "/home2/l_felsimsvn/krb5.conf"
        commands.getoutput("/usr/kerberos/bin/kinit -V -k -t ~/.krb5.keytab.D.PSI.CH l_felsimsvn@D.PSI.CH")
        #commands.getoutput("kinit -k -t ~/.krb5.keytab l_felsimsvn")

    #check if user also wants build test
    if "--build" in argv:
        buildopal = Builder(srcdir, "OPAL", "src/opal", builddir)
        totalNrTests += 1
        if buildopal.build():
            totalNrPassed += 1
        else:
            rep.appendReport("Build-test failed! Exiting..")
            bailout(runAsUser)
            return

    #check if user has already set an OPAL executable
    #if not use the one from the last build test
    env = os.getenv("OPAL_EXE_PATH")
    if env is None:
        os.environ["OPAL_EXE_PATH"] = builddir + "/src"

    rep.appendReport("\n")
    rep.appendReport("Start Regression Test on %s \n" % datetime.datetime.today())
    rep.appendReport("==========================================================\n")

    #only run regression tests if opal executable is valid
    if os.path.isfile(os.getenv("OPAL_EXE_PATH") + "/opal"):
        os.chdir(regdir)
        #first update all tests and regression test files
        #FIXME: detect SCM
        commands.getoutput("svn update")
        #walk the run dir tree
        arglist = [runtests]
        os.path.walk("./", callback, arglist)
    else:
        rep.appendReport("Error: OPAL_EXE_PATH is invalid")
        bailout(runAsUser)
        return

    rep.dumpXML("results.xml")
    
    #cp report to webdir and add entry in index.html
    if not runAsUser:
        commands.getoutput("/usr/kerberos/bin/kinit -V -k -t ~/.krb5.keytab.D.PSI.CH l_felsimsvn@D.PSI.CH")
        commands.getoutput("/usr/bin/aklog")
        failedtests = rep.NrFailed()
        brokentests = rep.NrBroken()
        webfilename = "results_%s_%s_%s.xml" % (d.day, d.month, d.year)
        commands.getoutput("cp results.xml " + www_folder + "/" + webfilename)
        commands.getoutput("mv plots_" + d.isoformat() + " " + www_folder + "/")
        indexhtml = open(www_folder + "/index.html").readlines()
        for line in range(len(indexhtml)):
            if "insert here" in indexhtml[line]:
                indexhtml.insert(line+1, "<a href=\"%s\">%s.%s.%s</a> [passed:%d | broken:%d | failed:%d | total:%d] <br/>\n" % (webfilename, d.day, d.month, d.year, totalNrPassed, brokentests, failedtests, totalNrTests))
                break
        indexhtmlout = open(www_folder + "/index.html", "w")
        indexhtmlout.writelines(indexhtml)
        indexhtmlout.close()
        #update xslt formating file
        commands.getoutput("cp " + rundir + "/results.xslt " + www_folder + "/")

        #update manual
        OpalDocumentation()

        #update doxygen
        OpalDoxygen()

    #move xml results to result-dir
    if os.path.isfile('results.xml'):
        resultdir = regdir + "/results/" + d.isoformat()
        if not os.path.isdir(resultdir):
            commands.getoutput("mkdir -p " + resultdir)
    
        commands.getoutput("mv " + "results.xml " + resultdir)
    
    if runAsUser:
        commands.getoutput("rm -rf " + regdir + "/plots_" + d.isoformat())
    else:
        commands.getoutput("kdestroy")

    bailout(runAsUser)

#call main
if __name__ == "__main__":
    main(sys.argv[1:])
