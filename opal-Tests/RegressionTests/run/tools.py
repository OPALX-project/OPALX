#!/usr/bin/python 
import commands
import os
import datetime
import subprocess

"""
parse linefile into list and strip newline 
char at the end
"""
def readfile(fname):
    list = []
    infile = open(fname,"r")
    #TODO: TRY
    #list = open(fname, "r").readlines()
    while infile:
        line = infile.readline()
        if not line: break
        list.append(line.rstrip('\n'))
    infile.close()
    return list

"""
sends the "protocol" of "nrtest" tests to all addresses in "addresses"
"""
def sendmails(addresses, protocol, nrtest):
    d = datetime.date.today()
    failedtests = protocol.count("failed")
    brokentests = protocol.count("broken")

    MAIL = "/usr/sbin/sendmail"
    p = os.popen("%s -t" % MAIL, "w")
    for adr in addresses:
        p.write("To: %s \n" % adr)
    p.write("Subject: [Regression Tests] BROKEN %d/%d | FAILED %d/%d tests on %s" % (brokentests,nrtest,failedtests,nrtest,d.isoformat()))
    p.write("\n")
    p.write(protocol)
    sts = p.close()
    if sts != 0:
        print "Mail Sent"


"""
simulates a "module load" by reading the new and to append environment variables
from the "module show" output
"""
def module_load(module):
    moduleout = commands.getoutput("source /opt/Modules/init/bash && module show " + module)
    lines = str.split(moduleout, "\n")
    for line in lines:
        if line.startswith("setenv"):
            sete = str.split(line, "\t")[2].lstrip().rstrip()
            env = str.split(sete, " ")
            #print "export " + env[0] + "=" + env[1]
            #os.putenv(env[0], env[1])
            os.environ[env[0]] = env[1]

        elif line.startswith("prepend-path"):
            sete = str.split(line, "\t")[1].lstrip().rstrip()
            envadd = str.split(sete, " ")
            #print "export " + envadd[0] + "=" + envadd[1] + ":" + envadd[0]
            if os.getenv(envadd[0]):
                #os.putenv(envadd[0], envadd[1] + ":" + os.getenv(envadd[0]))
                os.environ[envadd[0]] = envadd[1] + ":" + os.getenv(envadd[0])
            else:
                #os.putenv(envadd[0], envadd[1])
                os.environ[envadd[0]] = envadd[1]

"""
generate stat plot with gnuplot
returns fileame
"""
def genplot(simname, var):
    simnames = simname + ".stat"
    reference = "reference/" + simnames
    
    vars = []
    nrCol = -1
    varUnit = ''
    varParts = str.split(var, "_")
    prettyVar = varParts[0]
    if len(varParts) == 2:
        prettyVar = varParts[0] + "(" + varParts[1] + ")"
    numScalars = 0
    revLine = 0
    hasReadHeader = False
    readLines = 0
    opalRevision = ''
    refRevision = 'reference'
    lines = readfile(simname + ".stat")

    for line in lines:
        name = "name=" + var
        if line.find(name) != -1: #find offset in stat list
            param = str.split(line, "description=\"")[1]
            nrCol = int(str.split(param, " ")[0]) - 1
            unit = str.split(line, "units=")[1]
            varUnit = str.split(unit, " ")[0]

        elif "&parameter" in line:
            numScalars += 1
            if "name=revision" in line:
                revLine = numScalars

        elif "&data mode=ascii &end" in line:
            break

        readLines += 1

    opalRevision = lines[readLines + revLine]

    # a simple paste in gnuplot does not work since the
    # headers of the files do not necessarily have to have
    # the same length. therefore copy the data in temporary
    # files, delete them in the end
    if nrCol > -1:
        data1 = open('data1.dat','w')
        for line in lines[readLines + numScalars + 1:]:
            values = line.split()
            data1.write(values[1] + "\t" + values[nrCol] + "\n")
        data1.close()

    revLine = 0
    readLines = 0
    numScalars = 0
    
    lines = readfile(reference)

    for line in lines:
        if "&parameter" in line:
            numScalars += 1
            if "name=revision" in line:
                revLine = numScalars
        elif "&data mode=ascii &end" in line:
            break

        readLines += 1

    if revLine > 0:
        refRevision += ": " + lines[readLines + revLine]

    if nrCol > -1:
        data2 = open('data2.dat','w')
        for line in lines[readLines + numScalars + 1:]:
            values = line.split()
            data2.write(values[1] + "\t" + values[nrCol] + "\n")
        data2.close()
    
    filename = ""
    d = datetime.date.today()
    if nrCol > -1:
        filename = simname + "_" + var + "_" + str(d.day) + "_" + str(d.month) + "_" + str(d.year)
        plotcmd = "set terminal post enh col 20\n"
        plotcmd += "set output '" + filename + ".ps'\n"
        plotcmd += "set title '" + simname + "'\n"
        plotcmd += "set key below\n"
        plotcmd += "set ytics nomirror\n"
        plotcmd += "set y2tics\n"
        plotcmd += "set ylabel '" + prettyVar + " [" + varUnit + "]' font 'Helvetica-Bold,20'\n"
        plotcmd += "set y2label '{/Symbol D}" + prettyVar + " [" + varUnit + "]' font 'Helvetica-Bold,20'\n"
        plotcmd += "set xlabel 's [m]' font 'Helvetica-Bold,20'\n"
        plotcmd += "plot 'data1.dat' u 1:2 w l lw 2 t '" + opalRevision + "', "
        plotcmd += "'data2.dat' u 1:2 w l lw 2 t '" + refRevision + "', "
        plotcmd += "\"< paste data1.dat data2.dat\" u 1:($2-$4) w l lw 2 axis x1y2 t 'difference'" + ";\n"
        plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
        plot.communicate(plotcmd)

        commands.getoutput("convert -rotate 90 " + filename + ".ps " + filename + ".png")
        commands.getoutput("rm " + filename + ".ps data1.dat data2.dat")
    else:
        print "Error in genplot: Cannot find stat variable!" 

    return "plots_" + d.isoformat() + "/" + filename + ".png"
