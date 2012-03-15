#!/usr/bin/python 
import commands
import os
import datetime

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
    passedtests = protocol.count("passed")

    MAIL = "/usr/sbin/sendmail"
    p = os.popen("%s -t" % MAIL, "w")
    for adr in addresses:
        p.write("To: %s \n" % adr)
    p.write("Subject: [Regression Tests] passed %d/%d build tests on %s" % (passedtests,nrtest,d.isoformat()))
    p.write("\n")
    p.write(protocol)
    sts = p.close()
    if sts != 0:
        print "Mail Sent"

"""
later we fit this into class that has svn info as members
and pass to this function the desired revision number we 
want to build & test. With this it will be easy to backtrack
to the last working version
"""
def buildtest(dir, name, target):
    curdir = os.getcwd()
    os.chdir(dir)

    log = "";
    built = True;
    
    outtext = commands.getoutput('svn -uq status')
    info = commands.getoutput('svn info')

    #delete target to force relinking
    commands.getoutput("rm " + target)

    if outtext.startswith("Status") == 0:
        log += "%s has changed \n" % name
        outtext = commands.getoutput("svn update")
        info = commands.getoutput('svn info')
        log += "\t make -k maintainer-clean \n"
        outtext = commands.getoutput("make -k maintainer-clean")
        log += "\t Running autogen and configure \n"
        outtext = commands.getoutput("./autogen-regression.sh")
        if os.path.isfile("Makefile") == 0:
            log += "\t configure failed: \n"
            log += outtext
            built = False
    else:
        log += "%s up to date: no changes \n" % name
    
    #always make to force relinking
    log += "\t make \n"
    outtext = commands.getoutput("make -j 4")
    if os.path.isfile(target) == 0:
        log += "\t make failed: \n"
        log += outtext
        built = False

    infolines = str.split(info, "\n")
    revision = str.split(infolines[4], ":")[1].lstrip()
    if built == True:
        log += "%s (r%s) passed build test \n" % (name, revision)
    else:
        log += "%s (r%s) failed build test \n" % (name, revision)

    os.chdir(curdir)
    return (log, revision)

"""
####################### REG TEST ##########################
"""

class ResultType:
    unknown = 0
    last = 1
    average = 2
    all = 3
    error = 4

class TestType:
    unknown = 0
    stat = 5
    lbal = 6
    out = 7

"""
handler for comparison of various output files with reference files
"""
def checktest(test, simname, test_report, xml_report):
    nameparams = str.split(test,"\"")
    var = nameparams[1]
    params = str.split(nameparams[2].lstrip(), " ")
    if "stat" in test:
        return stattest(var, params[0], float(params[1]), simname, test_report, xml_report)
    elif "out" in test:
        return outtest(var, params[0], float(params[1]), simname, test_report, xml_report)
    elif "lbal" in test:
        return lbaltest(var, params[0], float(params[1]), simname, test_report, xml_report)
    else:
        return "Error: unknown test type %s " % testparams[0]

"""
method performs a test for a stat-file variable "var"
"""
def stattest(var, quant, eps, simname, test_report, xml_report):

    readvar_sim = readStatVariable(var, quant, simname)
    readvar_ref = readStatVariable(var, quant, "reference/" + simname)
    report = ""
    val = 0

    test_report.setAttribute("type", "stat")
    test_report.setAttribute("var", var)
    test_report.setAttribute("mode", quant)
    passed_report = xml_report.createElement("passed")
    eps_report = xml_report.createElement("eps")
    delta_report = xml_report.createElement("delta")
    
    if readvar_sim == [] or readvar_ref == []:
        report += "Error: unknown variable (%s) selected for stat test\n" % var
        report += "\t Test %s(%s) failed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("false")
        delta_d = xml_report.createTextNode("-")
        eps_d = xml_report.createTextNode("%s" % eps)

        passed_report.appendChild(passed_d)
        eps_report.appendChild(eps_d)
        delta_report.appendChild(delta_d)

        test_report.appendChild(passed_report)
        test_report.appendChild(eps_report)
        test_report.appendChild(delta_report)
        return report
        
    if len(readvar_sim) != len(readvar_ref):
        report += "Error: size of stat variables (%s) dont agree!\n" % var
        report += "\t Test %s(%s) failed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("false")
        delta_d = xml_report.createTextNode("-")
        eps_d = xml_report.createTextNode("%s" % eps)

        passed_report.appendChild(passed_d)
        eps_report.appendChild(eps_d)
        delta_report.appendChild(delta_d)

        test_report.appendChild(passed_report)
        test_report.appendChild(eps_report)
        test_report.appendChild(delta_report)
        return report

    if quant == "last":
        val = (readvar_sim[len(readvar_sim) -1] - readvar_ref[len(readvar_sim) -1])**2

    elif quant == "avg":
        sum = 0.0
        for i in range(len(readvar_sim)):
            sum += (readvar_sim[i] - readvar_ref[i])**2
        val = sum / len(readvar_sim)

    elif quant == "error":
        report += "TODO: error norm\n"

    elif quant == "all":
        report += "TODO: graph/all\n"

    else:
        report += "Error: unknown quantity %s \n" % quant
            
    #result generation
    if val < eps:
        report += "Test %s(%s) passed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("true")
    else:
        report += "Test %s(%s) failed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("false")

    delta_d = xml_report.createTextNode("%s" % val)
    eps_d = xml_report.createTextNode("%s" % eps)

    passed_report.appendChild(passed_d)
    eps_report.appendChild(eps_d) 
    delta_report.appendChild(delta_d)

    test_report.appendChild(passed_report)
    test_report.appendChild(eps_report)
    test_report.appendChild(delta_report)

    return report

"""
method parses a stat-file and returns found variable
"""
def readStatVariable(var, quant, simname):

    vars = []
    nrCol = 0
    numScalars = 0
    hasReadHeader = False
    lines = readfile(simname + ".stat")

    for line in lines:

        if var in line: #find offset in stat list
            param = str.split(line, "description=\"")[1]
            nrCol = int(str.split(param, " ")[0])-1

        elif "&parameter" in line:
            numScalars += 1

        elif "&data mode=ascii &end" in line:
            hasReadHeader = True

        #FIXME: this is very ugly
        elif hasReadHeader == True and numScalars > 0:
            numScalars -= 1
            continue

        elif hasReadHeader == True and numScalars == 0:
            values = str.split(line, "\t")
            vars.append(float(values[nrCol]))

    return vars

"""
method performs a test for "var" with reference file in a specific mode ("quant") for a specific accuracy ("eps")
"""
def outtest(var, quant, eps, simname, test_report, xml_report):

    #get ref and sim variable values
    readvar_sim = readOutVariable(var, quant, simname)
    readvar_ref = readOutVariable(var, quant, "reference/" + simname)
    report = ""
    val = list()
    passed = True
    
    #report stuff
    test_report.setAttribute("type", "out")
    test_report.setAttribute("var", var)
    test_report.setAttribute("mode", quant)
    passed_report = xml_report.createElement("passed")
    eps_report = xml_report.createElement("eps")
    delta_report = xml_report.createElement("delta")

    if len(readvar_sim) == 0 or len(readvar_ref) == 0:
        report += "Error: unknown variable (%s) selected for out test\n" % var
        report += "\t Test %s(%s) failed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("false")
        delta_d = xml_report.createTextNode("-")
        eps_d = xml_report.createTextNode("%s" % eps)

        passed_report.appendChild(passed_d)
        eps_report.appendChild(eps_d)
        delta_report.appendChild(delta_d)

        test_report.appendChild(passed_report)
        test_report.appendChild(eps_report)
        test_report.appendChild(delta_report)
        return report

    if len(readvar_sim) != len(readvar_ref):
        report += "Error: size of out variables (%s) dont agree!\n" % var
        report += "\t Test %s(%s) failed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("false")
        delta_d = xml_report.createTextNode("-")
        eps_d = xml_report.createTextNode("%s" % eps)

        passed_report.appendChild(passed_d)
        eps_report.appendChild(eps_d)
        delta_report.appendChild(delta_d)

        test_report.appendChild(passed_report)
        test_report.appendChild(eps_report)
        test_report.appendChild(delta_report)
        return report

    if quant == "last":
        for i in range (len(readvar_sim[0])):
            val.append( (readvar_sim[len(readvar_sim) -1][i] - readvar_ref[len(readvar_sim) -1][i])**2 )

        for i in range (len(readvar_sim[0])):
            passed = passed and (val[i] < eps)

    elif quant == "avg":
        if len(readvar_sim) != len(readvar_ref):
            report += "Error: size of stat variables dont agree!\n"
            return report

        for j in range (len(readvar_sim[0])): #number of components 
            sum = 0.0
            for i in range(len(readvar_sim)): #number of entries
                sum += (readvar_sim[i][j] - readvar_ref[i][j])**2
            
            val.append(sum / len(readvar_sim))

        for i in range (len(readvar_sim[0])):
            passed = passed and (val[i] < eps)

    elif quant == "error":
        report += "TODO: error norm\n"

    elif quant == "all":
        report += "TODO: graph/all\n"

    else:
        report += "Error: unknown quantity %s \n" % quant

    #result generation
    if passed:
        report += "Test %s(%s) passed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("true")
    else:
        report += "Test %s(%s) failed: %s (eps=%s) \n" % (var,quant,val,eps)
        passed_d = xml_report.createTextNode("false")
    
    if len(val) == 1:
        delta_d = xml_report.createTextNode("%s" % val[0])
    else:
        delta_d = xml_report.createTextNode("%s" % val)
    eps_d = xml_report.createTextNode("%s" % eps)

    passed_report.appendChild(passed_d)
    eps_report.appendChild(eps_d)
    delta_report.appendChild(delta_d)

    test_report.appendChild(passed_report)
    test_report.appendChild(eps_report)
    test_report.appendChild(delta_report)

    return report

"""
method parses an out-file and returns found variables as tuples
"""
def readOutVariable(var, quant, simname):
    vars = []
    nrCol = 0
    numScalars = 0
    lines = readfile(simname + ".out")

    for line in lines:

        if var in line: 

            varline = str.split(line, "=")
            value = ""
            for i in range (len(varline)):
                if var in varline[i]:
                    value = varline[i+1].lstrip()
                    if value.startswith("("): #vector
                        vec = str.split(value, " ")
                        retval = (float(vec[1]), float(vec[3]), float(vec[5]))
                    else: #scalar
                        retval = (float(str.split(value, " ")[0]),)

                    vars.append(retval)
                    break;
    return vars

"""
method performs a test for "var" with reference file in a specific mode ("quant") for a specific accuracy ("eps")
"""
def lbaltest(var, quant, eps, simname, test_report, xml_report):

    readvar_sim = readLbalFile(simname)
    readvar_ref = readLbalFile("reference/" + simname)

    if len(readvar_sim) != len(readvar_ref):
        print "Error: size does not agree!"

    if quant == "all":
        for i in range(len(readvar_sim)):
            #TODO: what delta?
            print "calc some delta"
        
    else:
        return "Error: please only use all for lbal files!"

def readLbalFile(simname):

    lbal = []
    lines = readfile(simname + ".lbal")
    #nrprocs = int(lines[0])
    lines = lines[1::]
    for line in lines:
        vals = str.split(line, "\t")
        lbal.append(vals)


