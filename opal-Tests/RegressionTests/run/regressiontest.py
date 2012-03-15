import commands
import datetime
import os
import time

from reporter import Reporter
from reporter import TempXMLElement

from tools import genplot
from tools import readfile

class RegressionTest:
    
    def __init__(self, dir, simname, resultdir):
        self.dirname = dir
        self.simname = simname
        self.resultdir = resultdir
        self.jobnr = -1
        self.totalNrTests = 0
        self.totalNrPassed = 0

    def submitToSGE(self):
        #FIXME: we could create a sge file on the fly if no sge is specified for a give test ("default sge")
        #FIXME queue
        qsub_command = "qsub " + self.simname + ".sge -v REG_TEST_DIR=" + self.dirname + ",OPAL_EXE_PATH=" + os.getenv("OPAL_EXE_PATH")
        submit_out = commands.getoutput(qsub_command)
        self.jobnr = str.split(submit_out, " ")[2]

    def waitUntilCompletion(self):
        username = commands.getoutput("whoami")
        qstatout = commands.getoutput("qstat -u " + username + " | grep \"" + self.jobnr + "\"")
        while len(qstatout) > 0:
            #we only check every 30 seconds if job has finished
            time.sleep(30)
            qstatout = commands.getoutput("qstat -u " + username + " | grep \"" + self.jobnr + "\"")

    def performTests(self, root):
        rep = Reporter()
        tests = readfile(self.simname + ".rt")
        root.addAttribute("description", tests[0].lstrip("\"").rstrip("\""))
        rep.appendChild(root)

        tests = tests[1::] #strip first line
        for test in tests:
            self.totalNrTests += 1
            test_root = TempXMLElement("Test")
            passed = self.checkTest(test, test_root)
            if passed:
                self.totalNrPassed += 1
            root.appendChild(test_root)


    """
    handler for comparison of various output files with reference files
    """
    def checkTest(self, test, root):
        nameparams = str.split(test,"\"")
        var = nameparams[1]
        params = str.split(nameparams[2].lstrip(), " ")
        rtest = 0
        if "stat" in test:
            rtest = StatTest(var, params[0], float(params[1]), self.simname)
        elif "out" in test:
            rtest = OutTest(var, params[0], float(params[1]), self.simname)
        elif "lbal" in test:
            rtest = LbalTest(var, params[0], float(params[1]), self.simname)
        else:
            rep = Reporter()
            rep.appendReport("Error: unknown test type %s " % testparams[0])
            return False

        return rtest.performTest(root)

    def run(self, root):
        curd = os.getcwd()
        os.chdir(self.dirname)

        isValid = self.validateReferenceFiles()
 
        if isValid:
            rep = Reporter()
            rep.appendReport("\t run simulation\n")
            os.environ["REG_TEST_DIR"] = self.dirname

            # cleanup all OLD job files if there are any
            commands.getoutput("rm " + self.simname + "-RT.*")

            self.submitToSGE()
            self.waitUntilCompletion()

            # copy o to out file
            commands.getoutput("cp " + self.simname + "-RT.o* " + self.simname + ".out")

            self.performTests(root)
            
            # move plots to plot dir
            d = datetime.date.today()
            plotdir = "plots_" + d.isoformat()
            commands.getoutput("mkdir " + curd + "/" + plotdir)
            commands.getoutput("mv *.png " + curd + "/" + plotdir)

            #move tests to result folder
            commands.getoutput("mv " + self.simname + ".stat " + curd + "/" + self.resultdir) 
            commands.getoutput("mv " + self.simname + ".lbal " + curd + "/" + self.resultdir) 
            commands.getoutput("mv " + self.simname + ".out " + curd + "/" + self.resultdir) 
            commands.getoutput("cp -R ../" + plotdir + " "  + curd + "/results/" + d.isoformat()) 

        os.chdir(curd)

    """
    This method checks if all files in the reference directory are present
    and if their md5 checksums still concure with the ones stored after
    the simulation run
    """
    def validateReferenceFiles(self):
        rep = Reporter()
        olddir = os.getcwd()
        os.chdir("reference")
        allok = False

        if os.path.isfile(self.simname + ".stat") and \
           os.path.isfile(self.simname + ".stat.md5") and \
           os.path.isfile(self.simname + ".out") and \
           os.path.isfile(self.simname + ".out.md5") and \
           os.path.isfile(self.simname + ".lbal") and \
           os.path.isfile(self.simname + ".lbal.md5"):

            statout = commands.getoutput("md5sum --check " + self.simname + ".stat.md5")
            outout = commands.getoutput("md5sum --check " + self.simname + ".out.md5")
            lbalout = commands.getoutput("md5sum --check " + self.simname + ".lbal.md5")

            rep.appendReport("\t Checksum for reference %s \n" % statout)
            rep.appendReport("\t Checksum for reference %s \n" % outout)
            rep.appendReport("\t Checksum for reference %s \n" % lbalout)
            allok = statout == self.simname + ".stat: OK" and outout == self.simname + ".out: OK" and lbalout == self.simname + ".lbal: OK"

        else:
            rep.appendReport("\t Error: reference dir is incomplete! \n")

        os.chdir(olddir)
        return allok


class StatTest:

    def __init__(self, var, quant, eps, simname):
        self.var = var
        self.quant = quant
        self.eps = eps
        self.simname = simname

    """
    method parses a stat-file and returns found variable
    """
    def readStatVariable(self, filename):
        vars = []
        nrCol = 0
        numScalars = 0
        hasReadHeader = False
        lines = readfile(filename + ".stat")

        for line in lines:
            name = "name=" + self.var
            if line.find(name) != -1: #find offset in stat list
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
    method performs a test for a stat-file variable "var"
    """
    def performTest(self, root):

        rep = Reporter()
        val = 0

        root.addAttribute("type", "stat")
        root.addAttribute("var", self.var)
        root.addAttribute("mode", self.quant)
        passed_report = TempXMLElement("passed")
        eps_report = TempXMLElement("eps")
        delta_report = TempXMLElement("delta")
        plot_report = TempXMLElement("plot")

        if not os.path.isfile(self.simname + ".stat"):
            rep.appendReport("ERROR: no statfile %s \n" % self.simname)
            rep.appendReport("\t Test %s(%s) broken \n" % (self.var,self.quant))
            passed_report.appendTextNode("false")
            delta_report.appendTextNode("-")
            eps_report.appendTextNode("%s" % self.eps)

            root.appendChild(passed_report)
            root.appendChild(eps_report)
            root.appendChild(delta_report)
            return False

        readvar_sim = self.readStatVariable(self.simname)
        readvar_ref = self.readStatVariable("reference/" + self.simname)

        plotfilename = genplot(self.simname, self.var)
    
        if readvar_sim == [] or readvar_ref == []:
            rep.appendReport("Error: unknown variable (%s) selected for stat test\n" % self.var)
            rep.appendReport("\t Test %s(%s) broken: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("false")
            delta_report.appendTextNode("-")
            eps_report.appendTextNode("%s" % self.eps)

            root.appendChild(passed_report)
            root.appendChild(eps_report)
            root.appendChild(delta_report)
            return False
        
        if len(readvar_sim) != len(readvar_ref):
            rep.appendReport("Error: size of stat variables (%s) dont agree!\n" % self.var)
            rep.appendReport("\t Test %s(%s) broken: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("false")
            delta_report.appendTextNode("-")
            eps_report.appendTextNode("%s" % self.eps)

            root.appendChild(passed_report)
            root.appendChild(eps_report)
            root.appendChild(delta_report)
            return False

        if self.quant == "last":
            val = abs(readvar_sim[len(readvar_sim) -1] - readvar_ref[len(readvar_sim) -1])

        elif self.quant == "avg":
            sum = 0.0
            for i in range(len(readvar_sim)):
                sum += (readvar_sim[i] - readvar_ref[i])**2
            val = (sum)**(0.5) / len(readvar_sim)

        elif self.quant == "error":
            rep.appendReport("TODO: error norm\n")

        elif self.quant == "all":
            rep.appendReport("TODO: graph/all\n")

        else:
            rep.appendReport("Error: unknown quantity %s \n" % self.quant)
            
        #result generation
        passed = False
        if val < self.eps:
            rep.appendReport("Test %s(%s) passed: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("true")
            passed = True
        else:
            rep.appendReport("Test %s(%s) failed: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("false")

        delta_report.appendTextNode("%s" % val)
        eps_report.appendTextNode("%s" % self.eps)

        root.appendChild(passed_report)
        root.appendChild(eps_report)
        root.appendChild(delta_report)

        if plotfilename != "":
            plot_report.appendTextNode(plotfilename)
            root.appendChild(plot_report)

        return passed


class OutTest:
    
    def __init__(self, var, quant, eps, simname):
        self.var = var
        self.quant = quant
        self.eps = eps
        self.simname = simname

    """
    method parses an out-file and returns found variables as tuples
    """
    def readOutVariable(self, filename):
        vars = []
        nrCol = 0
        numScalars = 0
        lines = readfile(filename + ".out")

        for line in lines:
            if self.var in line: 
                varline = str.split(line, "=")
                value = ""
                for i in range (len(varline)):
                    if self.var in varline[i]:
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
    def performTest(self, root):
        rep = Reporter()
        val = list()
        passed = True
    
        #report stuff
        root.addAttribute("type", "out")
        root.addAttribute("var", self.var)
        root.addAttribute("mode", self.quant)
        passed_report = TempXMLElement("passed")
        eps_report = TempXMLElement("eps")
        delta_report = TempXMLElement("delta")

        if not os.path.isfile(self.simname + ".out"):
            rep.appendReport("ERROR: no outfile %s \n" % self.simname)
            rep.appendReport("\t Test %s(%s) broken\n" % (self.var,self.quant))
            passed_report.appendTextNode("false")
            delta_report.appendTextNode("-")
            eps_report.appendTextNode("%s" % self.eps)

            root.appendChild(passed_report)
            root.appendChild(eps_report)
            root.appendChild(delta_report)
            return False

        #get ref and sim variable values
        readvar_sim = self.readOutVariable(self.simname)
        readvar_ref = self.readOutVariable("reference/" + self.simname)

        if len(readvar_sim) == 0 or len(readvar_ref) == 0:
            rep.appendReport("Error: unknown variable (%s) selected for out test\n" % self.var)
            rep.appendReport("\t Test %s(%s) broken: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("false")
            delta_report.appendTextNode("-")
            eps_report.appendTextNode("%s" % self.eps)

            root.appendChild(passed_report)
            root.appendChild(eps_report)
            root.appendChild(delta_report)
            return False

        if len(readvar_sim) != len(readvar_ref):
            rep.appendReport("Error: size of out variables (%s) dont agree!\n" % self.var)
            rep.appendReport("\t Test %s(%s) broken: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("false")
            delta_report.appendTextNode("-")
            eps_report.appendTextNode("%s" % self.eps)

            root.appendChild(passed_report)
            root.appendChild(eps_report)
            root.appendChild(delta_report)
            return False

        if self.quant == "last":
            for i in range (len(readvar_sim[0])):
                val.append( abs(readvar_sim[len(readvar_sim) -1][i] - readvar_ref[len(readvar_sim) -1][i]) )

            for i in range (len(readvar_sim[0])):
                passed = passed and (val[i] < self.eps)

        elif self.quant == "avg":
            if len(readvar_sim) != len(readvar_ref):
                rep.appendReport("Error: size of stat variables dont agree!\n")
                return

            for j in range (len(readvar_sim[0])): #number of components 
                sum = 0.0
                for i in range(len(readvar_sim)): #number of entries
                    sum += (readvar_sim[i][j] - readvar_ref[i][j])**2
            
                val.append((sum)**(0.5) / len(readvar_sim))

            for i in range (len(readvar_sim[0])):
                passed = passed and (val[i] < self.eps)

        elif self.quant == "error":
            rep.appendReport("TODO: error norm\n")

        elif self.quant == "all":
            rep.appendReport("TODO: graph/all\n")

        else:
            rep.appendReport("Error: unknown quantity %s \n" % self.quant)

        #result generation
        if passed:
            rep.appendReport("Test %s(%s) passed: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("true")
        else:
            rep.appendReport("Test %s(%s) failed: %s (eps=%s) \n" % (self.var,self.quant,val,self.eps))
            passed_report.appendTextNode("false")
    
        if len(val) == 1:
            delta_report.appendTextNode("%s" % val[0])
        else:
            delta_report.appendTextNode("%s" % val)
        eps_report.appendTextNode("%s" % self.eps)

        root.appendChild(passed_report)
        root.appendChild(eps_report)
        root.appendChild(delta_report)

        return passed


class LbalTest:
    
    def __init__(self, var, quant, eps, simname):
        self.var = var
        self.quant = quant
        self.eps = eps
        self.simname = simname


    def readLbalFile(self, filename):

        lbal = []
        lines = readfile(filename + ".lbal")
        #nrprocs = int(lines[0])
        lines = lines[1::]
        for line in lines:
            vals = str.split(line, "\t")
            lbal.append(vals)

    """
    method performs a test for "var" with reference file in a specific mode ("quant") for a specific accuracy ("eps")
    """
    def performTest(self, root):

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

        return False


