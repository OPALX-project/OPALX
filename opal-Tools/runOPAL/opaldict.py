"""
OpalDictionary class

@author: Andreas Adelmann <andreas.adelmann@psi.ch>
@author: Yves Ineichen <yves.ineichen@psi.ch>
@version: 0.1
"""

class OpalDict:

    def __init__(self, template):
        self.dict = {}
        self.rangevars = {}
        self.uservars = []
        self.numRanges = 0
        self.fillDictionary(template)

    def __iter__(self):
        return self.dict.__iter__()

    def __setitem__(self, key, value):
        scalevars = {}
        scalevars['GUNSOLB'] =  1.E-3
        scalevars['TRISE'] = 1.E-12
        scalevars['TFALL'] = 1.E-12
        scalevars['TFWHM'] = 1.E-12
        scalevars['QBUNCH'] = 1.E-12
        scalevars['SIGX'] = 1.E-6
        scalevars['SIGY'] = 1.E-6
        scalevars['DT'] = 1.E-12
        scalevars['DTGUN'] = 1.E-12
        scalevars['DTAP'] = 1.E-12

        try:
            self.dict[key] = value * scalevars[key]
        except KeyError:
            self.dict[key] = value

    def __getitem__(self, key):
        return self.dict[key]

    def fillDictionary(self, fileName):
        fp = open(fileName,"r")
        for line in fp:
            li = line.strip()
            if not li.startswith("#"):
                aline = line.split("#")[0]
                name,val = aline.split()
                self.dict[name.rstrip()] = val.lstrip().rstrip()	
        fp.close()

    def generateDirectoryName(self):
        dirname = ""
        for p in self.uservars:
            dirname += "_" + str(p[0]) + "=" + str(p[1])
        for (k, v) in self.rangevars.items():
            dirname += "_" + str(k) + "=" + str(self.dict[k])

        return dirname

    def scaleDictVar(self, var, scaleWith):
        if self.dict.has_key(var):
            self.dict[var] = float(self.dict[var])*scaleWith

    def getType(self, value):
        if value.find('.') > 0:
            return float(value)
        else:
            return int(value)
    
    def hasRanges(self):
        return self.numRanges > 0

    def Range(self):
        return self.rangevars

    def scale(self):
        self.scaleDictVar('GUNSOLB', 1.E-3)
        self.scaleDictVar('TRISE', 1.E-12)
        self.scaleDictVar('TFALL', 1.E-12)
        self.scaleDictVar('TFWHM', 1.E-12)
        self.scaleDictVar('QBUNCH', 1.E-12)

        self.scaleDictVar('SIGX', 1.E-6)
        self.scaleDictVar('SIGY', 1.E-6)

        self.scaleDictVar('DT', 1.E-12)
        self.scaleDictVar('DTGUN', 1.E-12)
        self.scaleDictVar('DTAP', 1.E-12)
    
    def addUserValues(self, argv):
        for arg in argv:
            if arg.find("=") > 0:
                var = arg.split("=")[0]
                num = arg.split("=")[1]
                if self.dict.has_key(var):
                    #check if we have a range
                    if num.find(':') > 0:
                        range = num.split(":")
                        if len(range) == 3:
                            rvar = []
                            for r in range:
                                rvar.append(self.getType(r))
                            self.rangevars[var] = rvar
                            self.numRanges = self.numRanges + 1
                        else:
                            print "OpalDict: Range has to be of the form from:to:step!"
                    else:
                        self.uservars.append( (var, num) )
                        self.dict[var] = self.getType(num)
                else: 
                    print 'OpalDict: Key (' + var + ')not found can not add to dictionary, check the OPAL template file'
