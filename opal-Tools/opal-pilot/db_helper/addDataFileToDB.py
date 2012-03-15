"""
This script adds all parameters found in a data file to the DB
Before you execute that script, better get a database backup.
"""

# export PYTHONPATH=~/svnwork/OPAL/opal-Tools/opal-pilot/src/ for this to work!

import os
import sys
import _mysql
import _mysql_exceptions

from util.database import DatabaseConnection

def addToDB(parameter, ui_label, ui_tooltip, ui_group_id, ui_unit_id, ui_default_value, acc_id):

    parameters = [(parameter, ui_label, ui_tooltip, ui_group_id, ui_unit_id), ]

    conn = DatabaseConnection()
    
    if 1:
        for tuple in parameters:
            try:
                pid = conn.insert('parameter', ['parameter', 'ui_label', 'ui_tooltip', 'ui_group_id', 'ui_unit_id'], [tuple])
                conn.insert('accelerator_parameter', ['accelerator_id', 'parameter_id', 'ui_default_value'], [(acc_id, int(pid), ui_default_value)])
            except _mysql_exceptions.IntegrityError:
                query = """
                UPDATE parameter SET ui_label = '%s', ui_tooltip = '%s', ui_group_id = %s, ui_unit_id = %s WHERE parameter = '%s'
                """ % (tuple[1], tuple[2], tuple[3], tuple[4], tuple[0])
                print query            
                conn.execute(query)

def createDictionary(filename):
    di = {}
    fp = open(filename,"r")
    for line in fp:
        li = line.strip()
        if not li.startswith("#"):
            name,val = li.partition("#")[0].split()
            di[name.strip()] = val.strip()	
    fp.close()
    return di


def createElementDictionary(datadict, tmplfilename):
    di = {}
    fp = open(tmplfilename,"r")
    for line in fp:
        li = line.strip()
        for element in datadict:
            elename = element.rsplit("_",1)[0]
            if li.startswith(elename + ":"):
                di[elename] =  (li.split(":")[1]).strip().split(",")[0]
    fp.close()
    return di


"""
main method
"""
def main(argv):
    
    inputfile = ""
    #if len(argv) < 1:
    #    print "Usage: ./addInputfileToDB.py inputfile.data"
    #else:
    #    inputfile = argv[0]

    datafilename = "/home2/amas/Simulations/SwissFELInjectorTestFacility/Phase3/FinPhase3.data"
    tmplfilename = "/home2/amas/Simulations/SwissFELInjectorTestFacility/Phase3/FinPhase3.tmpl"

    # find out what accelerator we are working on
    accelerator = (datafilename.split("/")[1]).split(".")[0]
    query = """
    SELECT acclerator_id FROM accelerator WHERE accelerator_name = '%s'""" % (accelerator)
    print query            
    #accelerator_id = conn.execute(query)
    accelerator_id = 4

    # ui group and units (TODO: later get from DB)
    ui_group = {'RFCAVITY' : 10, 'TRAVELINGWAVE' : 11, 'SOLENOID' : 9, 'QUADRUPOLE' : 8}
    ui_unit = {'lag' : 1, 'volt' : 2, 'k1' : 11, 'ks' : 10, 'dx' : 4, 'dy' : 4, 'dz' : 4, 'ds' : 4}

    # maps element to default value (all parameters we need to add are in this dictionary)
    dict = createDictionary(datafilename)
    # maps element name to element type
    elementdict = createElementDictionary(dict, tmplfilename)

    # finally add every element in dict to the DB
    for entry in dict:
        elename = entry.rsplit("_",1)[0]
        parameter = entry
        ui_label = entry
        ui_tooltip = ""
        
        try:
            ui_group_id = ui_group[elementdict[elename]]
        except KeyError:
            ui_group_id = ""

        try:
            ui_unit_id = ui_unit[entry.split('_')[-1]]
        except KeyError:
            ui_unit_id = ""

        try:
            ui_default_value = dict[entry]
        except KeyError:
            ui_default_value = 0

        addToDB(parameter, ui_label, ui_tooltip, ui_group_id, ui_unit_id, ui_default_value, accelerator_id)


#call main
if __name__ == "__main__":
    main(sys.argv[1:])
