"""
Before you execute that script, better get a database backup.                                                                                                                                                                                                                        
"""

# export PYTHONPATH=~/svnwork/OPAL/opal-Tools/opal-pilot/src/ for this to work!

import os
import sys
import _mysql
import _mysql_exceptions

from util.database import DatabaseConnection


def addToDB(parameter, ui_group, ui_unit): 

    parameters = [(parameter, parameter, "", ui_group, ui_unit), ]
    conn = DatabaseConnection()
    
    if 1:
        for tuple in parameters:
            try:
                pid = conn.insert('parameter', ['parameter', 'ui_label', 'ui_tooltip', 'ui_group_id', 'ui_unit_id'], [tuple])
                conn.insert('accelerator_parameter', ['accelerator_id', 'parameter_id', 'ui_default_value'], [(3, int(pid), 0.0)])
            except _mysql_exceptions.IntegrityError:
                query = """
                UPDATE parameter SET ui_label = '%s', ui_tooltip = '%s', ui_group_id = %s, ui_unit_id = %s WHERE parameter = '%s'
                """ % (tuple[1], tuple[2], tuple[3], tuple[4], tuple[0])
                print query            
                conn.execute(query)


"""
main method
"""
def main(argv):
    
    inputfile = ""
    #if len(argv) < 1:
    #    print "Usage: ./addInputfileToDB.py inputfile.in"
    #else:
    #    inputfile = argv[0]

    elements = ['RFCAVITY', 'TRAVELINGWAVE', 'SOLENOID', 'QUADRUPOLE']
    params = ['_ds', '_dx', '_dy', '_dz']
    baseNames = []
    filename = "/home2/amas/Simulations/SwissFELInjectorTestFacility/Phase3/FinPhase3.tmpl"
    
    ui_group = {'RFCAVITY' : 10, 'TRAVELINGWAVE' : 11, 'SOLENOID' : 9, 'QUADRUPOLE' : 8} #TODO: later read from DB!
    ui_unit = {'_lag' : 1, '_volt' : 2, '_k1' : 11, '_ks' : 10, '_dx' : 4, '_dy' : 4, '_dz' : 4, '_ds' : 4} #TODO: later read from DB!

    list = []
    infile = open(filename, "r")
    while infile:
        line = infile.readline()
        #TODO read until ';' reached!
        if not line: break
        for element in elements:
            if line.find(element) != -1 and not line.startswith("//"):
                element_name = line.split(":")[0]
                for param in params:
                    addToDB(element_name + param, ui_group[element], ui_unit[param])

                if element == 'RFCAVITY' or element == 'TRAVELINGWAVE':
                    addToDB(element_name + "_volt", ui_group[element], ui_unit["_volt"])
                    addToDB(element_name + "_lag", ui_group[element], ui_unit["_lag"])
                if element == 'SOLENOID':
                    addToDB(element_name + "_ks", ui_group[element], ui_unit["_ks"])
                if element == 'QUADRUPOLE':
                    addToDB(element_name + "_k1", ui_group[element], ui_unit["_k1"])


    infile.close()


#call main
if __name__ == "__main__":
    main(sys.argv[1:])
