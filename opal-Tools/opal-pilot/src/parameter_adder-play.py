"""
TODO: fix package issue, run script from root folder meanwhile

Before you execute that script, better get a database backup.                                                                                                                                                                                                                        
"""

I_REALLY_WANT_TO_ADD_STUFF_TO_THE_DATABASE = False

import os
import sys
import _mysql
import _mysql_exceptions

from util.database import DatabaseConnection

if (I_REALLY_WANT_TO_ADD_STUFF_TO_THE_DATABASE):
    parameters = [
        # (parameter, ui_label, ui_tooltip, ui_group_id, ui_unit_id)
        ('FINXB', 'FINXB', 'sdsdsd', 1, 1),
    
    ]
    
    conn = DatabaseConnection()
    
    for tuple in parameters:
        try:
            conn.insert('parameter', ['parameter', 'ui_label', 'ui_tooltip', 'ui_group_id', 'ui_unit_id'], [tuple])
        except _mysql_exceptions.IntegrityError:
            query = """
            UPDATE parameter SET ui_label = '%s', ui_tooltip = '%s', ui_group_id = %s, ui_unit_id = %s WHERE parameter = '%s'
            """ % (tuple[1], tuple[2], tuple[3], tuple[4], tuple[0])
            print query            
#           conn.execute(query)

# -------------------------------------
#
# Quadrupoles
#
# -------------------------------------
baseName=['F10D1-MQUA']

for j in range(len(baseName)):
    for i in range(10, 100, 5):
        print baseName[j]+str(i), baseName[j]+str(i)+'_k1', baseName[j]+str(i)+'_ds', baseName[j]+str(i)+'_dx', baseName[j]+str(i)+'_dy', baseName[j]+str(i)+'_dz'
        
baseName=['FINXB-MQUA']

for j in range(len(baseName)):
    for i in range(10, 50, 10):
        print baseName[j]+str(i), baseName[j]+str(i)+'_k1', baseName[j]+str(i)+'_ds', baseName[j]+str(i)+'_dx', baseName[j]+str(i)+'_dy', baseName[j]+str(i)+'_dz'

# -------------------------------------
#
# Solenoids
#
# -------------------------------------
baseName=['FINSB01_MSOL','FINSB02_MSOL']

for j in range(len(baseName)):
    for i in range(10, 50, 10):
        print baseName[j]+str(i), baseName[j]+str(i)+'_ks', baseName[j]+str(i)+'_ds', baseName[j]+str(i)+'_dx', baseName[j]+str(i)+'_dy', baseName[j]+str(i)+'_dz'
