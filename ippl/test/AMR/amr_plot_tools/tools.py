##
# @file tools.py
# @author Matthias Frey
# @date 22. Dec. 2016
# @version 1.0
# @brief This module contains several helper functions
#        for data analysis
#

import yt
import os
import re
import numpy as np
import matplotlib.pyplot as plt

##
# @param ds is the data
# @param direct is the direction 'x', 'y' or 'z' (normal)
# @param field to plot
# @param unit the data should be converted to (otherwise it
#        takes the default given by the data)
# @param zfactor is the zoom factor (default: 1, i.e. no zoom)
# @param col is the color for the time stamp and scale annotation
def doSlicePlot(ds, direct, field, unit, zfactor = 1, col = 'white'):
    slc = yt.SlicePlot(ds, normal=direct, fields=field)
        
    if unit is not None:
        slc.set_unit(field, unit)
    
    slc.zoom(zfactor)
    slc.annotate_grids()
    slc.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    slc.annotate_scale(corner='upper_right', size_bar_args={'color':col})
    slc.save()
    
##
# @param ds is the data
# @param direct is the direction 'x', 'y' or 'z' (normal)
# @param field to plot
# @param unit the data should be converted to (otherwise it
#        takes the default given by the data)
# @param zfactor is the zoom factor (default: 1, i.e. no zoom)
# @param col is the color for the time stamp and scale annotation
def doProjectionPlot(ds, direct, field, unit, zfactor = 1, col = 'white'):
    slc = yt.ProjectionPlot(ds, direct, fields=field)
        
    if unit is not None:
        slc.set_unit(field, unit)
    
    slc.zoom(zfactor)
    slc.annotate_grids()
    slc.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    slc.annotate_scale(corner='upper_right', size_bar_args={'color':col})
    slc.save()


##
# Take an integer and transform it
# to a string of four characters, e.g.\n
# \f$ 1 \rightarrow 0001 \f$ \n
# \f$ 12 \rightarrow 0012 \f$ \n
# \f$ 586 \rightarrow 0586 \f$ \n
# @param step is an integer
def concatenate(step):
    res = str(step)
    while len(res) < 4:
        res = '0' + res
    return res

##
# Count subdirectories
# @param parent is the path to the parent directory
# @param substr specifies the substring that should be contained
# @returns the number of subdirectories containing
# a given substring
def countSubdirs(parent, substr):
    nDirs = 0
    
    # 22. Dec. 2016
    # http://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
    for filename in os.listdir(parent):
        if substr in filename:
            nDirs = nDirs + 1
    return nDirs

##
# Read field data written by OPAL (-DDBG_SCALARFIELD=1) using
# regular expression
# @param f is the file that is read line by line
# @param pattern is a string specifying rule for matching
# @returns a matrix where each row is matched line in the file
def match(f, pattern):
    
    regex = re.compile(pattern)
    
    data = np.empty((0, regex.groups))
    
    # go through file line by line and match
    with open(f) as ff:
        for line in ff:
            res = regex.match(line)
            row = []
            for i in range(1, regex.groups + 1):
                row = np.append(row, [float(res.group(i))])
            data = np.append(data, np.array([row]), axis=0)
    return data

##
def doGridPlot(i, j, val, xlab, ylab, clab):
    imin = int(i[0])
    imax = int(i[-1])
    
    jmin = int(j[0])
    jmax = int(j[-1])
    
    i = np.reshape(i, (imax, jmax))
    j = np.reshape(j, (imax, jmax))
    val = np.reshape(val, (imax, jmax))
    
    plt.figure()
    plt.pcolor(i, j, val, cmap='YlGnBu')
    plt.xlim(imin, imax)
    plt.xlabel(xlab)
    plt.ylim(imin, jmax)
    plt.ylabel(ylab)
    plt.colorbar(label=clab)
    return plt
