##
# @file tools.py
# @author Matthias Frey
# @date 22. Dec. 2016
# @version 1.0
# @brief This module contains several helper functions
#        for data analysis
#

import yt

##
# @param ds is the data
# @param direct is the direction 'x', 'y' or 'z' (normal)
# @param field to plot
# @param unit the data should be converted to (otherwise it
#        takes the default given by the data)
# @param col is the color for the time stamp and scale annotation
def doSlicePlot(ds, direct, field, unit, col = 'white'):
    slc = yt.SlicePlot(ds, normal=direct, fields=field)
        
    if unit is not None:
        slc.set_unit(field, unit)
        
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
# @param col is the color for the time stamp and scale annotation
def doProjectionPlot(ds, direct, field, unit, col = 'white'):
    slc = yt.ProjectionPlot(ds, direct, fields=field)
        
    if unit is not None:
        slc.set_unit(field, unit)
        
    slc.annotate_grids()
    slc.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    slc.annotate_scale(corner='upper_right', size_bar_args={'color':col})
    slc.save()


##
# Take an integer and transform it
# to a string of for characters, e.g.\n
# \f$ 1 \rightarrow 0001 \f$ \n
# \f$ 12 \rightarrow 0012 \f$ \n
# \f$ 586 \rightarrow 0586 \f$ \n
# @param step is an integer
def concatenate(step):
    res = str(i)
    while len(res) < 4:
        res = '0' + res
    return res