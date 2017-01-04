##
# @file tracking_vis.py
# @author Matthias Frey
# @date 14. October 2016, LBNL
# @version 1.1 (22. Dec. 2016)
# @brief Plot the electric self-field, density and self-field
# potential using the yt framework.
# @details It reads in the all generated BoxLib plotfiles "plt_amr_*"
# (* are four digits) of a directory and saves the content as slice
# plots. It reads in the output of testReal.cpp.
# Make sure you sourced the yt directory $YT_DIR/yt-x86_64/bin/activate.
#

import os
import sys
import yt

from tools import countSubdirs, concatenate, doSlicePlot

try:
    opal = os.environ['OPAL_BUILD']
    
    # leading name of plotfiles
    substr = 'plt_amr_'
    
    nFiles = countSubdirs(opal + '/ippl/test/AMR/', substr)
    
    print ('Found ' + str(nFiles) + ' plotfiles')
    
    for i in range(0, nFiles):
        
        # do a string of four characters, e.g. 12 --> 0012
        res = concatenate(i)
        
        ds = yt.load(opal + '/ippl/test/AMR/' + substr + res, unit_system='mks')
        
        doSlicePlot(ds, 'z', 'rho', 'C/m**3', 'gray')
        
        doSlicePlot(ds, 'z', 'potential', 'V')
        
        doSlicePlot(ds, 'z', 'Ex', 'V/m')
        
        doSlicePlot(ds, 'z', 'Ey', 'V/m')
        
        doSlicePlot(ds, 'x', 'Ez', 'V/m')

except KeyError:
    print ("Please export the environment variable 'OPAL_BUILD'.")
except IOError as e:
    print (e.strerror)
except TypeError:
    print ("No subdirectories containing the given substring " +
           "found in " + opal + "'/ippl/test/AMR/'.")