# Matthias Frey
# 14. October 2016, LBNL
# 
# Plot the electric self-field, density and self-field
# potential using the yt framework.
# 1. mpirun -np #cores testSolver
# 2. python visualize.py (make sure you sourced the yt directory $YT_DIR/yt-x86_64/bin/activate)


import os
import yt
import sys

opal = os.environ['OPAL_BUILD']


#ds = yt.load(opal + "ippl/test/AMR/amr0000", unit_system='accel')
#slc = yt.SlicePlot(ds, normal='z', fields='rho')
#slc.annotate_grids()
#slc.save()

nSteps = int(sys.argv[1])

for i in range(0, nSteps):
    
    if i < 10:
        ds = yt.load(opal + "ippl/test/AMR/amr_000" + str(i), unit_system='accel')
    else:
        ds = yt.load(opal + "ippl/test/AMR/amr_00" + str(i), unit_system='accel')
    
    slc = yt.SlicePlot(ds, normal='z', fields='rho')
    #slc.set_unit('rho', 'e/m**3')
    slc.annotate_grids()
    slc.save()

#for i in range(10, 25):
    #ds = yt.load(opal + "ippl/test/AMR/amr_00" + str(i)), unit_system='accel')
    #slc = yt.SlicePlot(ds, normal='z', fields='rho')
    ##slc.set_unit('rho', 'e/m**3')
    #slc.annotate_grids()
    #slc.save()