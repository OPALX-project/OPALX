# Matthias Frey
# 14. October 2016, LBNL
# 
# Plot the electric self-field, density and self-field
# potential using the yt framework.
# 1. mpirun -np #cores testSolver
# 2. python visualize.py (make sure you sourced the yt directory $YT_DIR/yt-x86_64/bin/activate)


import os
import yt

opal = os.environ['OPAL_BUILD']
ds = yt.load(opal + "ippl/test/AMR/plt0000", unit_system='mks')

ds.print_stats()


print (ds.field_list)

print (ds.derived_field_list)

slc = yt.SlicePlot(ds, normal='z', fields='rho')
slc.annotate_grids()
slc.save()

slc = yt.SlicePlot(ds, normal='z', fields='Ex')
slc.annotate_grids()
slc.save()

slc = yt.SlicePlot(ds, normal='z', fields='Ey')
slc.annotate_grids()
slc.save()

slc = yt.SlicePlot(ds, normal='x', fields='Ez')
slc.annotate_grids()
slc.save()

slc = yt.SlicePlot(ds, normal='z', fields='potential')
slc.annotate_grids()
slc.save()

ad = ds.all_data()

phi = ad['potential']
print ( phi.max() )
print ( phi.min() )

#dx = ad['dx']
#print ( "dx = ", dx )