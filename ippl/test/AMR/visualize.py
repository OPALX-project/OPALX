import yt

ds = yt.load("/home/matthias/Documents/projects/OPAL/build/ippl/test/AMR/plt0000", unit_system='mks')

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

slc = yt.SlicePlot(ds, normal='z', fields='Ez')
slc.annotate_grids()
slc.save()

slc = yt.SlicePlot(ds, normal='z', fields='phi')
slc.annotate_grids()
slc.save()