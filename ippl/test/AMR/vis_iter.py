##
# @file vis_iter.py
# @author Matthias Frey
# @date 19. Oct. 2016, LBNL
#
# @pre Environment variable OPAL_BUILD has to be set.
# @details Script for visualizing the output files of iterative.cpp.
#          It plots the potential (centered in longitudinal direction)
#          and the elecric field components (centered)
# @brief Plot potential and electric field of iterative.cpp

import matplotlib.pyplot as plt
import numpy as np
import os

def doPlot(var1, var2, col, xlab, ylab, clab, sub, tit,
           xmin, xmax, ymin, ymax):
    plt.subplot(sub)
    plt.subplots_adjust(hspace=0.4, wspace=0.5)
    plt.scatter(var1, var2, c=col, marker=',', s=100, edgecolors='face')
    plt.title(tit)
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    cb = plt.colorbar()
    cb.set_label(clab)


opal = os.environ['OPAL_BUILD']


# density
i, j, k, phi = np.loadtxt(opal + "/ippl/test/AMR/phi.dat", unpack=True)

h = int(0.5 * max(k))

xi = np.extract(k == h, i)
yi = np.extract(k == h, j)
phi = np.extract(k == h, phi)

print max(phi)
print min(phi)



doPlot(xi, yi, phi, 'grid in x', 'grid in y', r'$\phi$', 221, 'Density Plot',
       0, max(xi)+1, 0, max(yi)+1)

# horizontal electric field
i, j, ex = np.loadtxt(opal + "/ippl/test/AMR/ex.dat", unpack=True)

doPlot(i, j, ex, 'grid in x', 'grid in y', r'$E_{x}$', 222, 'Horizontal Electric Field',
       -1, max(i)+1, -1, max(j)+1)

# vertical electric field
i, j, ey = np.loadtxt(opal + "/ippl/test/AMR/ey.dat", unpack=True)

doPlot(i, j, ey, 'grid in x', 'grid in y', r'$E_{y}$', 223, 'Vertical Electric Field',
       -1, max(i)+1, -1, max(j)+1)


# longitudinal electric field
j, k, ez = np.loadtxt(opal + "/ippl/test/AMR/ez.dat", unpack=True)

doPlot(j, k, ez, 'grid in y', 'grid in z', r'$E_{z}$', 224, 'Longitudinal Electric Field',
       -1, max(j)+1, -1, max(k)+1)

plt.show()