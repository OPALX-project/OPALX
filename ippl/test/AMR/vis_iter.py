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


# potential
i, j, k, phi = np.loadtxt(opal + "/ippl/test/AMR/phi.dat", unpack=True)

h = int(0.5 * max(k))

xi = np.extract(k == h, i)
yi = np.extract(k == h, j)
phi = np.extract(k == h, phi)

print max(phi)
print min(phi)



doPlot(xi, yi, phi, 'grid in x', 'grid in y', r'$\phi\quad [V]$', 221, 'Electrostatic Potential',
       0, max(xi)+1, 0, max(yi)+1)

# horizontal electric field
i, j, ex = np.loadtxt(opal + "/ippl/test/AMR/ex.dat", unpack=True)

doPlot(i, j, ex, 'grid in x', 'grid in y', r'$E_{x}\quad [V/m]$', 222, 'Horizontal Electric Field',
       -1, max(i)+1, -1, max(j)+1)


print 'max Ex = ', max(ex)
print 'min Ex = ', min(ex)

# vertical electric field
i, j, ey = np.loadtxt(opal + "/ippl/test/AMR/ey.dat", unpack=True)

doPlot(i, j, ey, 'grid in x', 'grid in y', r'$E_{y}\quad [V/m]$', 223, 'Vertical Electric Field',
       -1, max(i)+1, -1, max(j)+1)

print 'max Ey = ', max(ey)
print 'min Ey = ', min(ey)


# longitudinal electric field
j, k, ez = np.loadtxt(opal + "/ippl/test/AMR/ez.dat", unpack=True)

doPlot(j, k, ez, 'grid in y', 'grid in z', r'$E_{z}\quad [V/m]$', 224, 'Longitudinal Electric Field',
       -1, max(j)+1, -1, max(k)+1)

print 'max Ez = ', max(ez)
print 'min Ez = ', min(ez)

# 21. Dec. 2016
# http://stackoverflow.com/questions/6541123/improve-subplot-size-spacing-with-many-subplots-in-matplotlib
left  = 0.05
right = 0.95
bottom = 0.05
top = 0.95
wspace = 0.1
hspace = 0.15
plt.subplots_adjust(left, bottom, right, top, wspace, hspace)

plt.show()