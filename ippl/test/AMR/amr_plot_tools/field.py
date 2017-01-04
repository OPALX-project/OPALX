##
# @file field.py
# @author Matthias Frey
# @date 4. Jan. 2017
# @version 1.0
# @brief Plot vector and scalar fields
# @details When configuring OPAL with the flag
#          -DDBG_SCALARFIELD=1 it dumps the
#          scalar and vector fields of every
#          time step to "data" directory of the
#          simulation. These are visualized by
#          this script.\n
#          Call: python field.py [file] \n
#          where [file] either *field*
#          or *scalar*
#

from tools import match, doGridPlot
import sys
import numpy as np

# regular expression pattern, where parentheses () represent groups that can be accessed later
# group(0) is whole match

# format: i j k ( x-component , y-component , z-component )
# where i, j and k are the grid points
vector_field_pattern = r'(\d*) (\d*) (\d*) \( (-*.*\d*) \, (-*.*\d*) \, (-*.*\d*) \)'

# format: i j k value
# where i, j and k are the grid points
scalar_field_pattern = r'(\d*) (\d*) (\d*) (-*.*\d*)'


try:
    
    f = sys.argv[1]
    
    if 'field' in f:
        data = match(f, vector_field_pattern)
        
        # electric field in horizontal direction
        mid = int(0.5 * max(data[:, 2]))
        i = np.extract(data[:, 2] == mid, data[:, 0])
        j = np.extract(data[:, 2] == mid, data[:, 1])
        val = np.extract(data[:, 2] == mid, data[:, 3])
        plt = doGridPlot(i, j, val, "horizontal grid", "vertical grid", r'$E_x$ [V/m]')
        plt.savefig(f + '_slice_z_Ex.png')#, bbox_inches='tight')
        
        print ("max. Ex: ", max(val))
        print ("min. Ex: ", min(val))
        
        
        # electric field in vertical direction
        val = np.extract(data[:, 2] == mid, data[:, 4])
        plt = doGridPlot(i, j, val, "horizontal grid", "vertical grid", r'$E_y$ [V/m]')
        plt.savefig(f + '_slice_z_Ey.png')#bbox_inches='tight')
        
        print ("max. Ey: ", max(val))
        print ("min. Ey: ", min(val))
        
        # center in z
        mid = int(0.5 * max(data[:, 0]))
        i = np.extract(data[:, 0] == mid, data[:, 1])
        j = np.extract(data[:, 0] == mid, data[:, 2])
        val = np.extract(data[:, 0] == mid, data[:, 5])
        plt = doGridPlot(i, j, val, "vertical grid", "longitudinal grid", r'$E_z$ [V/m]')
        plt.savefig(f + '_slice_x_Ez.png')# bbox_inches='tight')
        
        print ("max. Ez: ", max(val))
        print ("min. Ez: ", min(val))
        
    elif 'scalar' in f:
        data = match(f, scalar_field_pattern)
        
        clab = ''
        if 'phi' in f:
            clab = r'$\phi$ [V]'
        elif 'rho' in f:
            clab = r'$\rho$ [C/m**3]'
        
        # center in x
        mid = int(0.5 * max(data[:, 0]))
        i = np.extract(data[:, 0] == mid, data[:, 1])
        j = np.extract(data[:, 0] == mid, data[:, 2])
        val = np.extract(data[:, 0] == mid, data[:, 3])
        plt = doGridPlot(i, j, val, "vertical grid", "longitudinal grid", clab)
        plt.savefig(f + '_slice_x.png')#, bbox_inches='tight')
        
        # center in y
        mid = int(0.5 * max(data[:, 1]))
        i = np.extract(data[:, 1] == mid, data[:, 0])
        j = np.extract(data[:, 1] == mid, data[:, 2])
        val = np.extract(data[:, 1] == mid, data[:, 3])
        plt = doGridPlot(i, j, val, "horizontal grid", "longitudinal grid", clab)
        plt.savefig(f + '_slice_y.png')#bbox_inches='tight')
        
        # center in z
        mid = int(0.5 * max(data[:, 2]))
        i = np.extract(data[:, 2] == mid, data[:, 0])
        j = np.extract(data[:, 2] == mid, data[:, 1])
        val = np.extract(data[:, 2] == mid, data[:, 3])
        plt = doGridPlot(i, j, val, "horizontal grid", "vertical grid", clab)
        plt.savefig(f + '_slice_z.png')# bbox_inches='tight')
        
        print ("max. value (slice z): ", max(val))
        print ("min. value (slice z): ", min(val))
        
    else:
        raise RuntimeError('This type of file is not supported.')
    
except Exception as e:
    print (e.strerror)