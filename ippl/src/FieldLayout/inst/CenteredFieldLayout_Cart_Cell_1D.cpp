// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "FieldLayout/CenteredFieldLayout.h"
#include "Meshes/Cartesian.h"
#include "Meshes/Centering.h"

// 1D Cartesian Cell instantiations
template class CenteredFieldLayout<1U,Cartesian<1U,IPPL_PRECISION_TYPE>,Cell>;


/***************************************************************************
 * $RCSfile: CenteredFieldLayout_Cart_Cell_1D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: CenteredFieldLayout_Cart_Cell_1D.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
