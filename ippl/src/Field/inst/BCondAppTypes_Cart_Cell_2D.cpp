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
#include "Field/BCond.h"
#include "Meshes/Cartesian.h"
#include "Meshes/Centering.h"
#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"
#include "AppTypes/AntiSymTenzor.h"

// 2D Cartesian Cell instantiations
template class BConds<Vektor<IPPL_PRECISION_TYPE,2U>,2U,
                      Cartesian<2U,IPPL_PRECISION_TYPE>,Cell>;
template class BConds<Tenzor<IPPL_PRECISION_TYPE,2U>,2U,
                      Cartesian<2U,IPPL_PRECISION_TYPE>,Cell>;
template class BConds<SymTenzor<IPPL_PRECISION_TYPE,2U>,2U,
                      Cartesian<2U,IPPL_PRECISION_TYPE>,Cell>;
template class BConds<AntiSymTenzor<IPPL_PRECISION_TYPE,2U>,2U,
                      Cartesian<2U,IPPL_PRECISION_TYPE>,Cell>;


/***************************************************************************
 * $RCSfile: BCondAppTypes_Cart_Cell_2D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: BCondAppTypes_Cart_Cell_2D.cpp,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/
