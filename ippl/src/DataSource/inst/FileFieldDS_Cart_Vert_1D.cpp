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
#include "DataSource/FileFieldDataSource.h"
#include "Meshes/Cartesian.h"
#include "Meshes/Centering.h"
#include "AppTypes/dcomplex.h"

// 1D Cartesian Vert instantiations
template class FileFieldDataSource<
   int, 1U, Cartesian<1U,IPPL_PRECISION_TYPE>, Vert>;

template class FileFieldDataSource<
   IPPL_PRECISION_TYPE, 1U, Cartesian<1U,IPPL_PRECISION_TYPE>, Vert>;

template class FileFieldDataSource<
   dcomplex, 1U, Cartesian<1U,IPPL_PRECISION_TYPE>, Vert>;


/***************************************************************************
 * $RCSfile: FileFieldDS_Cart_Vert_1D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: FileFieldDS_Cart_Vert_1D.cpp,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/

