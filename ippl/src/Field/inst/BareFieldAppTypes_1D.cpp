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
#include "Field/BareField.h"
#include "Field/Assign.h"
#include "Field/AssignDefs.h"
#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"
#include "AppTypes/AntiSymTenzor.h"

// 1D instantiations
template class BareField<Vektor<IPPL_PRECISION_TYPE,1U>,1U>;
template class BareField<Tenzor<IPPL_PRECISION_TYPE,1U>,1U>;
template class BareField<SymTenzor<IPPL_PRECISION_TYPE,1U>,1U>;
template class BareField<AntiSymTenzor<IPPL_PRECISION_TYPE,1U>,1U>;


/***************************************************************************
 * $RCSfile: BareFieldAppTypes_1D.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: BareFieldAppTypes_1D.cpp,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/

