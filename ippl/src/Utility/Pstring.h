// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef PSTRING_H
#define PSTRING_H

// include files

#if ( defined(IPPL_USE_STANDARD_HEADERS) || \
      defined(IPPL_USE_STANDARD_STRING) )

// use provided string header
#include <string>
// put string class in namespace
#else

#ifdef IPPL_STDSTL
#include <vector>
#else
#include <vector.h>
#endif // IPPL_STDSTL

// indicate that "bool" is a defined type
#define __BOOL_DEFINED 

// include locally modified version of old Modena string header file
#include "Utility/bstring.h"

#endif /* IPPL_USE_STANDARD_HEADERS */


#endif // PSTRING_H

/***************************************************************************
 * $RCSfile: Pstring.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: Pstring.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
