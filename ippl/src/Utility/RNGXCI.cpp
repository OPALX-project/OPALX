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
 * Visit www.amas.web.psi for more details
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
#include "Utility/RNGXCI.h"

// initialize static variables for RNGXCI

const RNlong RNGXCI::RN_MULT   =  19073486328125LL;        // 5^19 
const RNlong RNGXCI::RN_MOD    = 281474976710656LL;        // 2^48 
const RNlong RNGXCI::RN_PERIOD = RNGXCI::RN_MOD/4; // period
// normalize to (0,1)
const double RNGXCI::RN_NORM   = 1./281474976710656.;
// 48-bit mask
const RNlong RNGXCI::RN_MASK   = RNGXCI::RN_MOD - 1L;


/***************************************************************************
 * $RCSfile: RNGXCI.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: RNGXCI.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
