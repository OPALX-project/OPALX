// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef RANDOM_NUMBER_GEN_H
#define RANDOM_NUMBER_GEN_H

/***********************************************************************
 *
 * class RandomNumberGen
 *
 * RandomNumberGen is actually just a typedef for one of the many
 * RNGSequence classes available in IPPL.  It is selected by examining the
 * values of #define's.  Other RNGSequence classes may be used in a program,
 * but using RandomNumberGen give you the 'default' RNGSequence type as
 * selected at compile time.
 *
 * When using Random or Distributed number sequences, include this file.
 * 
 ***********************************************************************/

// include files
#include "Utility/RNGRand.h"
#include "Utility/RNGSimple.h"
#include "Utility/RNGXDiv.h"
#include "Utility/RNGXCI.h"

// typedef for the default RNG ... selected by #define
#if   defined(IPPL_USE_XDIV_RNG)
typedef RNGXDivSequence RandomNumberGen;
#elif defined(IPPL_USE_RAND_RNG)
typedef RNGRandSequence RandomNumberGen;
#elif defined(IPPL_USE_SIMPLE_RNG)
typedef RNGSimpleSequence RandomNumberGen;
#elif defined(IPPL_USE_XCI_RNG)
typedef RNGXCISequence RandomNumberGen;
#else
typedef RNGSimpleSequence RandomNumberGen;
#endif


// a default RandomNumberGen object for use in the Framework.  When
// running in parallel, the Ippl object will advance this by the
// node number so as to have different RNG sequences on each node.  If the
// same RNG sequence is needed on each node, the user must instantiate their
// own RNG sequence object and use that.
extern RandomNumberGen IpplRandom;


#endif // RANDOM_NUMBER_GEN_H

/***************************************************************************
 * $RCSfile: RandomNumberGen.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: RandomNumberGen.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
