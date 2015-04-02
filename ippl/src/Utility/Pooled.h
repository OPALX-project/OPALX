// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef POOLED_H
#define POOLED_H

// include files
#include "Utility/Pool.h"

//////////////////////////////////////////////////////////////////////

template <class T>
class Pooled
{
public:
#ifndef IPPL_DONT_POOL
public:

  // Here is the reason for the class: 
  // Defining new and delete for the classes that inherit from it.

  // Get some memory from the pool.
  inline void* operator new(size_t) { return MyPool.alloc(); }

  // Move a piece of memory back to the pool.
  inline void operator delete(void *p, size_t) { if (p) MyPool.free(p); }

private: 

  static Pool MyPool;	// Storage for objects of type T
#endif
};

//////////////////////////////////////////////////////////////////////

#include "Utility/Pooled.hpp"

#endif // POOLED_H

/***************************************************************************
 * $RCSfile: Pooled.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: Pooled.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
