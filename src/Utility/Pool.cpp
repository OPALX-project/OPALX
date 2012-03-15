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
#include "Utility/Pool.h"
#include "Utility/PAssert.h"
#include "Profile/Profiler.h"

#include <cstdlib>

// Null ctor creates an invalid Pool.
Pool::Pool()
  : outstandingAllocs(0),          // Number of blocks given to user
    bsize(0),                      // The size of each block
    nblock(0),                     // Number of blocks
    head(0)                        // the first one.
{
}

// Make a new Pool for objects of a given size.
Pool::Pool(size_t sz)
  : outstandingAllocs(0),          // Number of blocks given to user
    bsize(round_to_align(sz)),     // The size of each block
    nblock(blocks_in_page(bsize)), // Number of blocks
    head(0)                        // the first one.
{
}

// Destroy the Pool.
Pool::~Pool()
{
  PInsist(outstandingAllocs==0,"Not all of the pooled memory was freed!");

  // Loop over the allocated chunks and free them.
  std::vector<char*>::iterator p, pend = chunks.end();
  for (p = chunks.begin(); p != pend; ++p)
    delete [] *p;
}

// Grow a Pool.
void Pool::grow()
{
  TAU_PROFILE("Pool::grow()", "void ()", TAU_UTILITY);
  size_t alloc_this;
  if ( bsize>page() )
    alloc_this = bsize;
  else
    alloc_this = page();

  char *start = new char[alloc_this];     // Allocate aligned space.
  chunks.push_back(start);                // Put it in list of chunks.
  char *last  = start + (nblock-1)*bsize; // Get a pointer to the last one.

  for (char *p=start; p!=last; p+=bsize)  // For all but the last one
    ((Link*)p)->next = (Link*)(p+bsize);  //   point to the next.
  ((Link*)last)->next = head;		  // The last points to the head
  head = (Link*)start;			  // Reset the head to the first.
}

/***************************************************************************
 * $RCSfile: Pool.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: Pool.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
