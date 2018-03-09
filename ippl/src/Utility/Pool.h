// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef POOL_H
#define POOL_H

// include files
#include <cstddef>
#include <vector>
#include <memory.h>
class Pool 
{
private:
  struct Link { Link *next; };	// Point to the next pointer.
  static inline size_t page()       { return 4096-8; }
  static inline int blocks_in_page(size_t sz)
    {
      return (page()>sz)?(page()/sz):1;
    }

public: 

  static inline size_t log2_align() { return 3; }
  static inline size_t align()      { return (1<<log2_align()); }
  static inline size_t align_mask() { return align()-1; }

  static inline size_t round_to_align(size_t s)
    {
      if (s)
	s = (s & ~align_mask()) + ((s&align_mask())?align():0);
      else
	s = align();
      return s;
    }


  // Null ctor creates an invalid Pool.
  Pool();

  // Make a new Pool for objects of a given size.
  Pool(size_t);

  // Destroy the Pool.
  ~Pool();

  // Allocate a block from the pool.
  inline void* alloc()
  {
    ++outstandingAllocs;        // Record an alloc.
    if ( head==0 ) grow();	// If we're out of space, get more.
    Link *p = head;		// get the top of the stack.

    // Make the next one the new head of the list.
    // We can't do head = p->next since p will soon be treated
    // as something other than a Link.  By doing this assignment
    // with memcpy, we ensure that p->next will be read before
    // it is clobbered. 
    memcpy(&head, &p->next, sizeof(head));
      
    // Return the requested block.
    return p;
  }

  // Release a block to the pool.
  inline void free(void *b)
  {
    --outstandingAllocs;        // Record a free.
    Link *p = (Link*)b;		// Cast this to a Link
    p->next = head;		// Makew it point to top of stack.
    head = p;			// Make it the new top.
  }

private:
  int outstandingAllocs;        // Number of blocks given to user.
  size_t bsize;			// How big is each block.
  size_t nblock;		// How many to allocate at once.
  Link *head;			// The first one.
  std::vector<char*> chunks;         // Currently allocated chunks
  void grow();                  // Allocate another chunk and put
                                // its blocks in the free list
};

//////////////////////////////////////////////////////////////////////

#endif // POOL_H

/***************************************************************************
 * $RCSfile: Pool.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: Pool.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
