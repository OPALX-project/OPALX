// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

/***********************************************************************

This little gadget lets you allocate an array of elements of type T
and reference count the whole block.  As long as you refer to the block
using objects of the type RefBlockP it will keep the block around.  As
soon as the last RefBlockP that refers to the block is deleted, it
deletes the block.

You create a block like: RefBlockP<T> p(size)
Then you use p like a pointer. Pointer operations are as efficient as 
with a bare pointer as long as bounds checking is off.

It has a second template argument for whether to do bounds checking.
It defaults to the value of the preprocessor symbol
REFBLOCK_BOUNDS_CHECK_DEFAULT.  If you don't define that symbol it 
is set to false, for no bounds checking.

To turn bounds checking on for your whole code define
REFBLOCK_BOUNDS_CHECK_DEFAULT as true.

To turn bounds checking on for a particular RefBlockP you 
declare it like:

RefBlockP<T,true> p;

To turn it off you declare it as:

RefBlockP<T,false> p;

***********************************************************************/

#ifndef REF_BLOCK_H
#define REF_BLOCK_H

// include files
#include "Utility/RefCounted.h"

#include <cassert>
#include <cstddef>

//----------------------------------------------------------------------
#ifndef REFBLOCK_BOUNDS_CHECK_DEFAULT
#  ifndef BOUNDS_CHECK_DEFAULT
#    define REFBLOCK_BOUNDS_CHECK_DEFAULT false
#  else
#    define REFBLOCK_BOUNDS_CHECK_DEFAULT BOUNDS_CHECK_DEFAULT
#  endif
#endif

template<class T, bool BoundsCheck=REFBLOCK_BOUNDS_CHECK_DEFAULT>
class RefBlockP
{
  // 
  // Allow the other bounds checking polarity to look in here.
  // 
// friend class RefBlockP<T,!BoundsCheck>;

private:
  //
  // A simple class that has the pointer to the beginning of the block.
  // Because it is refcounted, when the last reference to the
  // block goes away, this dies and takes the block with it.
  //
  class RefBlockController : public RefCounted
    {
    public: 
      T* P;
      T* Pend;
      bool Dealloc;
      RefBlockController(T *p, size_t size)
	: P(p),Pend(p+size),Dealloc(false)
	  {
	  }
      RefBlockController(size_t size)
	: P(new T[size]),Pend(P+size),Dealloc(true)
	  {
	  }
      ~RefBlockController()
	{
	  if (Dealloc)
	    delete P;
	}
      T* getBlock()
	{
	  return P;
	}
      bool checkDeref(T* p)
	{
	  return (P<=p) && (p<Pend) ;
	}
    };

  //
  // We have two data items in each RefBlockP:
  //    A bare pointer that points somewhere in the block.
  //    A refcounted pointer to the block controller.
  //
  T* P;
  RefCountedP< RefBlockController > Controller;

public:

  //
  // Null ctor for arrays.  Initialize w/ null ptr.
  //
  RefBlockP() : P(0) {}

  //
  // Initialize a block of a given size.
  //
  RefBlockP(size_t size)
    : Controller(new RefBlockController(size))
      {
	P = Controller->getBlock();
      }

  //
  // Initialize with a user allocated pointer.
  // This turns off the deallocation but not the bounds checking.
  //
  RefBlockP(T* p, size_t size)
    : P(p), Controller( new RefBlockController(p,size) )
      {
      }
  
  // Copy ctor and assignment
  RefBlockP(const RefBlockP<T,BoundsCheck>& b)
    : P(b.P), Controller(b.Controller)
      {
      }
  RefBlockP<T,BoundsCheck>& operator=(const RefBlockP<T,BoundsCheck>& rhs)
    {
      P = rhs.P;
      Controller = rhs.Controller;
      return *this;
    }

  // Copy ctor and assignment from a RefBlockP opposite bounds checking.
//  RefBlockP(const RefBlockP<T,!BoundsCheck>& b)
//    : P(b.P), Controller(b.Controller)
//      {
//      }
//  RefBlockP<T,BoundsCheck>& operator=(const RefBlockP<T,!BoundsCheck>& rhs)
//    {
//      P = rhs.P;
//      Controller = rhs.Controller;
//      return *this;
//    }

  //
  // Provide all the usual pointer manipulation functions.
  // Each one makes sure it is legal to do the operation if
  // it is marked for bounds checking.
  // Because BoundsCheck is a compile time parameter, if you don't use it
  // the asserts will be optimized away.
  //
  RefBlockP<T,BoundsCheck>& operator++()
    {
      ++P;
      return *this;
    }
  RefBlockP<T,BoundsCheck> operator++(int)
    {
      RefBlockP<T> save(*this);
      ++P;
      return save;
    }
  RefBlockP<T,BoundsCheck>& operator--()
    {
      --P;
      return *this;
    }
  RefBlockP<T,BoundsCheck> operator--(int)
    {
      RefBlockP<T,BoundsCheck> save(*this);
      --P;
      return save;
    }

  T& operator*() const
    {
      if ( BoundsCheck ) 
	assert( Controller->checkDeref(P) );
      return *P;
    }
  T& operator[](int i) const
    {
      T* p = P+i;
      if ( BoundsCheck )
	assert( Controller->checkDeref(p) );
      return *p;
    }
  T* operator->() const
    {
      if ( BoundsCheck )
	assert( P==0 || Controller->checkDeref(P) );
      return P;
    }

  void operator+=(int i)
    {
      P += i;
    }
  void operator-=(int i)
    {
      P -= i;
    }

  RefBlockP<T,BoundsCheck> operator+(int i)
    {
      RefBlockP<T,BoundsCheck> ret(*this);
      ret += i;
      return ret;
    }
  RefBlockP<T,BoundsCheck> operator-(int i)
    {
      RefBlockP<T,BoundsCheck> ret(*this);
      ret -= i;
      return ret;
    }

  bool operator==(const RefBlockP<T,BoundsCheck>& a) const 
    {
      return P == a.P;
    }
//  bool operator==(const RefBlockP<T,!BoundsCheck>& a) const 
//    {
//      return P == a.P;
//    }

  void invalidate()
    {
      Controller.invalidate();
      P = 0;
    }
  bool valid()
    {
      return P!=0;
    }

  //
  // If something needs a T*, convert it.
  //
  operator T*() const 
    {
      return P;
    }
};

//////////////////////////////////////////////////////////////////////

#endif // REF_BLOCK_H

/***************************************************************************
 * $RCSfile: RefBlock.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: RefBlock.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
