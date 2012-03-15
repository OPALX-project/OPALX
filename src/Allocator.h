// ------------------------------------------------------------------------
// $RCSfile: Allocator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Allocator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <new>
#include <iostream>
#include <malloc.h>


/// Track "new" and "delete" calls.
//  Class "Allocator" replaces the global operators "new" and "delete".
//  The class has four static counters, each incremented by one of
//  the calls to "::new", "::new[]", "::delete", or "::delete[]".
//  {p}
//  A global object of class "Allocator" must be declared in the main 
//  program.  Its destructor prints the values of these counters.
//  This mechanism accumulates statistics for memory allocation and
//  deallocation.  Note that all global constructors are caught, since
//  the counters are set to zero statically.  However some global
//  destructors may not be caught, when the "Allocator" object goes
//  out of scope before these destructors are called.


struct Allocator {
  static int count1;
  static int count2;
  static int count3;
  static int count4;
  
  /// Constructor.
  //  Do-nothing.
  Allocator()
    {}

  /// Destructor.
  //  Print statistics on allocation and deallocation.
  ~Allocator()
    {
      cerr << "allocated: " << count1 << endl;
      cerr << "freed:     " << count2 << endl;
      cerr << "vector allocated: " << count3 << endl;
      cerr << "vector freed:     " << count4 << endl;
    }
} allocator;
  

int Allocator::count1 = 0;
int Allocator::count2 = 0;
int Allocator::count3 = 0;
int Allocator::count4 = 0;


void *operator new[](std::size_t i) throw(std::bad_alloc)
{
  Allocator::count3++;
  return malloc(i);
}


void operator delete[](void *p) throw()
{
  if (p) {
    Allocator::count4++;
    free(p);
  }
}


void *operator new(std::size_t i) throw(std::bad_alloc)
{
  Allocator::count1++;
  void *p = malloc(i);
  return p;
}


void operator delete(void *p) throw()
{
  if (p) {
    Allocator::count2++;
    free(p);
  }
}
