// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef SUB_FIELD_ASSIGN_DEFS_H
#define SUB_FIELD_ASSIGN_DEFS_H

// include files
#include "Field/AssignDefs.h"

// forward references
template<class T, unsigned int D, class S> class SubFieldIter;


//////////////////////////////////////////////////////////////////////
//
// Is the domain specification object compressed?
//
//////////////////////////////////////////////////////////////////////

struct DomainCompressed {
  typedef bool PETE_Return_t;
};

template<class T, class S, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,S> &p, DomainCompressed, C)
{
  return p.DomainCompressed();
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,SIndex<D> > &p, DomainCompressed, C)
{
  return p.DomainCompressed();
}

template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,NDIndex<D> > &p, DomainCompressed, C)
{
  return p.DomainCompressed();
}

template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,SOffset<D> > &p, DomainCompressed, C)
{
  return p.DomainCompressed();
}
#endif // __MWERKS__

template<class T, class C, unsigned int D>
inline bool
for_each(typename BareField<T,D>::iterator&, DomainCompressed, C)
{
  return true;
}

template<class C>
inline bool
for_each(Index::cursor&, DomainCompressed, C)
{
  return false;
}

template<class T, class C>
inline bool
for_each(PETE_Scalar<T>&, DomainCompressed, C)
{
  return true;
}


//////////////////////////////////////////////////////////////////////
//
// Do the terms all use the same kind of subset object?
//
//////////////////////////////////////////////////////////////////////

struct SameSubsetType {
  typedef bool PETE_Return_t;
  int fID;
  SameSubsetType(int id) : fID(id) {}
};

template<class T, class S, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,S> &p, SameSubsetType s, C)
{
  return p.matchType(s.fID);
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,SIndex<D> > &p, SameSubsetType s, C)
{
  return p.matchType(s.fID);
}

template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,NDIndex<D> > &p, SameSubsetType s, C)
{
  return p.matchType(s.fID);
}

template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,SOffset<D> > &p, SameSubsetType s, C)
{
  return p.matchType(s.fID);
}
#endif // __MWERKS__

template<class T, class C, unsigned int D>
inline bool
for_each(typename BareField<T,D>::iterator&, SameSubsetType s, C)
{
  return false;
}

template<class C>
inline bool
for_each(Index::cursor&, SameSubsetType, C)
{
  return true;
}

template<class T, class C>
inline bool
for_each(PETE_Scalar<T>&, SameSubsetType, C)
{
  return true;
}


//////////////////////////////////////////////////////////////////////
//
// Initialize all subset objects in an expression before the loop starts
//
//////////////////////////////////////////////////////////////////////

struct SubsetInit {
  typedef int PETE_Return_t;
};

template<class T, class S, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,S> &p, SubsetInit, C) 
{
  p.initialize();
  return 0;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SIndex<D> > &p, SubsetInit, C) 
{
  p.initialize();
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,NDIndex<D> > &p, SubsetInit, C) 
{
  p.initialize();
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SOffset<D> > &p, SubsetInit, C) 
{
  p.initialize();
  return 0;
}
#endif // __MWERKS__

template<class T, class C, unsigned int D>
inline int
for_each(typename BareField<T,D>::iterator &p, SubsetInit, C)
{
  return 0;
}

template<class C>
inline int
for_each(Index::cursor&, SubsetInit, C)
{
  return 0;
}

template<class T, class C>
inline int
for_each(PETE_Scalar<T>&, SubsetInit, C)
{
  return 0;
}


//////////////////////////////////////////////////////////////////////
//
// Set a subfield iterator to point to the next lfield
//
//////////////////////////////////////////////////////////////////////

struct SubsetNextLField {
  typedef int PETE_Return_t;
};

template<class T, class S, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,S> &p, SubsetNextLField, C)
{
  p.nextLField();
  return 0;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SIndex<D> > &p, SubsetNextLField, C)
{
  p.nextLField();
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,NDIndex<D> > &p, SubsetNextLField, C)
{
  p.nextLField();
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SOffset<D> > &p, SubsetNextLField, C)
{
  p.nextLField();
  return 0;
}
#endif // __MWERKS__

template<class T, class C, unsigned int D>
inline int
for_each(typename BareField<T,D>::iterator&, SubsetNextLField, C)
{
  return 0;
}

template<class C>
inline int
for_each(Index::cursor&, SubsetNextLField, C)
{
  return 0;
}

template<class T, class C>
inline int
for_each(PETE_Scalar<T>&, SubsetNextLField, C)
{
  return 0;
}


//////////////////////////////////////////////////////////////////////
//
// Do any of the terms in an expression have an ID equal to a given one?
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,S> &p, SameFieldID s, C)
{
  return p.getBareField().get_Id() == s.fID;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline bool
for_each(const SubFieldIter<T,D,SIndex<D> > &p, const SameFieldID &s, C)
{
  return p.getBareField().get_Id() == s.fID;
}

template<class T, class C, unsigned int D>
inline bool
for_each(const SubFieldIter<T,D,NDIndex<D> > &p, const SameFieldID &s, C)
{
  return p.getBareField().get_Id() == s.fID;
}

template<class T, class C, unsigned int D>
inline bool
for_each(const SubFieldIter<T,D,SOffset<D> > &p, const SameFieldID &s, C)
{
  return p.getBareField().get_Id() == s.fID;
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Plugbase.
//
//////////////////////////////////////////////////////////////////////

template<class T, unsigned D, class S, class C>
inline bool
for_each(SubFieldIter<T,D,S> &p, const PlugBase<D>& f, C)
{
  return p.plugBase(f.Domain);
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned D, class C>
inline bool
for_each(SubFieldIter<T,D,SIndex<D> > &p, const PlugBase<D>& f, C)
{
  return p.plugBase(f.Domain);
}

template<class T, unsigned D, class C>
inline bool
for_each(SubFieldIter<T,D,NDIndex<D> > &p, const PlugBase<D>& f, C)
{
  return p.plugBase(f.Domain);
}

template<class T, unsigned D, class C>
inline bool
for_each(SubFieldIter<T,D,SOffset<D> > &p, const PlugBase<D>& f, C)
{
  return p.plugBase(f.Domain);
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Check for compression.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,S> &p, IsCompressed, C)
{
  return p.IsCompressed();
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,SIndex<D> > &p, IsCompressed, C)
{
  return p.IsCompressed();
}

template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,NDIndex<D> > &p, IsCompressed, C)
{
  return p.IsCompressed();
}

template<class T, class C, unsigned int D>
inline bool
for_each(SubFieldIter<T,D,SOffset<D> > &p, IsCompressed, C)
{
  return p.IsCompressed();
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Evaluation functors.
// First, no arguments.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,S> &p, EvalFunctor_0)
{
  return *p;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SIndex<D> > &p, EvalFunctor_0)
{
  return *p;
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,NDIndex<D> > &p, EvalFunctor_0)
{
  return *p;
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SOffset<D> > &p, EvalFunctor_0)
{
  return *p;
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Evaluation functors.
// One argument.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,S> &p, const EvalFunctor_1& e)
{
  return p.offset(e.I);
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SIndex<D> > &p, const EvalFunctor_1& e)
{
  return p.offset(e.I);
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,NDIndex<D> > &p, const EvalFunctor_1& e)
{
  return p.offset(e.I);
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SOffset<D> > &p, const EvalFunctor_1& e)
{
  return p.offset(e.I);
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Evaluation functors.
// Two arguments.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,S> &p, const EvalFunctor_2& e)
{
  return p.offset(e.I,e.J);
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SIndex<D> > &p, const EvalFunctor_2& e)
{
  return p.offset(e.I,e.J);
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,NDIndex<D> > &p, const EvalFunctor_2& e)
{
  return p.offset(e.I,e.J);
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SOffset<D> > &p, const EvalFunctor_2& e)
{
  return p.offset(e.I,e.J);
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Evaluation functors.
// Three arguments.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,S> &p, const EvalFunctor_3& e)
{
  return p.offset(e.I,e.J,e.K);
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SIndex<D> > &p, const EvalFunctor_3& e)
{
  return p.offset(e.I,e.J,e.K);
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,NDIndex<D> > &p, const EvalFunctor_3& e)
{
  return p.offset(e.I,e.J,e.K);
}

template<class T, unsigned int D>
inline T&
for_each(SubFieldIter<T,D,SOffset<D> > &p, const EvalFunctor_3& e)
{
  return p.offset(e.I,e.J,e.K);
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Step in some dimension.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,S> &p, StepFunctor s, C)
{
  p.step(s.D);
  return 0;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SIndex<D> > &p, StepFunctor s, C)
{
  p.step(s.D);
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,NDIndex<D> > &p, StepFunctor s, C)
{
  p.step(s.D);
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SOffset<D> > &p, StepFunctor s, C)
{
  p.step(s.D);
  return 0;
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Rewind in some dimension.
//
//////////////////////////////////////////////////////////////////////

template<class T, class S, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,S> &p, RewindFunctor s, C)
{
  p.rewind(s.D);
  return 0;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SIndex<D> > &p, RewindFunctor s, C)
{
  p.rewind(s.D);
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,NDIndex<D> > &p, RewindFunctor s, C)
{
  p.rewind(s.D);
  return 0;
}

template<class T, class C, unsigned int D>
inline int
for_each(SubFieldIter<T,D,SOffset<D> > &p, RewindFunctor s, C)
{
  p.rewind(s.D);
  return 0;
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Does an iterator reference something with unit stride?
// Don't worry about it for now.
//
//////////////////////////////////////////////////////////////////////

template<class T, unsigned int D, class S, class C>
inline bool
for_each(SubFieldIter<T,D,S> &p, HasUnitStride, C)
{
  return false;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned int D, class C>
inline bool
for_each(SubFieldIter<T,D,SIndex<D> > &p, HasUnitStride, C)
{
  return false;
}

template<class T, unsigned int D, class C>
inline bool
for_each(SubFieldIter<T,D,NDIndex<D> > &p, HasUnitStride, C)
{
  return false;
}

template<class T, unsigned int D, class C>
inline bool
for_each(SubFieldIter<T,D,SOffset<D> > &p, HasUnitStride, C)
{
  return false;
}
#endif // __MWERKS__

//////////////////////////////////////////////////////////////////////
//
// Ask each term to fill guard cells and compress itself
//
//////////////////////////////////////////////////////////////////////

template<class T, unsigned int D, class S, class C, class T1>
inline int
for_each(SubFieldIter<T,D,S> &p, const FillGCIfNecessaryTag<D,T1> &f, C)
{
  //tjw3/3/99  p.FillGCIfNecessary(f.I, f.I);
  p.FillGCIfNecessary();
  return 0;
}

#ifdef __MWERKS__
// Work around partial-specialization template matching bug in CW4
template<class T, unsigned int D, class C, class T1>
inline int
for_each(SubFieldIter<T,D,SIndex<D> > &p, const FillGCIfNecessaryTag<D,T1> &f, C)
{
  //tjw3/3/99  p.FillGCIfNecessary(f.I, f.I);
  p.FillGCIfNecessary();
  return 0;
}

template<class T, unsigned int D, class C, class T1>
inline int
for_each(SubFieldIter<T,D,NDIndex<D> > &p, const FillGCIfNecessaryTag<D,T1> &f, C)
{
  //tjw3/3/99  p.FillGCIfNecessary(f.I, f.I);
  p.FillGCIfNecessary();
  return 0;
}

template<class T, unsigned int D, class C, class T1>
inline int
for_each(SubFieldIter<T,D,SOffset<D> > &p, const FillGCIfNecessaryTag<D,T1> &f, C)
{
  //tjw3/3/99  p.FillGCIfNecessary(f.I, f.I);
  p.FillGCIfNecessary();
  return 0;
}

#endif // __MWERKS__


#endif // SUB_FIELD_ASSIGN_DEFS_H

/***************************************************************************
 * $RCSfile: SubFieldAssignDefs.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: SubFieldAssignDefs.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
