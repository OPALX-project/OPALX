// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef MY_AUTO_PTR_H
#define MY_AUTO_PTR_H

//////////////////////////////////////////////////////////////////////
/*
  A simple compliant implementation of auto_ptr.
  This is from Greg Colvin's implementation posted to comp.std.c++.

  Instead of using mutable this casts away const in release.

  We have to do this because we can't build containers of these
  things otherwise.
  */
//////////////////////////////////////////////////////////////////////

template<class X>
class my_auto_ptr
{
  X* px;
public:
  my_auto_ptr() : px(0) {}
  my_auto_ptr(X* p) : px(p) {}
  my_auto_ptr(const my_auto_ptr<X>& r) : px(r.release()) {}
  my_auto_ptr& operator=(const my_auto_ptr<X>& r)
  {
    if (&r != this)
      {
	delete px;
	px = r.release();
      }
    return *this;
  }
  ~my_auto_ptr() { delete px; }
  X& operator*()  const { return *px; }
  X* operator->() const { return px; }
  X* get()        const { return px; }
  X* release()    const { X *p=px; ((my_auto_ptr<X>*)(this))->px=0; return p; }
};

#ifdef UNDEFINED
template<class X>
class my_auto_ptr
{
  bool owner;
  X* px;
public:
  my_auto_ptr()
    : owner(false), px(0) {}
  my_auto_ptr(X* p) 
    : owner(p), px(p) {}
  my_auto_ptr(const my_auto_ptr<X>& r) 
    : owner(r.owner), px(r.release()) {}
  my_auto_ptr& operator=(const my_auto_ptr<X>& r)
  {
    if (&r != this)
      {
	if (owner) 
	  delete px;
	owner = r.owner; 
	px = r.release();
      }
    return *this;
  }
  ~my_auto_ptr()        { if (owner) delete px; }
  X& operator*()  const { return *px; }
  X* operator->() const { return px; }
  X* get()        const { return px; }
  X* release()    const { ((my_auto_ptr<X>*)(this))->owner = false; return px;}
};
#endif // UNDEFINED

#endif // MY_AUTO_PTR_H

/***************************************************************************
 * $RCSfile: my_auto_ptr.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:34 $
 * IPPL_VERSION_ID: $Id: my_auto_ptr.h,v 1.1.1.1 2003/01/23 07:40:34 adelmann Exp $ 
 ***************************************************************************/
