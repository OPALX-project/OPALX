// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef INT_CIC_H
#define INT_CIC_H

/* IntCIC.h -- Definition of simple class to perform cloud-in-cell
   interpolation of data for a single particle to or from a IPPL Field.
   This interpolation scheme is also referred to as linear interpolation,
   area-weighting (in 2D), or volume-weighting (in 3D).                    */

// include files
#include "Particle/Interpolator.h"
#include "Field/Field.h"


// forward declaration
class IntCIC;

// specialization of InterpolatorTraits
template <class T, unsigned Dim>
struct InterpolatorTraits<T,Dim,IntCIC> {
  typedef CacheDataCIC<T,Dim> Cache_t;
};



// IntCICImpl class definition
template <unsigned Dim>
class IntCICImpl : public Interpolator {

public:
  // constructor/destructor
  IntCICImpl() {}
  ~IntCICImpl() {}

  // gather/scatter functions

  // scatter particle data into Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,Dim,M,C>& f,
	       const Vektor<PT,Dim>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,Dim> gpos, dpos, delta;
    NDIndex<Dim> ngp;
    CompressedBrickIterator<FT,Dim> fiter;
    int lgpoff[Dim];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<Dim; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    ERRORMSG("IntCIC::scatter: not implemented for Dim>3!!"<<endl);
    return;
  }

  // scatter particle data into Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,Dim,M,C>& f,
	       const Vektor<PT,Dim>& ppos, const M& mesh,
               NDIndex<Dim>& ngp, int lgpoff[Dim], Vektor<PT,Dim>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,Dim> gpos, delta;
    CompressedBrickIterator<FT,Dim> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<Dim; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    ERRORMSG("IntCIC::scatter: not implemented for Dim>3!!"<<endl);
    return;
  }

  // scatter particle data into Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,Dim,M,C>& f,
	       const NDIndex<Dim>& ngp, const int lgpoff[Dim],
               const Vektor<PT,Dim>& dpos) {
    CompressedBrickIterator<FT,Dim> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into local elements
    ERRORMSG("IntCIC::scatter: not implemented for Dim>3!!"<<endl);
    return;
  }

  // gather particle data from Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,Dim,M,C>& f,
	      const Vektor<PT,Dim>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,Dim> gpos, dpos, delta;
    NDIndex<Dim> ngp;
    CompressedBrickIterator<FT,Dim> fiter;
    int lgpoff[Dim];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<Dim; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    ERRORMSG("IntCIC::gather: not implemented for Dim>3!!"<<endl);
    return;
  }

  // gather particle data from Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,Dim,M,C>& f,
	      const Vektor<PT,Dim>& ppos, const M& mesh,
              NDIndex<Dim>& ngp, int lgpoff[Dim], Vektor<PT,Dim>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,Dim> gpos, delta;
    CompressedBrickIterator<FT,Dim> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<Dim; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    ERRORMSG("IntCIC::gather: not implemented for Dim>3!!"<<endl);
    return;
  }

  // gather particle data from Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,Dim,M,C>& f,
	      const NDIndex<Dim>& ngp, const int lgpoff[Dim],
              const Vektor<PT,Dim>& dpos) {
    CompressedBrickIterator<FT,Dim> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into particle attrib
    ERRORMSG("IntCIC::gather: not implemented for Dim>3!!"<<endl);
    return;
  }

};


template <>
class IntCICImpl<1U> : public Interpolator {

public:
  // constructor/destructor
  IntCICImpl() {}
  ~IntCICImpl() {}

  // gather/scatter functions

  // scatter particle data into Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,1U,M,C>& f,
	       const Vektor<PT,1U>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,1U> gpos, dpos, delta;
    NDIndex<1U> ngp;
    CompressedBrickIterator<FT,1U> fiter;
    int lgpoff[1U];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<1U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * pdata;
    fiter.offset(1) += dpos(0) * pdata;
    return;
  }

  // scatter particle data into Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,1U,M,C>& f,
	       const Vektor<PT,1U>& ppos, const M& mesh,
               NDIndex<1U>& ngp, int lgpoff[1U], Vektor<PT,1U>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,1U> gpos, delta;
    CompressedBrickIterator<FT,1U> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<1U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * pdata;
    fiter.offset(1) += dpos(0) * pdata;
    return;
  }

  // scatter particle data into Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,1U,M,C>& f,
	       const NDIndex<1U>& ngp, const int lgpoff[1U],
               const Vektor<PT,1U>& dpos) {
    CompressedBrickIterator<FT,1U> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * pdata;
    fiter.offset(1) += dpos(0) * pdata;
    return;
  }

  // gather particle data from Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,1U,M,C>& f,
	      const Vektor<PT,1U>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,1U> gpos, dpos, delta;
    NDIndex<1U> ngp;
    CompressedBrickIterator<FT,1U> fiter;
    int lgpoff[1U];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<1U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (*fiter) +
            dpos(0)       * fiter.offset(1);
    return;
  }

  // gather particle data from Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,1U,M,C>& f,
	      const Vektor<PT,1U>& ppos, const M& mesh,
              NDIndex<1U>& ngp, int lgpoff[1U], Vektor<PT,1U>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,1U> gpos, delta;
    CompressedBrickIterator<FT,1U> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<1U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (*fiter) +
            dpos(0)       * fiter.offset(1);
    return;
  }

  // gather particle data from Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,1U,M,C>& f,
	      const NDIndex<1U>& ngp, const int lgpoff[1U],
              const Vektor<PT,1U>& dpos) {
    CompressedBrickIterator<FT,1U> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (*fiter) +
            dpos(0)       * fiter.offset(1);
    return;
  }

};


template <>
class IntCICImpl<2U> : public Interpolator {

public:
  // constructor/destructor
  IntCICImpl() {}
  ~IntCICImpl() {}

  // gather/scatter functions

  // scatter particle data into Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,2U,M,C>& f,
	       const Vektor<PT,2U>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,2U> gpos, dpos, delta;
    NDIndex<2U> ngp;
    CompressedBrickIterator<FT,2U> fiter;
    int lgpoff[2U];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<2U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * (1 - dpos(1)) * pdata;
    fiter.offset(1,0) += dpos(0) * (1 - dpos(1)) * pdata;
    fiter.offset(0,1) += (1 - dpos(0)) * dpos(1) * pdata;
    fiter.offset(1,1) += dpos(0) * dpos(1) * pdata;
    return;
  }

  // scatter particle data into Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,2U,M,C>& f,
	       const Vektor<PT,2U>& ppos, const M& mesh,
               NDIndex<2U>& ngp, int lgpoff[2U], Vektor<PT,2U>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,2U> gpos, delta;
    CompressedBrickIterator<FT,2U> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<2U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * (1 - dpos(1)) * pdata;
    fiter.offset(1,0) += dpos(0) * (1 - dpos(1)) * pdata;
    fiter.offset(0,1) += (1 - dpos(0)) * dpos(1) * pdata;
    fiter.offset(1,1) += dpos(0) * dpos(1) * pdata;
    return;
  }

  // scatter particle data into Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,2U,M,C>& f,
	       const NDIndex<2U>& ngp, const int lgpoff[2U],
               const Vektor<PT,2U>& dpos) {
    CompressedBrickIterator<FT,2U> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * (1 - dpos(1)) * pdata;
    fiter.offset(1,0) += dpos(0) * (1 - dpos(1)) * pdata;
    fiter.offset(0,1) += (1 - dpos(0)) * dpos(1) * pdata;
    fiter.offset(1,1) += dpos(0) * dpos(1) * pdata;
    return;
  }

  // gather particle data from Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,2U,M,C>& f,
	      const Vektor<PT,2U>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,2U> gpos, dpos, delta;
    NDIndex<2U> ngp;
    CompressedBrickIterator<FT,2U> fiter;
    int lgpoff[2U];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<2U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (1 - dpos(1)) * (*fiter) +
            dpos(0)       * (1 - dpos(1)) * fiter.offset(1,0) +
            (1 - dpos(0)) * dpos(1)       * fiter.offset(0,1) +
            dpos(0)       * dpos(1)       * fiter.offset(1,1);
    return;
  }

  // gather particle data from Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,2U,M,C>& f,
	      const Vektor<PT,2U>& ppos, const M& mesh,
              NDIndex<2U>& ngp, int lgpoff[2U], Vektor<PT,2U>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,2U> gpos, delta;
    CompressedBrickIterator<FT,2U> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<2U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (1 - dpos(1)) * (*fiter) +
            dpos(0)       * (1 - dpos(1)) * fiter.offset(1,0) +
            (1 - dpos(0)) * dpos(1)       * fiter.offset(0,1) +
            dpos(0)       * dpos(1)       * fiter.offset(1,1);
    return;
  }

  // gather particle data from Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,2U,M,C>& f,
	      const NDIndex<2U>& ngp, const int lgpoff[2U],
              const Vektor<PT,2U>& dpos) {
    CompressedBrickIterator<FT,2U> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (1 - dpos(1)) * (*fiter) +
            dpos(0)       * (1 - dpos(1)) * fiter.offset(1,0) +
            (1 - dpos(0)) * dpos(1)       * fiter.offset(0,1) +
            dpos(0)       * dpos(1)       * fiter.offset(1,1);
    return;
  }

};


template <>
class IntCICImpl<3U> : public Interpolator {

public:
  // constructor/destructor
  IntCICImpl() {}
  ~IntCICImpl() {}

  // gather/scatter functions

  // scatter particle data into Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,3U,M,C>& f,
	       const Vektor<PT,3U>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,3U> gpos, dpos, delta;
    NDIndex<3U> ngp;
    CompressedBrickIterator<FT,3U> fiter;
    int lgpoff[3U];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<3U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * (1 - dpos(1)) * (1 - dpos(2)) * pdata;
    fiter.offset(1,0,0) += dpos(0) * (1 - dpos(1)) * (1 - dpos(2)) * pdata;
    fiter.offset(0,1,0) += (1 - dpos(0)) * dpos(1) * (1 - dpos(2)) * pdata;
    fiter.offset(1,1,0) += dpos(0) * dpos(1) * (1 - dpos(2)) * pdata;
    fiter.offset(0,0,1) += (1 - dpos(0)) * (1 - dpos(1)) * dpos(2) * pdata;
    fiter.offset(1,0,1) += dpos(0) * (1 - dpos(1)) * dpos(2) * pdata;
    fiter.offset(0,1,1) += (1 - dpos(0)) * dpos(1) * dpos(2) * pdata;
    fiter.offset(1,1,1) += dpos(0) * dpos(1) * dpos(2) * pdata;
    return;
  }

  // scatter particle data into Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,3U,M,C>& f,
	       const Vektor<PT,3U>& ppos, const M& mesh,
               NDIndex<3U>& ngp, int lgpoff[3U], Vektor<PT,3U>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,3U> gpos, delta;
    CompressedBrickIterator<FT,3U> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<3U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * (1 - dpos(1)) * (1 - dpos(2)) * pdata;
    fiter.offset(1,0,0) += dpos(0) * (1 - dpos(1)) * (1 - dpos(2)) * pdata;
    fiter.offset(0,1,0) += (1 - dpos(0)) * dpos(1) * (1 - dpos(2)) * pdata;
    fiter.offset(1,1,0) += dpos(0) * dpos(1) * (1 - dpos(2)) * pdata;
    fiter.offset(0,0,1) += (1 - dpos(0)) * (1 - dpos(1)) * dpos(2) * pdata;
    fiter.offset(1,0,1) += dpos(0) * (1 - dpos(1)) * dpos(2) * pdata;
    fiter.offset(0,1,1) += (1 - dpos(0)) * dpos(1) * dpos(2) * pdata;
    fiter.offset(1,1,1) += dpos(0) * dpos(1) * dpos(2) * pdata;
    return;
  }

  // scatter particle data into Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,3U,M,C>& f,
	       const NDIndex<3U>& ngp, const int lgpoff[3U],
               const Vektor<PT,3U>& dpos) {
    CompressedBrickIterator<FT,3U> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into local elements
    *fiter += (1 - dpos(0)) * (1 - dpos(1)) * (1 - dpos(2)) * pdata;
    fiter.offset(1,0,0) += dpos(0) * (1 - dpos(1)) * (1 - dpos(2)) * pdata;
    fiter.offset(0,1,0) += (1 - dpos(0)) * dpos(1) * (1 - dpos(2)) * pdata;
    fiter.offset(1,1,0) += dpos(0) * dpos(1) * (1 - dpos(2)) * pdata;
    fiter.offset(0,0,1) += (1 - dpos(0)) * (1 - dpos(1)) * dpos(2) * pdata;
    fiter.offset(1,0,1) += dpos(0) * (1 - dpos(1)) * dpos(2) * pdata;
    fiter.offset(0,1,1) += (1 - dpos(0)) * dpos(1) * dpos(2) * pdata;
    fiter.offset(1,1,1) += dpos(0) * dpos(1) * dpos(2) * pdata;
    return;
  }

  // gather particle data from Field using particle position and mesh
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,3U,M,C>& f,
	      const Vektor<PT,3U>& ppos, const M& mesh) {
    CenteringTag<C> ctag;
    Vektor<PT,3U> gpos, dpos, delta;
    NDIndex<3U> ngp;
    CompressedBrickIterator<FT,3U> fiter;
    int lgpoff[3U];
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<3U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (1 - dpos(1)) * (1 - dpos(2)) * (*fiter)
          + dpos(0)       * (1 - dpos(1)) * (1 - dpos(2)) * fiter.offset(1,0,0)
          + (1 - dpos(0)) * dpos(1)       * (1 - dpos(2)) * fiter.offset(0,1,0)
          + dpos(0)       * dpos(1)       * (1 - dpos(2)) * fiter.offset(1,1,0)
          + (1 - dpos(0)) * (1 - dpos(1)) * dpos(2)       * fiter.offset(0,0,1)
          + dpos(0)       * (1 - dpos(1)) * dpos(2)       * fiter.offset(1,0,1)
          + (1 - dpos(0)) * dpos(1)       * dpos(2)       * fiter.offset(0,1,1)
          + dpos(0)       * dpos(1)       * dpos(2)       * fiter.offset(1,1,1);
    return;
  }

  // gather particle data from Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,3U,M,C>& f,
	      const Vektor<PT,3U>& ppos, const M& mesh,
              NDIndex<3U>& ngp, int lgpoff[3U], Vektor<PT,3U>& dpos) {
    CenteringTag<C> ctag;
    Vektor<PT,3U> gpos, delta;
    CompressedBrickIterator<FT,3U> fiter;
    unsigned int d;
    // find nearest grid point for particle position, store in NDIndex obj
    ngp = FindNGP(mesh, ppos, ctag);
    // get position of ngp
    FindPos(gpos, mesh, ngp, ctag);
    // get mesh spacings
    FindDelta(delta, mesh, ngp, ctag);
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // Now we find the offset from ngp to next-lowest grip point (lgp)
    for (d=0; d<3U; ++d) {
      if (gpos(d) > ppos(d)) {
        lgpoff[d] = -1;                // save the offset
        gpos(d) = gpos(d) - delta(d);  // adjust gpos to lgp position
      }
      else {
        lgpoff[d] = 0;                 // save the offset
      }
    }
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // get distance from ppos to gpos
    dpos = ppos - gpos;
    // normalize dpos by mesh spacing
    dpos /= delta;
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (1 - dpos(1)) * (1 - dpos(2)) * (*fiter)
          + dpos(0)       * (1 - dpos(1)) * (1 - dpos(2)) * fiter.offset(1,0,0)
          + (1 - dpos(0)) * dpos(1)       * (1 - dpos(2)) * fiter.offset(0,1,0)
          + dpos(0)       * dpos(1)       * (1 - dpos(2)) * fiter.offset(1,1,0)
          + (1 - dpos(0)) * (1 - dpos(1)) * dpos(2)       * fiter.offset(0,0,1)
          + dpos(0)       * (1 - dpos(1)) * dpos(2)       * fiter.offset(1,0,1)
          + (1 - dpos(0)) * dpos(1)       * dpos(2)       * fiter.offset(0,1,1)
          + dpos(0)       * dpos(1)       * dpos(2)       * fiter.offset(1,1,1);
    return;
  }

  // gather particle data from Field using cached mesh information
  template <class FT, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,3U,M,C>& f,
	      const NDIndex<3U>& ngp, const int lgpoff[3U],
              const Vektor<PT,3U>& dpos) {
    CompressedBrickIterator<FT,3U> fiter;
    // Try to find ngp in local fields and get iterator
    fiter = getFieldIter(f,ngp);
    // adjust position of Field iterator to lgp position
    fiter.moveBy(lgpoff);
    // accumulate into particle attrib
    pdata = (1 - dpos(0)) * (1 - dpos(1)) * (1 - dpos(2)) * (*fiter)
          + dpos(0)       * (1 - dpos(1)) * (1 - dpos(2)) * fiter.offset(1,0,0)
          + (1 - dpos(0)) * dpos(1)       * (1 - dpos(2)) * fiter.offset(0,1,0)
          + dpos(0)       * dpos(1)       * (1 - dpos(2)) * fiter.offset(1,1,0)
          + (1 - dpos(0)) * (1 - dpos(1)) * dpos(2)       * fiter.offset(0,0,1)
          + dpos(0)       * (1 - dpos(1)) * dpos(2)       * fiter.offset(1,0,1)
          + (1 - dpos(0)) * dpos(1)       * dpos(2)       * fiter.offset(0,1,1)
          + dpos(0)       * dpos(1)       * dpos(2)       * fiter.offset(1,1,1);
    return;
  }

};


// IntCIC class -- what the user sees
class IntCIC {

public:
  // constructor/destructor
  IntCIC() {}
  ~IntCIC() {}

  // gather/scatter functions

  // scatter particle data into Field using particle position and mesh
  template <class FT, unsigned Dim, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,Dim,M,C>& f,
	       const Vektor<PT,Dim>& ppos, const M& mesh) {
    IntCICImpl<Dim>::scatter(pdata,f,ppos,mesh);
  }

  // scatter particle data into Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, unsigned Dim, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,Dim,M,C>& f,
	       const Vektor<PT,Dim>& ppos, const M& mesh,
               CacheDataCIC<PT,Dim>& cache) {
    IntCICImpl<Dim>::scatter(pdata,f,ppos,mesh,cache.Index_m,
                             cache.Offset_m,cache.Delta_m);
  }

  // scatter particle data into Field using cached mesh information
  template <class FT, unsigned Dim, class M, class C, class PT>
  static
  void scatter(const FT& pdata, Field<FT,Dim,M,C>& f,
	       const CacheDataCIC<PT,Dim>& cache) {
    IntCICImpl<Dim>::scatter(pdata,f,cache.Index_m,
                             cache.Offset_m,cache.Delta_m);
  }

  // gather particle data from Field using particle position and mesh
  template <class FT, unsigned Dim, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,Dim,M,C>& f,
	      const Vektor<PT,Dim>& ppos, const M& mesh) {
    IntCICImpl<Dim>::gather(pdata,f,ppos,mesh);
  }

  // gather particle data from Field using particle position and mesh
  // and cache mesh information for reuse
  template <class FT, unsigned Dim, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,Dim,M,C>& f,
	      const Vektor<PT,Dim>& ppos, const M& mesh,
              CacheDataCIC<PT,Dim>& cache) {
    IntCICImpl<Dim>::gather(pdata,f,ppos,mesh,cache.Index_m,
                            cache.Offset_m,cache.Delta_m);
  }

  // gather particle data from Field using cached mesh information
  template <class FT, unsigned Dim, class M, class C, class PT>
  static
  void gather(FT& pdata, const Field<FT,Dim,M,C>& f,
	      const CacheDataCIC<PT,Dim>& cache) {
    IntCICImpl<Dim>::gather(pdata,f,cache.Index_m,
                            cache.Offset_m,cache.Delta_m);
  }

};

#endif // INT_CIC_H

/***************************************************************************
 * $RCSfile: IntCIC.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:28 $
 * IPPL_VERSION_ID: $Id: IntCIC.h,v 1.1.1.1 2003/01/23 07:40:28 adelmann Exp $ 
 ***************************************************************************/

