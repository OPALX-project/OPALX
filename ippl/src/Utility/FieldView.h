// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef FIELD_VIEW_H
#define FIELD_VIEW_H

//----------------------------------------------------------------------

// class FieldView
// 
// The FieldView class produces object which provide a viewing port for the
// contents of a 2D Field through the GDL interface. A FielView object must
// be constructed with a 2D Field object. Thie FieldView instantiation must
// have the same template parameters as the Field being used to construct
// the FieldView object. The FieldView runs in parallel by coalescing all
// the Field data onto a parent node and redering from that node. The
// parent node may be chosen with a constructor argument (the default is
// process zero). For small Fields, the data may be scaled up by repeating
// values across pixels in the GDL window, the number of pixels to repeat
// for a given data point in the field can be set with the scaleX and
// scaleY arguments in the constructor. The minSizeX and minSizeY arguments
// to the constructor specifies a minimum number of pixels in the x and y
// direction. If the data set is small, the scale factors will be increased
// to the smallest integer satisfying the requirements specifiec by
// minSizeX and minSizeY.
// 
// FieldView objects can also be constructed from 3D Fields given
// an axis which is considered perpendicular to the desired 2D plane 
// of view. For 3D, the view member function requires an unsigned to
// specify which slice along the perpendicular axis to slice.
//
// J.V.W. Reynders - ACL/LANL July 15, 1996

// forward declarations
template<class T, unsigned D> class LField;
template<class T, unsigned D, class M, class C> class Field;
template<unsigned D, class T> class UniformCartesian;


//----------------------------------------------------------------------
template<class T, unsigned Dim, 
         class Mesh=UniformCartesian<Dim,double>, 
         class Centering=typename Mesh::DefaultCentering>
class FieldView {

public:

  // attach a 2D Field to a FieldView
  FieldView(Field<T,Dim,Mesh,Centering>& f, 
	    unsigned scaleX = 4, unsigned scaleY = 4,
	    unsigned minSizeX = 200, 
	    unsigned minSizeY = 200, 
	    unsigned parent = 0);

  // attach a 3D Field to a FieldView 
  FieldView(unsigned sliceDim, Field<T,Dim,Mesh,Centering>& f, 
	    unsigned scaleX = 4, unsigned scaleY = 4,
	    unsigned minSizeX = 200, 
	    unsigned minSizeY = 200, 
	    unsigned parent = 0);

  ~FieldView();
  void void_view(int& r);
  void void_view(unsigned, int& r);
  int view() { int r = 0; void_view(r); return r; }
  int view(unsigned s) { int r = 0; void_view(s,r); return r; }

private:
  T* Data;
  int* MapX;
  int* MapY;
  int dummy;
  Field<T,Dim,Mesh,Centering>& MyField;
  LField<T,2U>* MyLField; 

  unsigned SliceDim;   // for 3D - select the axis
  unsigned ScaleX, ScaleY;
  unsigned MinSizeX, MinSizeY;
  unsigned Parent;
  unsigned SizeX, SizeY;

  // draw all the data together onto the Parent process for viewing
  void update_2D_data();
  void update_3D_data(unsigned slice);

  // form map to fit the data into the viewing port
  void init_map();

  // form map to fit the data into the viewing port
  void void_apply_map(int& r);
  int apply_map() { int r; void_apply_map(r); return r; }

};
//----------------------------------------------------------------------

#include "Utility/FieldView.hpp"

#endif // FIELD_VIEW_H

/***************************************************************************
 * $RCSfile: FieldView.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: FieldView.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
