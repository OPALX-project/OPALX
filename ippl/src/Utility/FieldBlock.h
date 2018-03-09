// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef FIELD_BLOCK_H
#define FIELD_BLOCK_H

// class FieldBlock
// 
// The FieldBlock object is used to store the data from several Field
// variables in a single netcdf file. A record dimension is available to
// allow storage of field histories in the same netCDF file. Fields read
// and write to a netcdf file attached to the FieldBlock object.  The write
// and read member functions take a single Field and require a Field
// Variable ID and record number to access data from the netCDF file. The
// record number is defaulted to zero for the case where there is only a
// single record.
//
// J.V.W. Reynders - ACL/LANL July 28, 1996

#define MAX_FNAME_SIZE 30


// forward declarations
template<class T, unsigned Dim> class LField;
template<class T, unsigned Dim, class Mesh, class Centering> class Field;
template<unsigned Dim, class T> class UniformCartesian;
template<unsigned Dim> class FieldLayout;

//----------------------------------------------------------------------
template<class T, unsigned Dim, 
         class Mesh=UniformCartesian<Dim,double>, 
         class Centering=typename Mesh::DefaultCentering>
class FieldBlock {

public:

  // make a FieldBlock for writing
  FieldBlock(char* fname, FieldLayout<Dim>& layout, unsigned numFields);

  // make a FieldBlock for reading
  FieldBlock(char* fname, FieldLayout<Dim>& layout);

  ~FieldBlock(void) { }

  int get_NumRecords(void) { return NumRecords; };
  int get_NumFields(void) { return NumFields; };
  void write(Field<T,Dim,Mesh,Centering>&f, 
	    unsigned varID, unsigned record = 0);
  void read (Field<T,Dim,Mesh,Centering>&f, 
	     unsigned varID, unsigned record = 0);

private:

  char FName[MAX_FNAME_SIZE];
  FieldLayout<Dim>& Layout;
  unsigned NumFields;
  unsigned NumRecords;

  // don't allow copy or assign
  FieldBlock(const FieldBlock&) { };
  FieldBlock& operator=(const FieldBlock&) { return *this; }

};
//----------------------------------------------------------------------

#include "Utility/FieldBlock.hpp"

#endif // FIELD_BLOCK_H

/***************************************************************************
 * $RCSfile: FieldBlock.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: FieldBlock.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
