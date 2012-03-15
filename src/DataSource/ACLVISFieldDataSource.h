// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef ACLVIS_FIELD_DATA_SOURCE_H
#define ACLVIS_FIELD_DATA_SOURCE_H

/***********************************************************************
 * 
 * class ACLVISFieldDataSource
 *
 * ACLVISFieldDataSource is a specific version of FieldDataSource which takes
 * the data for a given Field and formats it properly for connection to the
 * ACL visualization code.  The initial version
 * collects all the data onto node 0, and then formats it for ACLVIS use.
 * Future versions will properly redistribute the data based on
 * the needs of the recipient.
 *
 ***********************************************************************/

// include files
#include "DataSource/FieldDataSource.h"


// forward declarations
class ACLVISDataConnect;
class ReadFieldTool;


template<class T, unsigned Dim, class M, class C>
class ACLVISFieldDataSource : public FieldDataSource<T,Dim,M,C> {

public:
  // constructor: the name, the connection, the transfer method,
  // the field to connect, and the parent node.
  ACLVISFieldDataSource(const char *, DataConnect *, int, Field<T,Dim,M,C>&);

  // destructor
  virtual ~ACLVISFieldDataSource();

  //
  // virtual function interface.
  //

  // Update the object, that is, make sure the receiver of the data has a
  // current and consistent snapshot of the current state.  Return success.
  virtual bool update();

  // Indicate to the receiver that we're allowing them time to manipulate the
  // data (e.g., for a viz program, to rotate it, change representation, etc.)
  // This should only return when the manipulation is done.
  virtual void interact(const char * = 0);

protected:
  // copy the data out of the given LField iterator (which is occupying the
  // given domain) and into the ACLVIS structure
  virtual void insert_data(const NDIndex<Dim>&,CompressedBrickIterator<T,Dim>);

private:
  // the ACLVIS data structure used to hold the data
  ReadFieldTool *LocalData;

  // recast of DataConnect to ACL-specific object
  char *object_name;
  ACLVISDataConnect *ACLVISConnection;

  // a function to set the mesh characteristics of the ReadFieldTool
  void set_mesh(Field<T,Dim,M,C>&);
};

#include "DataSource/ACLVISFieldDataSource.cpp"

#endif // ACLVIS_FIELD_DATA_SOURCE_H

/***************************************************************************
 * $RCSfile: ACLVISFieldDataSource.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISFieldDataSource.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
