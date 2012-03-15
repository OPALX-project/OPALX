// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef SCALAR_DATA_SOURCE_H
#define SCALAR_DATA_SOURCE_H

/***********************************************************************
 * 
 * ScalarDataSource is a DataSource for scalars.
 * It makes a scalar available for another program or data processing 
 * API (e.g., visualization) through the DataSource API. 
 *
***********************************************************************/

// include files
#include "DataSource/DataSource.h"

// forward declarations
class DataSourceObject;
class DataConnect;

// A DataSource class for handling scalars
template<class T>
class ScalarDataSource : public DataSource {

public:
  // constructor
  ScalarDataSource(T& S) : MyScalar(S) {}; 

  // destructor
  virtual ~ScalarDataSource() { }

  // Return ptr to Scalar
  T& scalarRef() { return MyScalar; }

protected:
  // a virtual function which is called by this base class to get a
  // specific instance of DataSourceObject based on the type of data
  // and the connection method (the argument to the call).
  DataSourceObject *createDataSourceObject(const char *,
						   DataConnect *,
						   int);

private:
  // The scalar
  T& MyScalar;
};

#include "DataSource/ScalarDataSource.cpp"

#endif // SCALAR_DATA_SOURCE_H

/***************************************************************************
 * $RCSfile: ScalarDataSource.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: ScalarDataSource.h,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/
