// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef STRING_DATA_SOURCE_H
#define STRING_DATA_SOURCE_H

/***********************************************************************
 * 
 * StringDataSource is a DataSource for strings.
 * It makes a string available for another program or data processing 
 * API (e.g., visualization) through the DataSource API. 
 *
***********************************************************************/

// include files
#include "DataSource/DataSource.h"

// forward declarations
class DataSourceObject;
class DataConnect;


// A DataSource class for handling strings
// A string can be a char* (T=char) or std::string (T=string)
template <class T>
class StringDataSource : public DataSource {

public:
  // constructor
  StringDataSource(T* S, int mlen) : MyString(S), StringLen(mlen) {}; 

  // destructor
  virtual ~StringDataSource() { }

  // Return ptr to Scalar
  T *stringPtr() { return MyString; }
  const T *stringPtr() const { return MyString; }
  
  // Return max length of string
  int stringLen() const { return StringLen; }

protected:
  // a virtual function which is called by this base class to get a
  // specific instance of DataSourceObject based on the type of data
  // and the connection method (the argument to the call).
  virtual DataSourceObject *createDataSourceObject(const char *,
						   DataConnect *,
						   int);

private:
  // The string
  T* MyString;
  int StringLen;
};

#include "DataSource/StringDataSource.cpp"

#endif // STRING_DATA_SOURCE_H

/***************************************************************************
 * $RCSfile: StringDataSource.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: StringDataSource.h,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/
