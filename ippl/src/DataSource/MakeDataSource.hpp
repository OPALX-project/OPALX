// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "DataSource/MakeDataSource.h"
#include "DataSource/DataSourceObject.h"
#include "DataSource/DataConnect.h"

//
// include the connection-method-specific class headers
//

// file connection method
#include "DataSource/FileFieldDataSource.h"
#include "DataSource/FilePtclBaseDataSource.h"
#include "DataSource/FilePtclAttribDataSource.h"

////////////////////////////////////////////////////////////////////////////
// a version of make_DataSourceObject for Field's.
// arguments: name, connection type, transfer metohd, Field
template<class T, unsigned Dim, class M, class C>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
		      Field<T,Dim,M,C>& F) {

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  // find what method it is, and make the appropriate DataSourceObject
  DataSourceObject *dso = 0;
  if (method == "file") {
    // create a DataSourceObject for this Field which will connect to a file
    dso = new FileFieldDataSource<T,Dim,M,C>(nm, dc, t, F);
  }

  // make a default connection is nothing has been found
  if (dso == 0)
    dso = new DataSourceObject;

  return dso;
}


////////////////////////////////////////////////////////////////////////////
// a version of make_DataSourceObject for ParticleAttrib's
template<class T>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
		      ParticleAttrib<T>& P) {

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "file") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // a file
    dso = new FileParticleAttribDataSource<T>(nm, dc, t, P);
  }

  // make a default connection is nothing has been found
  if (dso == 0)
    dso = new DataSourceObject;

  return dso;
}


////////////////////////////////////////////////////////////////////////////
// a version of make_DataSourceObject for IpplParticleBase's
template<class PLayout>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
		      IpplParticleBase<PLayout>& P) {

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "file") {
    // create a DataSourceObject for this FILE which will connect to
    // a file
    dso = new FileIpplParticleBaseDataSource<PLayout>(nm, dc, t, P);
  }

  // make a default connection is nothing has been found
  if (dso == 0)
    dso = new DataSourceObject;

  return dso;
}

////////////////////////////////////////////////////////////////////////////
// a version of make_DataSourceObject for ScalarDataSource's
template<class T>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
                      ScalarDataSource<T>& S) {

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "file") {
    // create a DataSourceObject for this FILE which will connect to
    // a file
  }

  // make a default connection is nothing has been found
  if (dso == 0)
    dso = new DataSourceObject;

  return dso;
}


////////////////////////////////////////////////////////////////////////////
// a version of make_DataSourceObject for StringDataSource's
template<class T>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
                      StringDataSource<T>& S) {
  
  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "file") {
    // create a DataSourceObject for this FILE which will connect to
    // a file
  }

  // make a default connection is nothing has been found
  if (dso == 0)
    dso = new DataSourceObject;

  return dso;
}

/***************************************************************************
 * $RCSfile: MakeDataSource.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: MakeDataSource.cpp,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/
