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
#include "Utility/Pstring.h"
#include "Profile/Profiler.h"

//
// include the connection-method-specific class headers
//

// file connection method
#include "DataSource/FileFieldDataSource.h"
#include "DataSource/FilePtclBaseDataSource.h"
#include "DataSource/FilePtclAttribDataSource.h"

// aclvis connection method
#ifdef IPPL_ACLVIS
#include "DataSource/ACLVISFieldDataSource.h"
#include "DataSource/ACLVISPtclBaseDataSource.h"
#include "DataSource/ACLVISPtclAttribDataSource.h"
#endif

// paws connection method
#ifdef IPPL_PAWS
#include "IpplPaws/PawsFieldDataSource.h"
#include "IpplPaws/PawsPtclAttribDataSource.h"
#include "IpplPaws/PawsScalarDataSource.h"
#include "IpplPaws/PawsStringDataSource.h"
#endif


////////////////////////////////////////////////////////////////////////////
// a version of make_DataSourceObject for Field's.
// arguments: name, connection type, transfer metohd, Field
template<class T, unsigned Dim, class M, class C>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
		      Field<T,Dim,M,C>& F) {
  TAU_TYPE_STRING(taustr, "DataSourceObject * (char *, DataConnect *, " +
		  CT(F) + " )");
  TAU_PROFILE("make_DataSourceObject()", taustr, TAU_VIZ);

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  // find what method it is, and make the appropriate DataSourceObject
  DataSourceObject *dso = 0;
  if (method == "aclvis") {
    // create a DataSourceObject for this Field which will connect to
    // the ACL visualization code
#ifdef IPPL_ACLVIS
    dso = new ACLVISFieldDataSource<T,Dim,M,C>(nm, dc, t, F);
#endif

  } else if (method == "paws") {
    // create a DataSourceObject for this Field which will connect to
    // another PAWS application
#ifdef IPPL_PAWS
    dso = new PawsFieldDataSource<T,Dim,M,C>(nm, dc, t, F);
#endif

  } else if (method == "file") {
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
  TAU_TYPE_STRING(taustr, "DataSourceObject * (char *, DataConnect *, int, "
		  + CT(P) + " )");
  TAU_PROFILE("make_DataSourceObject()", taustr, TAU_VIZ);

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "aclvis") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // the ACL visualization code
#ifdef IPPL_ACLVIS
    dso = new ACLVISParticleAttribDataSource<T>(nm, dc, t, P);
#endif

  } else if (method == "paws") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // another PAWS application
#ifdef IPPL_PAWS
    dso = new PawsParticleAttribDataSource<T>(nm, dc, t, P);
#endif

  } else if (method == "file") {
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
// a version of make_DataSourceObject for ParticleBase's
template<class PLayout>
DataSourceObject *
make_DataSourceObject(const char *nm, DataConnect *dc, int t,
		      ParticleBase<PLayout>& P) {
  TAU_TYPE_STRING(taustr, "DataSourceObject * (char *, DataConnect *, int, "
		  + CT(P) + " )");
  TAU_PROFILE("make_DataSourceObject()", taustr, TAU_VIZ);

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "aclvis") {
    // create a DataSourceObject for this ParticleBase which will connect to
    // the ACL visualization code
#ifdef IPPL_ACLVIS
    dso = new ACLVISParticleBaseDataSource<PLayout>(nm, dc, t, P);
#endif

  } else if (method == "paws") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // another PAWS application

  } else if (method == "file") {
    // create a DataSourceObject for this FILE which will connect to
    // a file
    dso = new FileParticleBaseDataSource<PLayout>(nm, dc, t, P);
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
  TAU_TYPE_STRING(taustr, "DataSourceObject * (char *, DataConnect *, int, "
                  + CT(S) + " )");
  TAU_PROFILE("make_DataSourceObject()", taustr, TAU_VIZ);

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "aclvis") {
    // create a DataSourceObject for this ParticleBase which will connect to
    // the ACL visualization code

  } else if (method == "paws") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // another PAWS application
#ifdef IPPL_PAWS
    dso = new PawsScalarDataSource<T>(nm, dc, t, S);
#endif
  } else if (method == "file") {
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
  TAU_TYPE_STRING(taustr, "DataSourceObject * (char *, DataConnect *, int, "
                  + CT(S) + " )");
  TAU_PROFILE("make_DataSourceObject()", taustr, TAU_VIZ);

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "aclvis") {
    // create a DataSourceObject for this ParticleBase which will connect to
    // the ACL visualization code

  } else if (method == "paws") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // another PAWS application
#ifdef IPPL_PAWS
    dso = new PawsStringDataSource<T>(nm, dc, t, S);
#endif
  } else if (method == "file") {
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
