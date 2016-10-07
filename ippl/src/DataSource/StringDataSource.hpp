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
#include "DataSource/DataSourceObject.h"
#include "DataSource/DataConnect.h"
#include  


#ifdef IPPL_PAWS
#include "IpplPaws/PawsStringDataSource.h"
#endif


//////////////////////////////////////////////////////////////////////////
// a virtual function which is called by this base class to get a
// specific instance of DataSourceObject based on the type of data
// and the connection method (the argument to the call).
template <class T>
DataSourceObject *StringDataSource<T>::createDataSourceObject(const char *nm,
							   DataConnect *dc,
							   int tm) {
  
  

  // get the connection method name, and make a string out of it
  std::string method(dc->DSID());

  DataSourceObject *dso = 0;
  if (method == "aclvis") {
    // create a DataSourceObject for this IpplParticleBase which will connect to
    // the ACL visualization code

  } else if (method == "paws") {
    // create a DataSourceObject for this ParticleAttrib which will connect to
    // another PAWS application
#ifdef IPPL_PAWS
    dso = new PawsStringDataSource<T>(nm, dc, tm, *this);
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
 * $RCSfile: StringDataSource.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: StringDataSource.cpp,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/
