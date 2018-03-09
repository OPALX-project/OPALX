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
#include "DataSource/ScalarDataSource.h"
#include "DataSource/MakeDataSource.h"
#include "DataSource/DataConnect.h"


//////////////////////////////////////////////////////////////////////////
// a virtual function which is called by this base class to get a
// specific instance of DataSourceObject based on the type of data
// and the connection method (the argument to the call).
template<class T>
DataSourceObject *ScalarDataSource<T>::createDataSourceObject(const char *nm, DataConnect *dc, int tm) {
  return make_DataSourceObject(nm, dc, tm, *this);
}




/***************************************************************************
 * $RCSfile: ScalarDataSource.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: ScalarDataSource.cpp,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/
