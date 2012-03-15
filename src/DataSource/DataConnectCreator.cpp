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
#include "DataSource/DataConnectCreator.h"
#include "DataSource/DataConnect.h"
#include "DataSource/FileDataConnect.h"
#include "Utility/IpplInfo.h"
#include "Profile/Profiler.h"

// include ACLVIS code if available
#ifdef IPPL_ACLVIS
#include "DataSource/ACLVISDataConnect.h"
#endif

// include PAWS code if available
#ifdef IPPL_PAWS
#include "IpplPaws/PawsDataConnect.h"
#endif

#include <cstring>


// static data for this file
static const int CONNECTMETHODS = 4;   // includes "no connection" method
static const char *ConnectMethodList  = "aclvis, paws, file, or none";
static const char *ConnectMethodNames[CONNECTMETHODS] =
  { "aclvis", "paws", "file", "none" };
static bool  ConnectMethodSupported[CONNECTMETHODS] = {
#ifdef IPPL_ACLVIS
  true,
#else
  false,
#endif
#ifdef IPPL_PAWS
  true,
#else
  false,
#endif
  true,
  true
};


// a global instance of DataConnectCreator ... when this is destroyed
// at the end of the program, it will clean up the default connection
// object, if necessary
static DataConnectCreator GlobalDataConnectCreatorInstance;



/////////////////////////////////////////////////////////////////////////
// static member data for DataConnectCreator
int          DataConnectCreator::ConnectNodes      = 1;
int          DataConnectCreator::InstanceCount     = 0;
int          DataConnectCreator::DefaultMethod     = CONNECTMETHODS - 1;
DataConnect *DataConnectCreator::DefaultConnection = 0;



/////////////////////////////////////////////////////////////////////////
// constructor: increment the instance count
DataConnectCreator::DataConnectCreator() {
  TAU_PROFILE("DataConnectCreator::DataConnectCreator()", "void ()", TAU_VIZ);
  InstanceCount++;
}


/////////////////////////////////////////////////////////////////////////
// destructor: deccrement the instance count, and if it goes to zero,
// delete any static object if necessary
DataConnectCreator::~DataConnectCreator() {
  TAU_PROFILE("DataConnectCreator::~DataConnectCreator()", "void ()", TAU_VIZ);
  if (--InstanceCount == 0)
    if (DefaultConnection != 0)
      delete DefaultConnection;
}


/////////////////////////////////////////////////////////////////////////
// return the name of the Nth method
int DataConnectCreator::getNumMethods() {
  TAU_PROFILE("DataConnectCreator::getNumMethods()", "int ()", TAU_VIZ);
  return CONNECTMETHODS;
}


/////////////////////////////////////////////////////////////////////////
// return the name of the Nth method
const char *DataConnectCreator::getMethodName(int n) {
  TAU_PROFILE("DataConnectCreator::getMethodName()", "char * (int )", TAU_VIZ);
  if (n >= 0 && n < CONNECTMETHODS)
    return ConnectMethodNames[n];
  else
    return 0;
}


/////////////////////////////////////////////////////////////////////////
// return a list of all the methods, as a single string
const char *DataConnectCreator::getAllMethodNames() {
  TAU_PROFILE("DataConnectCreator::getAllMethodNames()", "char * ()", TAU_VIZ);
  return ConnectMethodList;
}


/////////////////////////////////////////////////////////////////////////
// check if the given connection method is supported
bool DataConnectCreator::supported(int cm) {
  TAU_PROFILE("DataConnectCreator::supported()", "bool (int )", TAU_VIZ);
  return (known(cm) ? ConnectMethodSupported[cm] : false);
}


/////////////////////////////////////////////////////////////////////////
// check if the given connection method is supported
bool DataConnectCreator::supported(const char *nm) {
  TAU_PROFILE("DataConnectCreator::supported()", "bool (char *)", TAU_VIZ);
  return supported(libindex(nm));
}


/////////////////////////////////////////////////////////////////////////
// check if the given connection method is known at all
bool DataConnectCreator::known(int cm) {
  TAU_PROFILE("DataConnectCreator::known()", "bool (int )", TAU_VIZ);
  return (cm >= 0 && cm < CONNECTMETHODS);
}


/////////////////////////////////////////////////////////////////////////
// check if the given connection method is known at all
bool DataConnectCreator::known(const char *nm) {
  TAU_PROFILE("DataConnectCreator::known()", "bool (char *)", TAU_VIZ);
  return known(libindex(nm));
}


/////////////////////////////////////////////////////////////////////////
// create a new connection.  Arguments = type, name, direction, nodes
// If n <= 0, use the "default" number of nodes set earlier
DataConnect *DataConnectCreator::create(int cm, const char *nm, int n) {
  TAU_PROFILE("DataConnectCreator::create()", "DataConnect * (int, char * )",
	      TAU_VIZ);
  
  // initially, we have a null pointer for the connection.  If everything
  // checks out, we'll have a non-null pointer at the end.  If we still
  // have a null pointer at the end of this routine, something went wrong
  // and we return NULL to indicate an error.
  DataConnect *dataconn = 0;

  // figure out how many nodes the connection should use, if it cares
  // at all about that.  For example, ACLVIS may use multipple nodes
  // in the visualization, or maybe just one.
  int nodes = n;
  if (n <= 0)
    nodes = getDefaultNodes();

  // based on the connection method, create a new DataConnect
  if (cm == 0) {
    // use the ACLVIS library for the connection, which will result in the
    // data objects being displayed in a visualization window
#ifdef IPPL_ACLVIS
    dataconn = new ACLVISDataConnect(nm, nodes);
#endif
  } else if (cm == 1) {
    // use the PAWS library for the connection, which is a general-purpose
    // method for transferring data to/from another application
#ifdef IPPL_PAWS
    dataconn = new PawsDataConnect(nm, nodes);
#endif
  } else if (cm == 2) {
    // transfer the data to/from a file using some form of parallel I/O
    dataconn = new FileDataConnect(nm, nodes);
  } else if (cm == 3) {
    // just make a dummy connect object, which does nothing
    dataconn = new DataConnect(nm, getMethodName(cm), DataSource::OUTPUT,
			       nodes);
  }

  // make sure we have something
  if (dataconn == 0) {
    ERRORMSG("DataConnectCreator: unknown connection method." << endl);
  }

  return dataconn;
}


/////////////////////////////////////////////////////////////////////////
// a final method for creating objects; this one provides a default name,
// and if the default connection object has already been created, this just
// returns that one.
DataConnect *DataConnectCreator::create() {
  TAU_PROFILE("DataConnectCreator::create()", "DataConnect * ()", TAU_VIZ);

  if (DefaultConnection == 0)
    DefaultConnection = create(getMethodName(DefaultMethod));
  return DefaultConnection;
}


/////////////////////////////////////////////////////////////////////////
// change the default connection method.  Return success.
bool DataConnectCreator::setDefaultMethod(int cm) {
  TAU_PROFILE("DataConnectCreator::setDefaultMethod()", "bool (int)", TAU_VIZ);
  if (supported(cm)) {
    DefaultMethod = cm;
    return true;
  }
  return false;
}


/////////////////////////////////////////////////////////////////////////
// return the index of the given named method, or (-1) if not found
int DataConnectCreator::libindex(const char *nm) {
  TAU_PROFILE("DataConnectCreator::libindex()", "int (char * )", TAU_VIZ);
  for (int i=0; i < CONNECTMETHODS; ++i) {
    if (strcmp(nm, getMethodName(i)) == 0)
      return i;
  }

  // if here, it was not found
  return (-1);
}


/////////////////////////////////////////////////////////////////////////
// change the default number of nodes to use for the connection
void DataConnectCreator::setDefaultNodes(int n) {
  ConnectNodes = n;
}


/////////////////////////////////////////////////////////////////////////
// return the default number of nodes to use in a connection
int DataConnectCreator::getDefaultNodes() {
  return (ConnectNodes >= 0 && ConnectNodes <= Ippl::getNodes() ?
	  ConnectNodes :
	  Ippl::getNodes());
}


/***************************************************************************
 * $RCSfile: DataConnectCreator.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: DataConnectCreator.cpp,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
