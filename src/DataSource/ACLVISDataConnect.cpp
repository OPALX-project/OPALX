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
#include "DataSource/ACLVISDataConnect.h"
#include "Utility/IpplInfo.h"

// ACLVIS library headers
#ifdef IPPL_LUX
# include "interface/script.h"
#else
# include "aclvis/aclvis.h"
#endif


// static data members of this class
bool ACLVISDataConnect::ACLVISInitialized = false;


/////////////////////////////////////////////////////////////////////////
// constructor: initialize ACLVIS window for display
ACLVISDataConnect::ACLVISDataConnect(const char *nm, int n)
  : DataConnect(nm, "aclvis", DataSource::OUTPUT, n)
{
  TAU_PROFILE("ACLVISDataConnect::ACLVISDataConnect()",
	      "void (char *)", TAU_VIZ);

  // Inform dbgmsg("ACLVISDataConnect", INFORM_ALL_NODES);
  // dbgmsg << "Setting up new ACLVISDataConnect '" << nm << "'..." << endl;

  if (onConnectNode()) {
    // do one-time ACLVIS initialization
    if ( ! ACLVISInitialized ) {
      // dbgmsg << "Doing aclvis initialization ..." << endl;
#ifdef IPPL_LUX
      Connection = lux_init("IPPL-R1");
#else
      aclvis_init(Ippl::getArgc(), Ippl::getArgv());
#endif
      ACLVISInitialized = true;
    }

    // create ACLVIS object, to initialize a new window
    // if successful, set proper flags
    // dbgmsg << "Creating new VizTool ..." << endl;
#ifndef IPPL_LUX
    Connection = new VizTool();
#endif

  } else {
    // not much to do; we just hope it worked on node 0
    // dbgmsg << "Simple initialization on node " << Ippl::myNode() << endl;
    Connection = 0;
  }
}


/////////////////////////////////////////////////////////////////////////
// destructor: disconnect all data objects, and close window
ACLVISDataConnect::~ACLVISDataConnect() {
  TAU_PROFILE("ACLVISDataConnect::~ACLVISDataConnect()", "void ()", TAU_VIZ);

  // disconnect any remaining data objects
  disconnectConnections();

  // delete the ACLVIS object
  if (Connection != 0)
    delete Connection;
}


/////////////////////////////////////////////////////////////////////////
// are we currently connected to a receiver?  The base-class default
// behavior for this is to indicate that we're not connected.  Here, we
// just say that we are connected.
bool ACLVISDataConnect::connected() const {
  return true;
}


/////////////////////////////////////////////////////////////////////////
// allow all connections to perform an interactive action
// For ACLVIS, since interact passes control on to the graphics window,
// we only want to call interact for the first connected item, instead
// of for all of them, since that one interact will affect all the connected
// (well, displayed) elements.
void ACLVISDataConnect::interact(const char *str, DataConnect *dc) {
  TAU_PROFILE("ACLVISDataConnect::interact()", "void (const char *)", TAU_VIZ);

  // just do interact for first connection
  if (size() > 0)
    (*(begin()))->interact(str, dc);
}


/***************************************************************************
 * $RCSfile: ACLVISDataConnect.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISDataConnect.cpp,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
