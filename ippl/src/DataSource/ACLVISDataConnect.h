// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef ACLVIS_DATA_CONNECT_H
#define ACLVIS_DATA_CONNECT_H

/***********************************************************************
 * 
 * ACLVISDataConnect represents a single connection to an ACLVIS display
 * window.  This window can have multipple data objects displayed within
 * it.  Only 'send' connections can be made to it.  When constructed, it
 * creates a new window, and readies it for new data.
 *
 ***********************************************************************/

// include files
#include "DataSource/DataConnect.h"

// forward declarations
class VizTool;


class ACLVISDataConnect : public DataConnect {

public:
  // constructor: the name for the window, and the number of nodes to use
  ACLVISDataConnect(const char *, int = 0);

  // destructor: shut down the connection
  virtual ~ACLVISDataConnect();

  //
  // methods specific to this type of DataConnect
  //

  // return the ACLVIS connection object
  VizTool *getConnection() const { return Connection; }

  //
  // DataConnect virtual methods
  //

  // are we currently connected to a receiver?
  virtual bool connected() const;

  // allow all connections to perform an interactive action.  An optional
  // command string can be supplied; if it is null, it will be ignored.
  // For ACLVIS, since interact passes control on to the graphics window,
  // we only want to call interact for the first connected item, instead
  // of for all of them, since that one interact will affect all the connected
  // (well, displayed) elements.
  virtual void interact(const char * = 0, DataConnect * = 0);

private:
  // the ACLVIS connection object
  VizTool *Connection;

  // have we initialized the ACLVIS code yet?
  static bool ACLVISInitialized;
};

#endif // ACLVIS_DATA_CONNECT_H

/***************************************************************************
 * $RCSfile: ACLVISDataConnect.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISDataConnect.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
