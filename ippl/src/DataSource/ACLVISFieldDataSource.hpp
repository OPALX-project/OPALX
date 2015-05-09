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
#include "DataSource/ACLVISFieldDataSource.h"
#include "DataSource/ACLVISOperations.h"
#include "DataSource/ACLVISDataConnect.h"
#include "Field/Field.h"
#include "Message/Communicate.h"
#include "Utility/IpplInfo.h"


////////////////////////////////////////////////////////////////////////////
// constructor: the name, the connection, the transfer method,
// the field to connect, and the parent node.
template<class T, unsigned Dim, class M, class C>
ACLVISFieldDataSource<T,Dim,M,C>::ACLVISFieldDataSource(const char *nm,
							DataConnect *dc,
							int tm,
							Field<T,Dim,M,C>& F)
  : FieldDataSource<T,Dim,M,C>(nm, dc, tm, F) {

  // Inform dbgmsg("ACLVISFieldDataSource", INFORM_ALL_NODES);
  // dbgmsg << "Setting up new ACLVISFieldDataSource '" << nm << "'...";
  // dbgmsg << endl;

  // do general initialization
  LocalData = 0;
  object_name = (char *)NULL;
  std::string filestring = "aclvis";
  if (std::string(dc->DSID()) != filestring) {
    ERRORMSG("Illegal DataConnect object for ACLVIS Data Object." << endl);
    ACLVISConnection = 0;
    Connection = 0;
  } else if (tm != DataSource::OUTPUT) {
    ERRORMSG("ACLVIS data connections may only be of type OUTPUT." << endl);
    ACLVISConnection = 0;
    Connection = 0;
  } else {    
    ACLVISConnection = (ACLVISDataConnect *)dc;

    // only do connection on parent node
    // dbgmsg << "Doing viz field initialization on node " << Ippl::myNode();
    // dbgmsg << endl;
    if (dc->onConnectNode()) {
      // create a new ACLVIS structure to hold the data
      LocalData = new ReadFieldTool;
      object_name = new char[strlen(nm)+2];
      strcpy(object_name,nm);

      // ... then call ACLVIS API to register this data
      ACLVISConnection->getConnection()->connect((void *)this, (char *)nm,
						 LocalData, FieldDataType);
    }
  }
}


////////////////////////////////////////////////////////////////////////////
// destructor
template<class T, unsigned Dim, class M, class C>
ACLVISFieldDataSource<T,Dim,M,C>::~ACLVISFieldDataSource() {

  // call ACLVIS API to unregister this data.
  if (ACLVISConnection->onConnectNode() && connected()) {
    // only do disconnection on parent node
#ifdef IPPL_LUX
    ACLVISConnection->getConnection()->disconnect(object_name);
#else  
    ACLVISConnection->getConnection()->disconnect((void *)this);
#endif

    // delete ACLVIS storage
    if (LocalData != 0) {
      delete LocalData;
      if(object_name)
        delete [] object_name;
    }
  }
}


////////////////////////////////////////////////////////////////////////////
// a function to set the mesh characteristics of the ReadFieldTool
template<class T, unsigned Dim, class M, class C>
void ACLVISFieldDataSource<T,Dim,M,C>::set_mesh(Field<T,Dim,M,C>& F) {
  
  

  // only change mesh on parent node
  if (ACLVISConnection->onConnectNode() && connected()) {
    // set the dimension of each axis
    int size[3];
    for (unsigned i=0; i < 3; i++)
      size[i] = (i < Dim ? F.getDomain()[i].length() : 1);
    LocalData->SetDimensions(size);

    // get the field mesh information
    unsigned int sdim;
    float spacing[3], origin[3];
    NDIndex<Dim> mo;
    for (sdim=0; sdim < Dim; ++sdim)
      mo[sdim] = Index(1);
    typename M::MeshVektor_t morigin  = F.get_mesh().get_origin();
    typename M::MeshVektor_t mspacing = F.get_mesh().getDeltaVertex(mo);
    for (sdim=0; sdim < 3; ++sdim) {
      origin[sdim]  = (sdim < Dim ? (float)(morigin[sdim])  : 0.0);
      spacing[sdim] = (sdim < Dim ? (float)(mspacing[sdim]) : 0.0);
    }

    // store mesh data in ReadFieldTool
    LocalData->SetAspectRatio(spacing);
    LocalData->SetOrigin(origin);
  }
}


////////////////////////////////////////////////////////////////////////////
// Update the object, that is, make sure the receiver of the data has a
// current and consistent snapshot of the current state.  Return success.
template<class T, unsigned Dim, class M, class C>
bool ACLVISFieldDataSource<T,Dim,M,C>::update() {
  
  

  if (ACLVISConnection->onConnectNode() && connected()) {
    // tell the ACLVIS library what type and how much data to expect
    LocalData->GetVizData()->InitData(MyField.getDomain().size(),
			      ACLVISTraits<ReadFieldTool,T>::getType());
  }

  // update the mesh data
  set_mesh(MyField);

  // gather all the data into the vtk object
  gather_data();

  if (ACLVISConnection->onConnectNode() && connected()) {
    // give the data to the vtk data set
    LocalData->PrepareFinishedData();

    // on parent node, call ACLVIS API to update data
#ifdef IPPL_LUX
    ACLVISConnection->getConnection()->update(object_name);
#else  
    ACLVISConnection->getConnection()->update((void *)this);
#endif
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////
// Indicate to the receiver that we're allowing them time to manipulate the
// data (e.g., for a viz program, to rotate it, change representation, etc.)
// This should only return when the manipulation is done.
template<class T, unsigned Dim, class M, class C>
void ACLVISFieldDataSource<T,Dim,M,C>::interact(const char *str) {
  
  

  // on parent node, hand off control to ACLVIS API.  If a command string
  // is given, instead call the viz API function to execute it as an
  // interactive command.
  if (ACLVISConnection->onConnectNode() && connected()) {
    if (str != 0 && *str != 0) {
#ifndef IPPL_LUX
      ACLVISConnection->getConnection()->InterpretCommand(str);
#endif
    } else {
#ifdef IPPL_LUX
      ACLVISConnection->getConnection()->Interact();
#else
      ACLVISConnection->getConnection()->Interact(1);
#endif
    }
  }
}


////////////////////////////////////////////////////////////////////////////
// copy the data out of the given LField iterator (which is occupying the
// given domain) and into the ACLVIS structure
template<class T, unsigned Dim, class M, class C >
void ACLVISFieldDataSource<T,Dim,M,C>::insert_data(
                   const NDIndex<Dim>& dom,
		   CompressedBrickIterator<T,Dim> rhs) { 

  // if we're not connected, just ignore the data
  if (! connected())
    return;

  // find the values needed to select where in the ACLVIS object to put the
  // data
  int Stride[Dim], Length[Dim], Count[Dim], i;
  int indx = 0, n = 1;
  for (i=0; i < Dim; ++i) {
    Count[i] = 0;
    Stride[i] = n;
    Length[i] = dom[i].length();
    indx += n*(dom[i].first() - (MyField.getDomain())[i].first());
    n *= (MyField.getDomain())[i].length();
  }

  // go through the data, and put it into the proper position
  for (unsigned items = dom.size(); items > 0; --items, ++rhs) {
    // copy the next element
    ACLVISTraits<ReadFieldTool,T>::setPoint(LocalData, indx, *rhs);

    // increment the index
    for (i=0; i < Dim; ++i) {
      indx += Stride[i];
      if (++Count[i] == Length[i]) {
	Count[i] = 0;
	indx -= (Stride[i] * Length[i]);
      } else {
	break;
      }
    }
  }
}


/***************************************************************************
 * $RCSfile: ACLVISFieldDataSource.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISFieldDataSource.cpp,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
