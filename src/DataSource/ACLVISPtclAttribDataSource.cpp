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
#include "DataSource/ACLVISPtclAttribDataSource.h"
#include "DataSource/ACLVISOperations.h"
#include "Particle/ParticleAttrib.h"
#include "Particle/ParticleLayout.h"
#include "Message/Message.h"
#include "Utility/IpplInfo.h"
#include "Utility/Pstring.h"
#include "Profile/Profiler.h"


////////////////////////////////////////////////////////////////////////////
// constructor: the name, the connection, the transfer method, the attrib
template<class T>
ACLVISParticleAttribDataSource<T>::ACLVISParticleAttribDataSource(const
								  char *nm,
      DataConnect *dc, int tm, ParticleAttrib<T>& pa)
  : ParticleAttribDataSource(nm, dc, tm, &pa, &pa), MyParticles(pa) {

  TAU_TYPE_STRING(taustr, "void (char *, DataConnect *, int, " +
		  CT(pa));
  TAU_PROFILE(
    "ACLVISParticleAttribDataSource::ACLVISParticleAttribDataSource()",
    taustr, TAU_VIZ);

  // do general initialization
  LocalData = 0;
  std::string filestring = "aclvis";
  if (std::string(dc->DSID()) != filestring) {
    ERRORMSG("Illegal DataConnect object for ACLVIS Data Object." << endl);
    Connection = 0;
  } else if (tm != DataSource::OUTPUT) {
    ERRORMSG("ACLVIS data connections may only be of type OUTPUT." << endl);
    Connection = 0;
  } else {    
    // check if our ParticleBase container has been connected already ...
    if (PBase != 0) {
      // yes it has ... set up our ACLVIS info.
      // make a Tool object to store the data on the parent node
      // We do not need to register this with the ACLVIS API, as this is
      // done by the ACLVISParticleBaseDataSource
      if (getConnection()->onConnectNode())
	LocalData = new ReadParticleTool;

      // tell the PBase we're connected ... we had to wait until now to do
      // this, because we had to have the ReadParticleTool object available
      PBase->connect_attrib(this);

    } else {
      // no luck; indicate we did not get connected
      Connection = 0;
    }
  }
}


////////////////////////////////////////////////////////////////////////////
// destructor
template<class T>
ACLVISParticleAttribDataSource<T>::~ACLVISParticleAttribDataSource() {

  TAU_TYPE_STRING(taustr, CT(*this) + " void ()");
  TAU_PROFILE("ACLVISParticleAttribDataSource::~ACLVISParticleAttribDataSource()", 
    taustr, TAU_VIZ);

  // disassociate ourselves from the ParticleBase container
  if (PBase != 0)
    PBase->disconnect_attrib(this);

  // delete local storage
  if (LocalData != 0)
    delete LocalData;
}


////////////////////////////////////////////////////////////////////////////
// Update the object, that is, make sure the receiver of the data has a
// current and consistent snapshot of the current state.  Return success.
// In fact, this version of update does nothing.  Only the update for
// particlebase actually does any work.  This version will just silently
// return (since it may be called if updateConnections is called for the
// DataConnect object).
template<class T>
bool ACLVISParticleAttribDataSource<T>::update() {
  TAU_TYPE_STRING(taustr, CT(*this) + " bool ()");
  TAU_PROFILE("ACLVISParticleAttribDataSource::update()", taustr, TAU_VIZ);

  // just return, the user should do an update for the ParticleBase
  return true;
}


////////////////////////////////////////////////////////////////////////////
// Indicate to the receiver that we're allowing them time to manipulate the
// data (e.g., for a viz program, to rotate it, change representation, etc.)
// This should only return when the manipulation is done.
template<class T>
void ACLVISParticleAttribDataSource<T>::interact(const char *str) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void ()");
  TAU_PROFILE("ACLVISParticleAttribDataSource::interact()", taustr, TAU_VIZ);

  // hand off control to ACLVIS API
  if (PBase != 0)
    PBase->interact(str);
}


////////////////////////////////////////////////////////////////////////////
// put the local particle data into a message
template<class T>
void ACLVISParticleAttribDataSource<T>::putMessage(Message *msg) {

  TAU_TYPE_STRING(taustr, CT(*this) + " void (Message *)");
  TAU_PROFILE("ACLVISParticleAttribDataSource::putMessage()", taustr, TAU_VIZ);

  unsigned N = MyParticles.size();

  //  Inform dbgmsg("attrib_putmessage", INFORM_ALL_NODES);
  //  dbgmsg << "Putting ParticleAttrib's in message: num = " << N << endl;

  // put in the data
  if (N > 0) {
    // WARNMSG("ACLVIS putMessage: on node " << Ippl::myNode() << ", N = ");
    // WARNMSG(N << endl);
    MyParticles.putMessage(*msg, N, 0);
    // WARNMSG("Message now: " << *msg << endl);
  }
}


////////////////////////////////////////////////////////////////////////////
// prepare the agency-specific data structures for update; this may
// require reallocation of storage memory, etc.
// Argument = are we at the start (true) or end (false) of the data update;
//            # of particles to prepare for.
template<class T>
void ACLVISParticleAttribDataSource<T>::prepare_data(bool start, unsigned N) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (bool, unsigned)");
  TAU_PROFILE("ACLVISParticleAttribDataSource::prepare_data()", taustr, TAU_VIZ);

  if (start) {
    // reallocate new storage for the attribute and coordinate data
    LocalData->GetVizData()->InitData(N,
	      ACLVISTraits<ReadParticleTool,T>::getType(), 1);
  } else {
    // tell data storage objects that we're done adding data for this step
    LocalData->PrepareFinishedData();
  }
}


////////////////////////////////////////////////////////////////////////////
// copy the data out of the given Message and into the proper vtk structure.
// If Message is 0, just put in the data from our local particles.
// Arguments: num particles, starting index for inserted particles,
// Message with data, total number of particles in entire system,
// iterators for position and ID data.
template<class T>
void ACLVISParticleAttribDataSource<T>::insert_data(unsigned N,unsigned sIndx,
						    Message* msg) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (unsigned, unsigned, Message *)");
  TAU_PROFILE("ACLVISParticleAttribDataSource::insert_data()", taustr, TAU_VIZ);

  // maximum index for the data
  unsigned maxIndx = sIndx + N;

  if (msg != 0) {
    // insert particles from a message into the vtk structure
    T* msgdata = (T *)(msg->remove());
    if (msgdata != 0) {
      T *data = msgdata;
      for ( ; sIndx < maxIndx; ++sIndx, ++data)
	ACLVISTraits<ReadParticleTool,T>::setPoint(LocalData, sIndx, *data);

      // delete unneeded storage
      free(static_cast<void *>(msgdata));
    }

  } else {
    // grab particles from our local storage and put them into the vtk struct
    typename ParticleAttrib<T>::iterator msgdata = MyParticles.begin();
    for ( ; sIndx < maxIndx; ++sIndx, ++msgdata)
      ACLVISTraits<ReadParticleTool,T>::setPoint(LocalData, sIndx, *msgdata);
  }
}


/***************************************************************************
 * $RCSfile: ACLVISPtclAttribDataSource.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISPtclAttribDataSource.cpp,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
