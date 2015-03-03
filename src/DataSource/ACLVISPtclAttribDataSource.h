// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef ACLVIS_PARTICLE_ATTRIB_DATA_SOURCE_H
#define ACLVIS_PARTICLE_ATTRIB_DATA_SOURCE_H

/***********************************************************************
 * 
 * class ACLVISParticleAttribDataSource
 *
 * ACLVISParticleAttribDataSource is a specific version of
 * ParticleAttribDataSource which takes the data for a given ParticleAttrib
 * and formats it properly for connection to the ACL visualization code,
 * which uses the VTK toolkit.  The initial version
 * collects all the data onto node 0, and then formats it for VTK use.
 * Future versions will properly redistribute the data based on
 * the needs of the recipient.
 *
 ***********************************************************************/

// include files
#include "DataSource/PtclAttribDataSource.h"


// forward declarations
template<class T> class ParticleAttrib;
class ReadParticleTool;


template<class T>
class ACLVISParticleAttribDataSource : public ParticleAttribDataSource {

public:
  // constructor: the name, the connection, the transfer method, the attrib
  ACLVISParticleAttribDataSource(const char *, DataConnect *, int,
				 ParticleAttrib<T>&);

  // destructor
  virtual ~ACLVISParticleAttribDataSource();

  //
  // DataSourceObject virtual function interface.
  //

  // Update the object, that is, make sure the receiver of the data has a
  // current and consistent snapshot of the current state.  Return success.
  virtual bool update();

  // Indicate to the receiver that we're allowing them time to manipulate the
  // data (e.g., for a viz program, to rotate it, change representation, etc.)
  // This should only return when the manipulation is done.
  virtual void interact(const char * = 0);

  //
  // ParticleAttribDataSource virtual function interface.
  //

  // retrieve the agency-specific data structure
  virtual void *getConnectStorage() { return (void *)LocalData; }

  // copy attrib data on the local processor into the given Message.
  virtual void putMessage(Message *);

  // prepare the agency-specific data structures for update; this may
  // require reallocation of storage memory, etc.
  // Argument = are we at the start (true) or end (false) of the data update;
  //            # of particles to prepare for.
  virtual void prepare_data(bool, unsigned);

  // take data for N particles out of the given message, and put it
  // into the proper agency-specific structure.
  // If the Message pointer is 0, put in the locally-stored particles.
  // Arguments: num particles,
  // starting index for inserted particles, Message containing
  // particles
  virtual void insert_data(unsigned, unsigned, Message *);

private:
  // the particle attrib object
  ParticleAttrib<T>& MyParticles;

  // our ACLVIS data structure
  ReadParticleTool *LocalData;
};

#include "DataSource/ACLVISPtclAttribDataSource.hpp"

#endif // ACLVIS_PARTICLE_ATTRIB_DATA_SOURCE_H

/***************************************************************************
 * $RCSfile: ACLVISPtclAttribDataSource.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISPtclAttribDataSource.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
