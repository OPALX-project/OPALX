// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef ACLVIS_PARTICLE_BASE_DATA_SOURCE_H
#define ACLVIS_PARTICLE_BASE_DATA_SOURCE_H

/***********************************************************************
 * 
 * class ACLVISParticleBaseDataSource
 *
 * ACLVISParticleBaseDataSource is a specific version of ParticleBaseDataSource
 * which takes position and attribute data for a given ParticleBase and
 * provides it to an external agency.  This is done by collecting data on
 * node 0, and formatting it for VTK use.
 *
 * ACLVISParticleBaseDataSource will take the position data, and all connected
 * attributes for the ParticleBase, and update their values.
 *
 ***********************************************************************/

// include files
#include "DataSource/PtclBaseDataSource.h"
#include "Utility/vmap.h"


// forward declarations
template<class PLayout> class ParticleBase;
class ACLVISDataConnect;


// class definition
template<class PLayout>
class ACLVISParticleBaseDataSource : public ParticleBaseDataSource {

public:
  // constructor: name, connection method, transfer method, pbase
  ACLVISParticleBaseDataSource(const char *, DataConnect *, int,
			       ParticleBase<PLayout>&);

  // destructor
  virtual ~ACLVISParticleBaseDataSource();

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
  // ParticleBaseDataSource virtual function interface
  //

  // make a connection using the given attribute.  Return success.
  virtual bool connect_attrib(ParticleAttribDataSource *);

  // disconnect from the external agency the connection involving this
  // particle base and the given attribute.  Return success.
  virtual bool disconnect_attrib(ParticleAttribDataSource *);

  // check to see if the given ParticleAttribBase is in this ParticleBase's
  // list of registered attributes.  Return true if this is so.
  virtual bool has_attrib(ParticleAttribBase *);

private:
  // the set of particles to connect
  ParticleBase<PLayout>& MyParticleBase;

  // recast of DataConnect to ACL-specific object
  ACLVISDataConnect *ACLVISConnection;

  // a map from particle ID's to their serialized index (0 ... totalnum-1).
  // This is updated at the end of each update, so that the following update
  // can tell where a particle moved to.  We actually keep two maps, one with
  // the current values, and the one used to store new values until an update
  // is complete.
  typedef vmap<typename PLayout::Index_t, unsigned int> IDMap_t;
  IDMap_t IDMapA;
  IDMap_t IDMapB;
  IDMap_t *IDMap;
  IDMap_t *NewIDMap;

  // copy data out of the given coordinates iterator into the data structure
  // for the given ParticleAttribDataSource which represents the agency-
  // specific storage.  Arguments = number of particles, starting index,
  // pointer to value, pointer to ID's, and ParticleAttribDataSource
  void insert_pos(unsigned, unsigned, typename PLayout::SingleParticlePos_t *,
		  typename PLayout::Index_t *,
		  ParticleAttribDataSource *);
};

#include "DataSource/ACLVISPtclBaseDataSource.cpp"

#endif // ACLVIS_PARTICLE_BASE_DATA_SOURCE_H

/***************************************************************************
 * $RCSfile: ACLVISPtclBaseDataSource.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISPtclBaseDataSource.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
