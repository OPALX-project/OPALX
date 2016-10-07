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
#include "Utility/DiscParticle.h"
#include "Utility/DiscConfig.h"
#include "Utility/PAssert.h"
#include "Particle/IpplParticleBase.h"
#include "Particle/ParticleAttrib.h"
#include "Message/Communicate.h"
#include "Message/Tags.h"

#include "Utility/IpplInfo.h"

// debugging macros; MWERKS: put these at top of .h files.

///////////////////////////////////////////////////////////////////////////
//                           READ METHODS
// read the specifed record in the file into the given IpplParticleBase or
// ParticleAttrib object, depending on how the DiscParticle was created.
// If the method is to read all the IpplParticleBase, this will delete all the
// existing particles in the given object, create new ones and store the
// values, and then do an update.  If an attribute is being read, this
// will only work if the number of particles in the attribute already
// matches the number in the file.
///////////////////////////////////////////////////////////////////////////
//MWERKS Moved into class definition (.h file).

///////////////////////////////////////////////////////////////////////////
//                           WRITE METHODS
// write the data from the given IpplParticleBase or ParticleAttrib into the
// file.  Data is appended as a new record, and the meta file is updated.
//
///////////////////////////////////////////////////////////////////////////
//MWERKS Moved into class definition (.h file).

/***************************************************************************
 * $RCSfile: DiscParticle.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: DiscParticle.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/

