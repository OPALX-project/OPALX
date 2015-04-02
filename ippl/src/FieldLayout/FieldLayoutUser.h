// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef FIELD_LAYOUT_USER_H
#define FIELD_LAYOUT_USER_H

/***********************************************************************
 * 
 * FieldLayoutUser is a base class for all classes which need to use
 * a FieldLayout - it is derived from User, which provides a virtual
 * function 'notifyUserOfDelete' which is called when the FieldLayout
 * is deleted, and the virtual function 'Repartition' which is called
 * when a Field needs to be redistributed between processors.
 *
 ***********************************************************************/

// include files
#include "Utility/User.h"
#include "Utility/UserList.h"


// class definition
class FieldLayoutUser : public User {

public:
  // constructor - the base class selects a unique ID value
  FieldLayoutUser();

  // destructor, nothing to do here
  virtual ~FieldLayoutUser();

  //
  // virtual functions for FieldLayoutUser's
  //

  // Repartition onto a new layout
  virtual void Repartition(UserList *) = 0;
};

#endif // FIELD_LAYOUT_USER_H

/***************************************************************************
 * $RCSfile: FieldLayoutUser.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: FieldLayoutUser.h,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
