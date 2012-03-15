#ifndef OPAL_Configure_HH
#define OPAL_Configure_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Configure.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Configure
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


// Namespace Configure
// ------------------------------------------------------------------------
/// The OPAL configurator.
//  This class must be modified to configure the commands to be contained
//  in an executable OPAL-9 program. For each command an exemplar object
//  is constructed and linked to the main directory. This exemplar is then
//  available to the OPAL parser for cloning.
//  This class could be part of the class OpalData.  It is separated from
//  that class and opale into a special module in order to reduce
//  dependencies between modules.

namespace Configure {

  /// Configure all commands.
  extern void configure();
};

#endif // OPAL_Configure_HH
