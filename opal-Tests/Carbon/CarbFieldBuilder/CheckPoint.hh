#ifndef OPAL_checkPoint_HH
#define OPAL_checkPoint_HH 

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

namespace checkPoint {

  /// Configure all commands.
  // extern int isLeft( Point P0, Point P1, Point P2 );
  extern int cn_PnPoly( Point P, Point* V, int n );
  extern int wn_PnPoly( Point P, Point* V, int n );
    
};

#endif // OPAL_Configure_HH
