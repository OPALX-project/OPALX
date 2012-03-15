#ifndef OPAL_MatchLimits_HH
#define OPAL_MatchLimits_HH 1

// ------------------------------------------------------------------------
// $RCSfile: MatchLimits.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Enum MatchLimits:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


// Enum MatchLimits
// Possible limits for matching variables and constraints.
// Note: LOWER_LIMIT + UPPER_LIMIT == BOTH_LIMITS.
// ------------------------------------------------------------------------
/// The possible limits used in a matching process.

enum {
  NO_LIMIT    = 0,      // No limits present.
  LOWER_LIMIT = 1,      // Lower limit only.
  UPPER_LIMIT = 2,      // Upper limit only.
  BOTH_LIMITS = 3,      // Interval limits.
  EQUALITY    = 4       // Equality prescribed.
};

#endif // OPAL_MatchLimits_HH

