#ifndef OPAL_MatchState_HH
#define OPAL_MatchState_HH 1

// ------------------------------------------------------------------------
// $RCSfile: MatchState.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Enum MatchState:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------


// Enum MatchState
// ------------------------------------------------------------------------
/// The possible states of a matching process.
//  Flag values transmitted to the match function.
//  [b]NOTE:[/b] The order of definition of the constants matters for
//  control of evaluation and printing.

enum MatchState {    // State                  print for
    INTERNAL,          // Internal call          level > 2
    PROGRESS,          // New minimum found      level > 1
    START,             // First call             level > 0
    RESTART,           // Algorithm restarted    level > 0
    CHECK,             // Check required         level > 0
    CONVERGED,         // Match has converged    level > 0
    FAILED,            // Match failed           level > 0
    CALL_LIMIT,        // Call limit exceeded    level > 0
    ACCURACY_LIMIT,    // Tolerance too small    level > 0
    TERMINATED         // Terminated by user     level > 0
};

#endif // OPAL_MatchState_HH
