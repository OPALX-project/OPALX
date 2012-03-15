#ifndef MAD_Timer_HH
#define MAD_Timer_HH 1
// ------------------------------------------------------------------------
// $RCSfile: GTTimer.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTUtilities/GTTimer.hh,v 1.1.1.1 2003/01/23 09:13:58 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// ------------------------------------------------------------------------
// Class category: 
// ------------------------------------------------------------------------
//
// $Date: 2003/01/23 09:13:58 $
// $Author: adelmann $
// $Log: GTTimer.hh,v $
// Revision 1.1.1.1  2003/01/23 09:13:58  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#include <ctime>
#include <string>
using namespace std;

// Class Timer
// ------------------------------------------------------------------------

class MTimer {

public:

  // Construction/destruction.
  MTimer();
  ~MTimer();

  // Return date.
  string date() const;

  // Return time.
  string time() const;

private:

  // Not implemented.
  MTimer(const MTimer &);
  void operator=(const MTimer &);

  time_t timer;
};

#endif // MAD_Timer_HH
