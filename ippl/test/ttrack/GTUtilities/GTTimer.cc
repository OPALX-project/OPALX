// ------------------------------------------------------------------------
// $RCSfile: GTTimer.cc,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Timer
//   Get Calendar date and time.
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2003/01/23 09:13:58 $
// $Author: adelmann $
// $Log: GTTimer.cc,v $
// Revision 1.1.1.1  2003/01/23 09:13:58  adelmann
// GenTrackE
//
// Revision 3.0  2001/08/22 14:47:04  adelmann
// The stable Version
//
// Revision 1.1.1.1  2000/11/30 20:34:18  adelmann
// g++ and KCC
//
// Revision 1.1.1.1  2000/07/14 07:22:00  adelmann
// linux version Fri Jul 14 09:15:27 CEST 2000
//
// Revision 1.1.1.1  2000/05/20 11:12:50  adelmann
// Initial working version without: field dump and collimators
//
// Revision 1.1.1.1  2000/01/06 07:32:32  adelmann
// linux version works with gcc 991007 V2.96
//
// Revision 2.1  1999/10/27 07:31:04  adelmann
// SGI-LINUX g++, with RF-Gap, REVSCATTER and read distribution from file
//
// Revision 1.1.1.1  1999/10/25 14:04:11  adelmann
// MAD 9.2p
//
// Revision 9.1  1999/06/04 04:26:54  adelmann
// Name clash with POOMA
//
// Revision 9.21  1999/03/08 06:18:30  adelmann
// *** empty log message ***
//
// Revision 1.1.1.1  1999/03/02 15:07:30  fci
// release 9.2/0
//
// ------------------------------------------------------------------------

#include "GTTimer.hh"
#include <ctime>


// Class Timer
// ------------------------------------------------------------------------

MTimer::MTimer()
{
  ::time(&timer);
}


MTimer::~MTimer()
{}


string MTimer::date() const
{
  char buffer[12];
  strftime(buffer, 12, "%d/%m/%Y", localtime(&timer));
  return string(buffer, 10);
}


string MTimer::time() const
{
  char buffer[12];
  strftime(buffer, 12, "%H.%M.%S", localtime(&timer));
  return string(buffer, 8);
}
