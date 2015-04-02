// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * CVS INFORMATION:
 * $Id
 ***************************************************************************
 * DESCRIPTION:
 *	The Timer class allows for easy timing of the program.  The timer
 * tracks real (clock) time elapsed, user time, and system time.
 *
 ***************************************************************************/

#ifndef TIMER_H
#define TIMER_H

#ifdef __sgi
// make sure this is defined for BSD time routines
#define _BSD_TYPES
// fix a glitch in ANSI compatibility with SGI headers
#define _STAMP_T
#endif

#ifndef __MWERKS__
// For now, stub out all Timer guts for MetroWerks

#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#endif // __MWERKS__

#ifdef __sgi
// fix a glitch in ANSI compatibility with SGI headers
#undef _STAMP_T
#endif


class Timer
{
public:
  Timer();			// Constructor
  ~Timer();                     // Destructor
  void clear();			// Set all accumulated times to 0
  void start();			// Start timer
  void stop();			// Stop timer

  double clock_time();		// Report clock time accumulated in seconds
  double user_time();		// Report user time accumlated in seconds
  double system_time();		// Report system time accumulated in seconds
  double cpu_time()
  {
    // Report total cpu_time which is just user_time + system_time
    return ( user_time() + system_time() );
  }		

#ifndef __MWERKS__
// For now, stub out all Timer guts for MetroWerks

  double calibration;		// Calibration time: time it takes to
                                // get in and out of timer functions
private:
  short timer_state;		// State of timer, either on or off
#ifdef IPPL_XT3
  double last_clock;            // Previous value from dclock()
  double current_clock;         // Current value from dclock()
#else
  long cpu_speed;		  // CPU speed for times() call

  unsigned long last_secs;	  // Clock seconds value when the
				  // timer was last started
  long last_usecs;		  // Clock useconds value when the
				  // timer was last started
  unsigned long last_user_time;   // User time when timer was last started
  unsigned long last_system_time; // System time when timer was last started

  long current_secs;		// Current accumulated clock seconds
  long current_usecs;		// Current accumulated clock useconds
  long current_user_time;	// Current accumulated user time
  long current_system_time;	// Current accumulated system time
#endif

#if ( defined(IPPL_T3E) || defined(IPPL_XT3) )
  // don't need these structs
#else
  struct tms tmsbuf;	        //  Values from call to times
  struct timeval tvbuf;	        //  Values from call to gettimeofday
  struct timezone tzbuf;        //  Timezone values from gettimeofday
	  		        //  These values aren't used for anything
#endif

#endif // __MWERKS__
};

#endif // TIMER_H

/***************************************************************************
 * $RCSfile: Timer.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: Timer.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
