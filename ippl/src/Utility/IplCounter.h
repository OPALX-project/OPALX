// -*- C++ -*-
/***************************************************************************
 *
 * The IPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef IPL_COUNTER_H
#define IPL_COUNTER_H 

/***************************************************************************
 * IplCounter - a simple megaflops counter that accesses hardware counters
 * for measureing megaflop performance.
 *
 * To use these counters:
 *   1. Create a counter.
 *       example:   Counter FFTcounter("FFT");
 *
 *   2. Locate the function which you would like to measure and start the 
 *       counter by placing the startCounter() method before it.  Then, stop
 *       the counter after it by using the stopCounter() method.  
 *       example:   FFTcounter.startCounter();
 *                  fft->transform(....);
 *                  FFTcounter.stopCounter(); 
 *  
 *   3. Use the printIt() method to print to the screen.
 *       example:   FFTcounter.printIt();
 *
 ***************************************************************************/

// include files
#include "Utility/Inform.h"
#include  


class IplCounter
{
public:
  // constructor
  IplCounter(const char *category);

  // destructor
  ~IplCounter();

  // counter operations
  void startCounter();
  void stopCounter();
  void printIt(); 

private:
#if defined(IPL_LONGLONG)
  typedef long long CounterLong;
#else
  typedef long CounterLong;
#endif

  CounterLong totalcyc_m;
  CounterLong totalinst_m;
  CounterLong c0_m;
  CounterLong c21_m;

  int e0_m;
  int e1_m;
  int gen_start_m;
  int gen_read_m;

  std::string category_m;
  Inform msg_m;
};

#endif

/***************************************************************************
 * $RCSfile: IplCounter.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 ***************************************************************************/

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

