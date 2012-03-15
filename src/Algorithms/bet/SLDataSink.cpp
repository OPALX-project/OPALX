#ifdef HAVE_ENVELOPE_SOLVER
// ------------------------------------------------------------------------
// $RCSfile: SLDataSink.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SLDataSink
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2004/06/02 19:38:54 $
// $Author: adelmann $
// $Log: SLDataSink.cpp,v $
// Revision 1.3  2004/06/02 19:38:54  adelmann


#include "Algorithms/bet/SLPartBunch.h"
#include "Algorithms/bet/SLDataSink.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "ValueDefinitions/RealVariable.h"
#include <string>

using namespace std;

/////////////////////////////////////////////////////////////////////////////

SLDataSink::SLDataSink()
{

  StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

  // Because we write only once a header into the statfile. This flag gets
  // reset on the first write to the statfile
  firstWriteToStat_m = true;
 
  string fn = OPAL.getInputFn(); 

  INFOMSG(fn<<endl);

  int pos=fn.find(string("."),0);  
  fn.erase(pos,fn.size()-pos);
  
  fn += string(".sdds");

}


SLDataSink::SLDataSink(int restartStep)
{

  StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

  firstWriteToStat_m = true;

  string fn;
  
  if (OPAL.hasRestartFile()){

  }else {
    fn = OPAL.getInputFn();
    int pos=fn.find(string("."),0);  
    fn.erase(pos,fn.size()-pos);
    fn += string(".sdds");
    
  }  
}

/////////////////////////////////////////////////////////////////////////////

SLDataSink::~SLDataSink()
{

    Ippl::Comm->barrier();
}


void SLDataSink::writePhaseSpace(SLPartBunch &beam) {

  beam.calcBeamParameters();

  /**
     
  Write stuff !
  
  
  */

}


void SLDataSink::writeStatData(SLPartBunch &beam,
			    const string fname, 
			    const string el) {
  ofstream os_statData;
  ofstream momentsOf_m;

  IpplTimings::startTimer(StatMarkerTimer_m); 

  unsigned int pwi = 10;

  beam.calcBeamParameters();

  if(Ippl::myNode() == 0) {
    if (firstWriteToStat_m) {
      if (OPAL.inRestartRun()) {
	os_statData.open(fname.c_str(),ios::app);
	os_statData.precision(15);
	os_statData.setf(ios::scientific, ios::floatfield);
      }
      else {
	os_statData.open(fname.c_str(),ios::out);
	os_statData.precision(15);
	os_statData.setf(ios::scientific, ios::floatfield);
      }
      firstWriteToStat_m = false;
    }
    else {
      os_statData.open(fname.c_str(),ios::app);
      os_statData.precision(15);
      os_statData.setf(ios::scientific, ios::floatfield);
    }
  }
  
  if(Ippl::myNode() == 0) {
    
    os_statData << endl;
    os_statData.close();
  }
  
  IpplTimings::stopTimer(StatMarkerTimer_m); 
}
/***************************************************************************
 * $RCSfile: SLDataSink.cpp,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2004/06/02 19:38:54 $
 ***************************************************************************/
#endif
