// ------------------------------------------------------------------------
// $RCSfile: DataSink.hh,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DataSink
//   Original, Observer in the Linac code written by  Tim Cleland, 
//             Julian Cummings, William Humphrey, and Graham Mark
//             Salman Habib and Robert Ryne
//             Los Alamos National Laboratory
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2003/01/23 13:29:44 $
// $Author: adelmann $
// $Log: DataSink.hh,v $
// Revision 1.1.1.1  2003/01/23 13:29:44  adelmann
// Classic
//
// Revision 3.1  2001/08/30 11:27:22  adelmann
// Clean up, remobe .madprogress stuff
//
// Revision 3.0  2001/08/22 14:41:33  adelmann
// The stable Version
//
// Revision 2.17  2001/08/11 05:30:22  adelmann
// Production Version
//
// Revision 1.1.1.1  2000/11/30 20:29:52  adelmann
// g++ and KCC
//
// Revision 1.3  2000/10/26 05:17:59  adelmann
// Remove DX stuff and add timestam and title
//
// Revision 1.2  2000/08/10 10:56:35  adelmann
// Some cleanup and add a option dx (data explorer) !
// The new printall skript extracts data from this stat file format
//
// Revision 1.1.1.1  2000/07/14 07:20:54  adelmann
// linux version Fri Jul 14 09:15:27 CEST 2000
//
// Revision 1.1.1.1  2000/05/20 11:13:58  adelmann
// Initial working version, thick elements, without: field dump and collimators
//
// Revision 1.3  2000/01/28 07:22:42  adelmann
// Fixrd some bugs with the Mad9pOutput
//
// Revision 1.2  2000/01/27 14:13:29  adelmann
// - Add  bunch->dataSink_m.saveStatDataGnuplotFormat( . )
//        DiscParticle write out
//        updateDotProgress
//
// Revision 1.1.1.1  2000/01/06 07:33:27  adelmann
// linux version works with gcc 991007 V2.96
//
// Revision 2.2  1999/10/29 05:02:07  adelmann
// *** empty log message ***
//
// Revision 2.1  1999/10/27 06:37:00  adelmann
// SGI-LINUX g++, with RF-Gap, REVSCATTER and read distribution from file
//
// Revision 1.1.1.1  1999/10/26 04:22:18  adelmann
// Classic 2.1 (with p)
//
// Revision 1.1.1.1  1999/10/26 04:14:36  adelmann
// Classic 2.1 (with p)
//
//
// ------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////////
// Class: DataSink
// Observer generates diagnostic output of the accelerator beam.
// The class computes statistical descriptors of particles' positions
// and momemta and writes those in files.
// It also writes particles' positions and momenta to files.
// prints initial and final particle positions into files
// accumulates beam statistics in Inform objects (files) during run
/////////////////////////////////////////////////////////////////////////////

#ifndef DataSink_H_
#define DataSink_H_

#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunch.h"
#include "Ippl.h"

#include <hdf5.h>
#include "H5Part.h"

/////////////////////////////////////////////////////////////////////////////

class DataSink
{
public:
  DataSink();    
  DataSink(int restartStep);    

  ~DataSink();

private:
 
  DataSink( const DataSink & ) { }
  DataSink & operator = ( const DataSink & ) { return *this; }

public:
 

  /** Write Stat Data 
  */
  void writeStatData(PartBunch &beam,
		     const string fname, 
		     const string el);

  //FDext = {BHead, EHead, BRef, ERef, BTail, ETail}
  /** Dumps Phase Space to H5 file
   *  \param beam the beam 
   *  \param FDext the external E and B field for head, ref and tail particle. The vector array has the following layout: 
   *    - FDext[0] = BHead
   *    - FDext[1] = EHead
   *    - FDext[2] = BRef
   *    - FDext[3] = ERef
   *    - FDext[4] = BTail
   *    - FDext[5] = ETail
   *  \param sposHead position of the head particle
   *  \param sposRef position of the reference particle 
   *  \param sposTail position of the reference particle
   *
   */
  void writePhaseSpace(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);
  void writePhaseSpace_cycl(PartBunch &beam, Vector_t FDext[]);
private:
      
  bool firstWriteToStat_m;
  bool firstWriteH5part_m;

  H5PartFile *H5file_m;
  h5part_int64_t H5call_m;

  IpplTimings::TimerRef StatMarkerTimer_m;
  IpplTimings::TimerRef H5PartTimer_m;
  
  // to compensate IMPACT-T wrong spos
  double sshift_m;

};

#endif // DataSink_H_

/***************************************************************************
 * $RCSfile: DataSink.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 13:29:44 $
 ***************************************************************************/

