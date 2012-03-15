#ifdef HAVE_ENVELOPE_SOLVER
#ifndef OPAL_SLPartBunch_HH
#define OPAL_SLPartBunch_HH

// ------------------------------------------------------------------------
// $Rfile: SLPartBunch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class SLPartBunch
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: Andreas Adelmann  and Co. $
//
// ------------------------------------------------------------------------

#include "Algorithms/Particle.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"

#include <iosfwd>
#include <vector>
#include "Ippl.h"

class SLPartBunch: public PartBunch {
  const PartData *reference;

 public:
  /// Default constructor.
  //  Construct empty bunch.
  SLPartBunch(const PartData *ref);

  /// Conversion.
  SLPartBunch(const std::vector<Particle> &,const PartData *ref);

  SLPartBunch(const SLPartBunch &);
  ~SLPartBunch();


  /**
     Comes the slice data
     
  */
  unsigned int LastSection[1000];  // last em-field section
  double                dt[1000];  // eigen-time of the slice
  Vector_t               Z[1000];  // centroid in z of the slice

  Inform &print(Inform &os);

  void calcBeamParameters() { INFOMSG("HELLO SL calcBeamParameters" << endl);}

  int getLocalNumSlices() { return 1;}
  int getTotalNumSlices() { return getLocalNumSlices();}

};

inline Inform &operator<<(Inform &os, SLPartBunch &p)
{
  return p.print(os);
}

#endif // OPAL_SLPartBunch_HH
#endif
