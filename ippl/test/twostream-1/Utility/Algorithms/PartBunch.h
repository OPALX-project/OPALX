#ifndef CLASSIC_PartBunch_HH
#define CLASSIC_PartBunch_HH

// ------------------------------------------------------------------------
// $RCSfile: PartBunch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class PartBunch
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/Particle.h"
#include <iosfwd>
#include <vector>

template <class T, int, int> class FMatrix;
template <class T, int> class FVector;

// Class PartBunch.
// ------------------------------------------------------------------------
//: Particle Bunch.
//  A representation of a particle bunch as a vector of particles.


class PartBunch: public std::vector<Particle> {

public:

  //: Default constructor.
  //  Construct empty bunch.
  PartBunch();

  //: Conversion.
  PartBunch(const std::vector<Particle> &);

  PartBunch(const PartBunch &);
  ~PartBunch();

  //: Return bunch distribution.
  //  Return the bunch centroid in [b]centroid[/b],
  //  and the second moments in [b]moment[/b].
  void beamEllipsoid(FVector<double,6>   &centroid,
		     FMatrix<double,6,6> &moment) const;

  //: Return maximum amplitudes.
  //  The matrix [b]D[/b] is used to normalise the first two modes.
  //  The maximum normalised amplitudes for these modes are stored
  //  in [b]axmax[/b] and [b]aymax[/b].
  void maximumAmplitudes(const FMatrix<double,6,6> &D,
			 double &axmax, double &aymax) const;
};


std::ostream &operator<<(std::ostream &os, const PartBunch &bunch);

#endif // CLASSIC_PartBunch_HH
