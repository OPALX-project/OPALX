// ------------------------------------------------------------------------
// $RCSfile: PartBunch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class PartBunch
//   Interface to a particle bunch.
//   Can be used to avoid use of a template in user code.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/PartBunch.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <iostream>


// Class PartBunch
// ------------------------------------------------------------------------

PartBunch::PartBunch():
  std::vector<Particle>()
{}


PartBunch::PartBunch(const PartBunch &rhs):
  std::vector<Particle>(rhs)
{}


PartBunch::PartBunch(const std::vector<Particle> &rhs):
  std::vector<Particle>(rhs)
{}


PartBunch::~PartBunch()
{}


void PartBunch::beamEllipsoid(FVector<double,6>   &centroid,
			      FMatrix<double,6,6> &moment) const
{
  for (int i = 0; i < 6; ++i) {
    centroid(i) = 0.0;
    for (int j = 0; j <= i; ++j) {
      moment(i,j) = 0.0;
    }
  }

  PartBunch::const_iterator last = end();
  for (PartBunch::const_iterator part = begin(); part != last; ++part) {
    for (int i = 0; i < 6; ++i) {
      centroid(i) += (*part)[i];
      for (int j = 0; j <= i; ++j) {
	moment(i,j) += (*part)[i] * (*part)[j];
      }
    }
  }
  
  double factor = 1.0 / double(size());
  for (int i = 0; i < 6; ++i) {
    centroid(i) *= factor;
    for (int j = 0; j <= i; ++j) {
      moment(j,i) = moment(i,j) *= factor;
    }
  }
}


void PartBunch::maximumAmplitudes(const FMatrix<double,6,6> &D,
				  double &axmax, double &aymax) const
{
  axmax = aymax = 0.0;
  PartBunch::const_iterator last = end();

  for (PartBunch::const_iterator part = begin(); part != last; ++part) {
    double xnor =
      D(0,0)*part->x()  + D(0,1)*part->px() + D(0,2)*part->y() +
      D(0,3)*part->py() + D(0,4)*part->t()  + D(0,5)*part->pt();
    double pxnor =
      D(1,0)*part->x()  + D(1,1)*part->px() + D(1,2)*part->y() +
      D(1,3)*part->py() + D(1,4)*part->t()  + D(1,5)*part->pt();
    double ynor =
      D(2,0)*part->x()  + D(2,1)*part->px() + D(2,2)*part->y() +
      D(2,3)*part->py() + D(2,4)*part->t()  + D(2,5)*part->pt();
    double pynor =
      D(3,0)*part->x()  + D(3,1)*part->px() + D(3,2)*part->y() +
      D(3,3)*part->py() + D(3,4)*part->t()  + D(3,5)*part->pt();

    axmax = std::max(axmax, (xnor*xnor + pxnor*pxnor));
    aymax = std::max(aymax, (ynor*ynor + pynor*pynor));
  }
}


std::ostream &operator<<(std::ostream &os, const PartBunch &bunch)
{
  std::streamsize old_prec = os.precision(15);
  std::ios::fmtflags old_flag = os.setf(std::ios::scientific, std::ios::floatfield);
  PartBunch::const_iterator last = bunch.end();

  for (PartBunch::const_iterator part = bunch.begin(); part != last; ++part) {
    os << part->x() << " " << part->px() << " "
       << part->y() << " " << part->py() << " "
       << part->t() << " " << part->pt() << std::endl;
  }
  
  os.precision(old_prec);
  os.flags(old_flag);
  return os;
}
