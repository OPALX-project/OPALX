// ------------------------------------------------------------------------
// $RCSfile: AbstractMapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractMapper
//   This abstract visitor class defines part of the interface for building
//   the transfer map for a beamline.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/AbstractMapper.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"

typedef FTps<double,6> Series;


// Class AbstractMapper
// ------------------------------------------------------------------------

AbstractMapper::AbstractMapper(const Beamline &beamline,
			       const PartData &reference,
			       bool backBeam, bool backTrack):
  DefaultVisitor(beamline, backBeam, backTrack),
  itsReference(reference)
{}


AbstractMapper::~AbstractMapper()
{}


Series AbstractMapper::
buildMultipoleVectorPotential(const BMultipoleField &field)
{
  int order = field.order();

  if (order > 0) {
    static const Series x = Series::makeVariable(X);
    static const Series y = Series::makeVariable(Y);
    Series kx = + field.normal(order) / double(order);
    Series ky = - field.skew(order)   / double(order);
    
    while (order > 1) {
      Series kxt = x * kx - y * ky;
      Series kyt = x * ky + y * kx;
      order--;
      kx = kxt + field.normal(order) / double(order);
      ky = kyt - field.skew(order)   / double(order);
    }

    Series As = x * kx - y * ky;
    As.setTruncOrder(As.getMaxOrder());
    return As;
  } else {
    return Series(0.0);
  }
}


void printSciForm(std::ostream &os, double num, int prec = 14, int fw = 22)
{
  std::streamsize old_prec = os.precision(prec);       // Save old,
  os.setf(std::ios::scientific, std::ios::floatfield); // and set new formats.
  os << std::setw(fw) << num;                          // Print number.
  os.precision(old_prec);                              // Restore old formats.
  os.setf(std::ios::fixed, std::ios::floatfield);
}

Series AbstractMapper::
buildSBendVectorPotential(const BMultipoleField &field, double h)
{
  //std::cerr << "==> In AbstractMapper::buildSBendVectorPotential(const BMultipoleField &field, double h)"
  //          << std::endl;
  int order = field.order();
  Series As; 

  //std::cerr << " h = "; printSciForm(std::cerr,h);
  //std::cerr << "\n" << order << std::endl;
  //std::cerr << std::endl;
  //for (int m = 1; m <= order; ++m) {
  //  std::cerr << "   B(" << m << ") = "; printSciForm(std::cerr, field.normal(m));
  //  std::cerr << "   A(" << m << ") = "; printSciForm(std::cerr, field.skew(m));
  //  std::cerr << std::endl;
  //}

  if (order > 0) {
    static const Series x = Series::makeVariable(X);
    static const Series y = Series::makeVariable(Y);

    // Construct terms constant and linear in y.
    Series Ae = + field.normal(order); // Term even in y.
    Series Ao = - field.skew(order);   // Term odd  in y.

    for (int i = order; --i >= 1; ) {
      Ae = Ae * x + field.normal(i);
      Ao = Ao * x - field.skew(i);
    }
    Ae.setTruncOrder(Ae.getMaxOrder());
    Ao.setTruncOrder(Ao.getMaxOrder());

    Series hx1 = 1. + h * x; // normalized radius
    Ae = + (Ae * hx1).integral(X);
    Ao = - (Ao * hx1);
    // Add terms up to maximum order.
    As = Ae + y * Ao;

    int k = 2;
    if (k <= order) {
      Series yp = y * y / 2.0;

      while (true) {
	// Terms even in y.
	Ae = Ae.derivative(X);
        Ae = h*Ae/hx1 - Ae.derivative(X);
	As += Ae * yp;
	if (++k > order) break;
	yp *= y / double(k);

        // Terms odd in y.
        Ao = Ao.derivative(X);
        Ao = h*Ao/hx1 - Ao.derivative(X);
        As += Ao * yp;
       if (++k > order) break;
       yp *= y / double(k);
      }
    }
  }
  //std::cerr << " As = " << As << std::endl;
  //std::cerr << "==> Leaving AbstractMapper::buildSBendVectorPotential(...)" << std::endl;

  return As;
}
