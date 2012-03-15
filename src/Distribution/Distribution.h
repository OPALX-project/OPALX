#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

// ------------------------------------------------------------------------
// $RCSfile: Distribution.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/bet/SLPartBunch.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Attributes/String.h"
#include "Attributes/Real.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ObjectFunction.h"
#include "Attributes/Real.h"
#include "Attributes/String.h"
#include "Attributes/RealArray.h"

#include "ValueDefinitions/StringConstant.h"
#include "ValueDefinitions/RealVariable.h"

#include "Ippl.h"

#include <hdf5.h>
#include "H5Part.h"

#include <string>

#include "Distribution/Bins.h"

#define RANLIBX
#include "ranlib.h"

#define sqr(x) pow(x,2.)


// Class Distribution
// ------------------------------------------------------------------------
/// The Distribution definition.
//  A Distribution definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class Distribution: public Definition {

public:

  /// Exemplar constructor.
  Distribution();

  virtual ~Distribution();

  /// Test if replacement is allowed.
  //  Can replace only by another Distribution.
  virtual bool canReplaceBy(Object *object);

  /// Make clone.
  virtual Distribution *clone(const string &name);

  /// Check the Distribution data.
  virtual void execute();

  /// Find named Distribution.
  static Distribution *find(const string &name);

  /// Return the embedded CLASSIC PartData.
  const PartData &getReference() const;

  /// Update the Distribution data.
  virtual void update();

  void create(PartBunch &p, size_t Np);
  void create(PartBunch &p, size_t Np, bool scan);

#ifdef HAVE_ENVELOPE_SOLVER
  void createSlicedBunch(SLPartBunch &p);
#endif

  void doRestart(PartBunch &p, size_t Np, size_t restartStep);

  void doRestart_cycl(PartBunch &p, size_t Np, size_t restartStep, const int specifiedNumBunch);

  /// Print the TFS descriptors for the beam.
  void tfsDescriptors(std::ostream &os) const;

  Inform &print(Inform &os) const; 

private:

  void createBinom(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma, 
		   Vector_t bincoef, PartBunch &beam, size_t particles,
		   bool isBinned);
  
  void createUniformTUniformL(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma, 
		                      PartBunch &beam, size_t particles, bool isBinned);
  
  void binnDistribution(PartBunch &beam, size_t Np, string distType);

  double eVtoBetaGamma(const double &valueIneV, const double &mass)
    {
      double tmp = 1. + valueIneV/mass;
      return sqrt(tmp*tmp - 1.);
    }

  // Not implemented.
  Distribution(const Distribution &);
  void operator=(const Distribution &);

  // Clone constructor.
  Distribution(const string &name, Distribution *parent);

  // The particle reference data.
  PartData reference;

  // The structure for particle binning
  PartBins *pbin_m;

  // If we scan we already have a particles and to not have to create them anymore  
  bool scan_m;

};

inline Inform &operator<<(Inform &os, const Distribution &d)
{
  return d.print(os);
}



#endif // OPAL_Distribution_HH
