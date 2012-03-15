#ifndef OPAL_FILTER_HH
#define OPAL_FILTER_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalFilter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalFilter
//
// ------------------------------------------------------------------------
//
// $Date: 2008/10/14 09:33:44 $
// $Author: Christof Kraus $
//
// ------------------------------------------------------------------------

#include <vector>
#include "AbstractObjects/Definition.h"
#include "Filters/Filter.h"

// Class OpalFilter
// ------------------------------------------------------------------------
/// The FILTER definition.
//  A FILTER definition is used to define a filter which can be applied
//  to a 1D histogram in order to get rid of noise.

class OpalFilter: public Definition {

public:

  /// Exemplar constructor.
  OpalFilter();

  virtual ~OpalFilter();

  /// Test if replacement is allowed.
  //  Can replace only by another OpalFilter
  virtual bool canReplaceBy(Object *object);

  /// Make clone.
  virtual OpalFilter *clone(const string &name);

  /// Check the OpalFilter data.
  virtual void execute();

  /// Find named FILTER.
  static OpalFilter *find(const string &name);

  /// Update the OpalFilter data.
  virtual void update();

  /// Print the TFS descriptors for the Filter
  void tfsDescriptors(std::ostream &os) const;

  Inform &print(Inform &os) const; 

  void initOpalFilter();

  inline void apply(vector<double> &histogram);
  inline void calc_derivative(vector<double> &histogram, const double &hz);

  Filter *filter_m;
private:

  // Not implemented.
  OpalFilter(const OpalFilter &);
  void operator=(const OpalFilter &);

  // Clone constructor.
  OpalFilter(const string &name, OpalFilter *parent);

};

void OpalFilter::apply(vector<double> &histogram)
{
  if (filter_m)
    filter_m->apply(histogram);
}

void OpalFilter::calc_derivative(vector<double> &histogram, const double &hz)
{
  if (filter_m)
    filter_m->calc_derivative(histogram, hz);
}

inline Inform &operator<<(Inform &os, const OpalFilter &b)
{
  return b.print(os);
}

#endif // OPAL_FILTER_HH
