// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef ACLVIS_OPERATIONS_H
#define ACLVIS_OPERATIONS_H

/***********************************************************************
 * 
 * This file contains various ACLVIS-specific helper classes and templated
 * functions, which are used to convert templated objects into type-specific
 * ACLVIS structures.
 *
 ***********************************************************************/

// include files
#ifdef IPPL_LUX
# include "luxvis.h"
#else
# include "aclvis/aclvis.h"
#endif
#include "AppTypes/Vektor.h"
#include "AppTypes/dcomplex.h"


////////////////////////////////////////////////////////////////////////////
// simple templated class with specializations to convert data to float
// if possible, or else to just return 0
struct ACLVISToFloat {
  template<class T>
  float operator()(const T&)            { return 0.0; }

  float operator()(const char& val)     { return static_cast<float>(val); }
  float operator()(const short& val)    { return static_cast<float>(val); }
  float operator()(const int& val)      { return static_cast<float>(val); }
  float operator()(const long& val)     { return static_cast<float>(val); }
  float operator()(const float& val)    { return static_cast<float>(val); }
  float operator()(const double& val)   { return static_cast<float>(val); }
  float operator()(const dcomplex& val) {
    return static_cast<float>(norm(val));
  }
};


////////////////////////////////////////////////////////////////////////////
// a templated traits class which is used to determine what type of data is
// being displayed, and to convert the data into the format needed by the
// ACLVIS tools.  By default, it works for scalar data.
template<class V, class T>
struct ACLVISTraits {
  // return what type of data this is
#ifdef IPPL_LUX
  static int getType() { return vizStructuredFieldDataType::LUXVIS_SCALAR; }
#else
  static int getType() { return vizStructuredFieldDataType::ACLVIS_SCALAR; }
#endif

  // set the value of the Nth point
  static void setPoint(V *vizdata, int indx, const T &val) {
    vizdata->GetVizData()->AddScalar(indx, ACLVISToFloat()(val));
  }

  // set the value of the Nth ID ... this should only be defined in this
  // version of the traits class
  static void setID(V *vizdata, int indx, int id) {
    vizdata->GetVizData()->AddIdInfoVal(indx, id);
  }
};


////////////////////////////////////////////////////////////////////////////
// a version of the Traits class for Vektors
template<class V, class T, unsigned int D>
struct ACLVISTraits<V, Vektor<T,D> > {
  // return what type of data this is
#ifdef IPPL_LUX
  static int getType() { return vizStructuredFieldDataType::LUXVIS_VECTOR; }
#else
  static int getType() { return vizStructuredFieldDataType::ACLVIS_VECTOR; }
#endif

  // set the value of the Nth point
  static void setPoint(V *vizdata, int indx, const Vektor<T,D> &val) {
    float vect[3];
    for (unsigned int d=0; d < 3; ++d)
      vect[d] = (d < D ? ACLVISToFloat()(val[d]) : 0.0);
    vizdata->GetVizData()->AddVector(indx, vect);
  }

  // set the value of the Nth particle coordinate ... this should only
  // be defined in this version of the traits class
  static void setCoord(V *vizdata, int indx, const Vektor<T,D> &pos) {
    float vect[3];
    for (unsigned int d=0; d < 3; ++d)
      vect[d] = (d < D ? ACLVISToFloat()(pos[d]) : 0.0);
    vizdata->GetVizData()->AddPoint(indx, vect);
  }
};


#endif // ACLVIS_OPERATIONS_H

/***************************************************************************
 * $RCSfile: ACLVISOperations.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: ACLVISOperations.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
