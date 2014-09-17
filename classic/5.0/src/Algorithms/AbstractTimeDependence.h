#ifndef _CLASSIC_SRC_ALGORITHMS_ABSTRACTTIMEDEPENDENCE_H_
#define _CLASSIC_SRC_ALGORITHMS_ABSTRACTTIMEDEPENDENCE_H_

#include <string>
#include <map>

#include "Utilities/OpalException.h"

/** Time dependence abstraction for field parameters that vary slowly with time
 */
class AbstractTimeDependence {
  // Unit tests for this class are in with PolynomialTimeDependence
  public:
    virtual ~AbstractTimeDependence() {}
    virtual AbstractTimeDependence* clone() = 0;
    virtual double getValue(double time) = 0;

    static AbstractTimeDependence* getTimeDependence(std::string name);
    static void setTimeDependence(std::string name,
                                  AbstractTimeDependence* time_dep);

    /** Reverse the mapping - slow
     */
    static std::string getName(AbstractTimeDependence* time_dep);

  private:
    static std::map<std::string, AbstractTimeDependence*> td_map;
};

#endif

