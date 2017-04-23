#ifndef OPAL_TRACERPARTICLES_HH
#define OPAL_TRACERPARTICLES_HH
/** 
 *  @file    TracerPartiles.h
 *  @author  Andreas Adelmann
 *  @date    27/4/2017  
 *  @version 
 *  
 *  @brief Holds a set of tracer particles on every core
 *
 *  @section DESCRIPTION
 *  
 *  A set of non self interacting particles, that are
 *  subject of external fields. Every core holds them.
 *  At the moment they are never lost either and only one
 *  particle is stored. 
 */
#include "Algorithms/PartBunch.h"
#include "Algorithms/Tracker.h"
#include "Algorithms/PartPusher.h"
#include "Structure/DataSink.h"

#include <vector>
#include <cassert>
class TracerParticles {

 public:

  TracerParticles() {  
    RefPartR_m.reserve(1);
    RefPartP_m.reserve(1);

    RefPartR_m[0] = Vector_t(0.0);
    RefPartP_m[0] = Vector_t(0.0);
  }
  inline size_t size() const { return RefPartR_m.size();}

  inline Vector_t getRefR() { return RefPartR_m[0]; }
  inline Vector_t getRefP() { return RefPartP_m[0]; }

  inline void setRefR(Vector_t r) { RefPartR_m[0] = r; }
  inline void setRefP(Vector_t p) { RefPartP_m[0] = p; }

  inline void operator = (const TracerParticles &p ) { 
    assert(p.size() != RefPartR_m.size());
    for (auto i=0; i<1; i++) { 
      RefPartR_m[i] = p.RefPartR_m[i];
      RefPartP_m[i] = p.RefPartP_m[i];
    }
  }

  ~TracerParticles() { }

 private:

  std::vector<Vector_t> RefPartR_m;
  std::vector<Vector_t> RefPartP_m;

};

#endif
