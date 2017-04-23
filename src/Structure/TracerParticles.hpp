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
#include <sstream>
#include <fstream>

class TracerParticles {

 public:

  TracerParticles() {  
    tracePartR_m.reserve(1);
    tracePartP_m.reserve(1);
    tracePartE_m.reserve(1);
    tracePartB_m.reserve(1);

    tracePartR_m[0] = Vector_t(0.0);
    tracePartP_m[0] = Vector_t(0.0);
    tracePartE_m[0] = Vector_t(0.0);
    tracePartB_m[0] = Vector_t(0.0);
  }

  TracerParticles(size_t n) {  
    tracePartR_m.reserve(n);
    tracePartP_m.reserve(n);
    tracePartE_m.reserve(n);
    tracePartB_m.reserve(n);
    for (auto i=0; i<1; i++) { 
      tracePartR_m[i] = Vector_t(0.0);
      tracePartP_m[i] = Vector_t(0.0);
      tracePartE_m[i] = Vector_t(0.0);
      tracePartB_m[i] = Vector_t(0.0);
    }
  }

  inline size_t size() const { return tracePartR_m.size();}

  inline Vector_t getRefR() { return tracePartR_m[0]; }
  inline Vector_t getRefP() { return tracePartP_m[0]; }
  inline Vector_t getRefE() { return tracePartE_m[0]; }
  inline Vector_t getRefB() { return tracePartB_m[0]; }

  inline void setRefR(Vector_t r) { tracePartR_m[0] = r; }
  inline void setRefP(Vector_t p) { tracePartP_m[0] = p; }
  inline void setRefE(Vector_t e) { tracePartE_m[0] = e; }
  inline void setRefB(Vector_t b) { tracePartB_m[0] = b; }

  inline Vector_t getTracePartR(size_t i) { return tracePartR_m[i]; }
  inline Vector_t getTracePartP(size_t i) { return tracePartP_m[i]; }
  inline Vector_t getTracePartE(size_t i) { return tracePartE_m[i]; }
  inline Vector_t getTracePartB(size_t i) { return tracePartB_m[i]; }

  inline void setTracePartR(Vector_t r, size_t i) { tracePartR_m[i] = r; }
  inline void setTracePartP(Vector_t p, size_t i) { tracePartP_m[i] = p; }
  inline void setTracePartE(Vector_t e, size_t i) { tracePartE_m[i] = e; }
  inline void setTracePartB(Vector_t b, size_t i) { tracePartB_m[i] = b; }


  inline void operator = (const TracerParticles &p ) { 
    assert(p.size() != tracePartR_m.size());
    for (auto i=0; i<1; i++) { 
      tracePartR_m[i] = p.tracePartR_m[i];
      tracePartP_m[i] = p.tracePartP_m[i];
      tracePartE_m[i] = p.tracePartE_m[i];
      tracePartB_m[i] = p.tracePartB_m[i];
    }
  }

  ~TracerParticles() { 
    if (of_m)
      of_m.close();
  }

  void openFile();
  void closeFile();
  void writeToFile();
  
 private:

  std::vector<Vector_t> tracePartR_m;
  std::vector<Vector_t> tracePartP_m;
  std::vector<Vector_t> tracePartE_m;
  std::vector<Vector_t> tracePartB_m;

  std::ofstream of_m;
};

#endif
