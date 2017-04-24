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

#include "Ippl.h"

#include <vector>
#include <cassert>
#include <sstream>
#include <fstream>

// define our particle layout type                                              
typedef ParticleUniformLayout<double,3> playout_t;

class TracerParticles { // : public IpplParticleBase< ParticleUniformLayout<double, 3> > {

 public:

  TracerParticles() {  
    R.create(1);
    P.create(1);
    Ef.create(1);
    Bf.create(1);

    R[0] = Vector_t(0.0);
    P[0] = Vector_t(0.0);
    Ef[0] = Vector_t(0.0);
    Bf[0] = Vector_t(0.0);
  }

  TracerParticles(size_t n) {  
    R.create(n);
    P.create(n);
    Ef.create(n);
    Bf.create(n);
    for (auto i=0; i<1; i++) { 
      R[i] = Vector_t(0.0);
      P[i] = Vector_t(0.0);
      Ef[i] = Vector_t(0.0);
      Bf[i] = Vector_t(0.0);
    }
  }

  inline size_t size() const { return R.size();}
  inline size_t getLocalNum() const { return size();}

  inline Vector_t getRefR() { return R[0]; }
  inline Vector_t getRefP() { return P[0]; }
  inline Vector_t getRefE() { return Ef[0]; }
  inline Vector_t getRefB() { return Bf[0]; }

  inline void setRefR(Vector_t r) { R[0] = r; }
  inline void setRefP(Vector_t p) { P[0] = p; }
  inline void setRefE(Vector_t e) { Ef[0] = e; }
  inline void setRefB(Vector_t b) { Bf[0] = b; }

  inline Vector_t getTracePartR(size_t i) { return R[i]; }
  inline Vector_t getTracePartP(size_t i) { return P[i]; }
  inline Vector_t getTracePartE(size_t i) { return Ef[i]; }
  inline Vector_t getTracePartB(size_t i) { return Bf[i]; }

  inline void setTracePartR(Vector_t r, size_t i) { R[i] = r; }
  inline void setTracePartP(Vector_t p, size_t i) { P[i] = p; }
  inline void setTracePartE(Vector_t e, size_t i) { Ef[i] = e; }
  inline void setTracePartB(Vector_t b, size_t i) { Bf[i] = b; }


  inline void operator = (const TracerParticles &p ) { 
    assert(p.size() != R.size());
    for (auto i=0; i<1; i++) { 
      R[i] = p.R[i];
      P[i] = p.P[i];
      Ef[i] = p.Ef[i];
      Bf[i] = p.Bf[i];
    }
  }

  ~TracerParticles() { 
    if (of_m)
      of_m.close();
  }

  void openFile();
  void closeFile();
  void writeToFile();
  
  /// make it public to mimic PartBunch
  ParticleAttrib<Vector_t> R;
  ParticleAttrib<Vector_t> P;
  ParticleAttrib<Vector_t> Ef;
  ParticleAttrib<Vector_t> Bf;

 private:

  std::ofstream of_m;
};

#endif
