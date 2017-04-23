#include "Structure/TracerParticles.hpp"

void TracerParticles::openFile() {
  /// only core 0 is writing
  if (Ippl::myNode()==0) {
    of_m.open("data/tracer.dat", std::ios::out);
    of_m.precision(15);
    of_m.setf(std::ios::scientific, std::ios::floatfield);
    of_m << "# id, x, px, Ex, Bx, y, py, Ey, By, z, pz, Ez, Bz " << std::endl;
  }
} 

void TracerParticles::closeFile() { 
  if (of_m)
    of_m.close();
}

void TracerParticles::writeToFile() { 
  /// only core 0 is writing
  if ((Ippl::myNode() == 0) && (of_m)) {
    for (unsigned i=0; i<tracePartR_m.size(); i++) {
      of_m << i << "\t";
      for (auto d=0; d<3; d++)
	of_m << tracePartR_m[i](d) << "\t" << tracePartP_m[i](d) << "\t"
	     << tracePartE_m[i](d) << "\t" << tracePartB_m[i](d) << "\t" ;
    }
    of_m << std::endl;
  } 
}
