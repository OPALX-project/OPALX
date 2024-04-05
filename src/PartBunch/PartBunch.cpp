#include "PartBunch/PartBunch.hpp"


template <>
void PartBunch<double,3>::computeSelfFields() {

        static IpplTimings::TimerRef SolveTimer = IpplTimings::getTimer("SolveTimer");    
        IpplTimings::startTimer(SolveTimer);
        //this->par2grid();
        this->fsolver_m->runSolver();        
        // gather E / B field
        //this->grid2par();
        IpplTimings::stopTimer(SolveTimer);
}
