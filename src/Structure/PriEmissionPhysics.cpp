#include "Structure/PriEmissionPhysics.h"
//FIXME: cleanup
#include "Algorithms/PartBunch.h"
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"
#include "Elements/OpalBeamline.h"
#include "Distribution/Distribution.h"
#include "Ippl.h"
#include "hdf5.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <sys/stat.h>



PriEmissionPhysics::PriEmissionPhysics() {

}


PriEmissionPhysics::~PriEmissionPhysics() {

}


/**
 * Field emission model.
 */


void PriEmissionPhysics::Fieldemission(PartBunch *itsBunch, const double &fa, const double &Enormal, const double &parameterFNB,
                                       const double &workFunction, const double &parameterFNVYZe, const double &parameterFNVYSe,
                                       const double &parameterFNY, const double &fieldEnhancement, const double &maxFNemission,
                                       const double &TriArea, const std::vector<Vector_t> &vertex, const Vector_t TriNormal, size_t &Nstp)  {

    int node_num = Ippl::getNodes();

    size_t *count = new size_t [node_num];

    for(int i = 0; i < node_num; i++) {

        count[i] = 0;

    }


    // Self-field is not considered at moment. Only 1D Child-Langmuir law is implemented for space charge limited current density.
   
    double Jtmp = fa * Enormal * Enormal * exp(- parameterFNB *  workFunction * sqrt(workFunction) * (parameterFNVYZe -  parameterFNVYSe * pow(parameterFNY * sqrt(fieldEnhancement * (-Enormal)) /  workFunction, 2)) /  fieldEnhancement / (-Enormal));    // Fowler-Nordheim model;
    double tmpd = 0.5 * 1.7588e11 * (-Enormal) * itsBunch->getdT() * itsBunch->getdT();// assume particle move 0.5*a*t*t
    double Jchild = 2.33395e-6 * pow((-Enormal), 1.5) / sqrt(tmpd);// Child-Langmuir model;

    if(Jchild < Jtmp) {
        /* We can ignore this if we have space charge solver which could deal with complicated boundary geometry.
         * The space charge effect will limit the surface field and get self-consistent
         * Fowler-Nordheim emission model and in that situation we will not need Child-Langmuir model anymore.
         */

        Jtmp = Jchild;
        double chargeScalFactor = 1.0;
        size_t N = static_cast<size_t>(Jtmp * TriArea / -itsBunch->getChargePerParticle() * itsBunch->getdT());

        if(N >  maxFNemission) {
            chargeScalFactor = static_cast<double>(N / maxFNemission);
            N =  maxFNemission;
        }
        Nstp += N;
        int pc = 0;
        for(size_t k = 0; k < N; k++) {

            if(pc == Ippl::myNode()) {

                if(count[pc] == 0) {
                    count[pc] = itsBunch->getLocalNum();

                }

                double r1 = IpplRandom();
                double r2 = IpplRandom();

                itsBunch->create(1);
                itsBunch->R[count[pc]] = (vertex[0] + r1 * (vertex[1] -  vertex[0]) + r2 * (1 - r1) * (vertex[2] -  vertex[0])) +  TriNormal * tmpd ;
                itsBunch->P[count[pc]] = Vector_t(0.0);
                itsBunch->Bin[count[pc]] = 0;
                itsBunch->PType[count[pc]] = 1;
                itsBunch->TriID[count[pc]] = 0;
                itsBunch->Q[count[pc]] = chargeScalFactor * itsBunch->getChargePerParticle();
                itsBunch->LastSection[count[pc]] = 0;
                itsBunch->Ef[count[pc]] = Vector_t(0.0);
                itsBunch->Bf[count[pc]] = Vector_t(0.0);
                itsBunch->dt[count[pc]] = itsBunch->getdT();
                itsBunch->LastSection[count[pc]] = 0;

                count[pc]++;

            }
            pc++;
            if(pc == Ippl::getNodes())//nodes number should be equal in BoundaryGeometry and PriemissionPhysics.
                pc = 0;

        }
    } else {
        double chargeScalFactor = 1.0;
        size_t N = static_cast<size_t>(Jtmp *  TriArea / -itsBunch->getChargePerParticle() * itsBunch->getdT());

        if(N >  maxFNemission) {

            chargeScalFactor = static_cast<double>(N /  maxFNemission);
            N =  maxFNemission;

        }
        Nstp += N;
        int pc = 0;

        for(size_t k = 0; k < N; k++) {

            if(pc == Ippl::myNode()) {

                if(count[pc] == 0) {
                    count[pc] = itsBunch->getLocalNum();
                }

                double r1 = IpplRandom();
                double r2 = IpplRandom();

                itsBunch->create(1);
                // Here r1+r2*(1-r1)<=1 should be true;
                itsBunch->R[count[pc]] =  vertex[0] + r1 * (vertex[1] -  vertex[0]) + r2 * (1 - r1) * (vertex[2] -  vertex[0]) +  TriNormal * tmpd;  // During our simulation, we find some emitted particles will not move if they just on the surface. So here we also add margins =  TriNormal * tmpd from the emission surface to make sure that the emitted particles will "feel" the RF field.
                itsBunch->P[count[pc]] = Vector_t(0.0);
                itsBunch->Bin[count[pc]] = 0;
                itsBunch->PType[count[pc]] = 1;
                itsBunch->TriID[count[pc]] = 0;
                itsBunch->Q[count[pc]] = chargeScalFactor * itsBunch->getChargePerParticle();
                itsBunch->LastSection[count[pc]] = 0;
                itsBunch->Ef[count[pc]] = Vector_t(0.0);
                itsBunch->Bf[count[pc]] = Vector_t(0.0);
                itsBunch->dt[count[pc]] = itsBunch->getdT();
                itsBunch->LastSection[count[pc]] = 0;
                count[pc]++;

            }

            pc++;

            if(pc == Ippl::getNodes())//nodes number should be equal in BoundaryGeometry and PriemissionPhysics.
                pc = 0;

        }
    }

    if(count) {
        delete count;
    }
}
