// include files
#include "Utility/PoomaTimings.h"
#include <iostream>
#include <Const.hh>
using namespace std;

////////////////////////////////////////////////
//constructor
template <class T, unsigned int Dim>
P2PPoissonSolver<T,Dim>::P2PPoissonSolver(T e):eps_m(e) {
    INFOMSG("Eps= " << eps_m << endl);
    setType(P2PSolver);
}


////////////////////////////////////////////////////////////////////////////
// destructor
template <class T, unsigned int Dim>
P2PPoissonSolver<T,Dim>::~P2PPoissonSolver() {}

///////////////////////////////////////////
// calculate electric field via PP summation of Field Energy.

template<class T, unsigned int Dim>
void P2PPoissonSolver<T,Dim>::computeElectricField(ChargedParticles<T,Dim> *beam,
        Vektor<T,Dim> hr)
{

    //beam->Epot_m = 0.0;
    beam->Ef = 0.0;

    //calculating electric field WITH SMOOTHING EPSILON


    for (unsigned long i=0; i<beam->getTotalNum(); i++) {
        for (unsigned long j=0; j<beam->getTotalNum(); j++) {
            if (i!=j) beam->Ef[i] +=  beam->getCoupling()*beam->qp_m[j] * ( beam->R[i] - beam->R[j] )
                / pow( (dot(beam->R[i] - beam->R[j] , beam->R[i] - beam->R[j]) + eps_m*eps_m),1.5);
        }
    }
}



