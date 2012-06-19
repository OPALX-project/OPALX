// include files
#include "Utility/PoomaTimings.h"
#include <iostream>
#include <Const.hh>

extern "C" {
#include "treecode.h"
}

using namespace std;

////////////////////////////////////////////////
//constructor
template <class T, unsigned int Dim>
TreePoissonSolver<T,Dim>::TreePoissonSolver(ChargedParticles<T,Dim> *beam, T t, T e, int r, bool u):
    theta_m(t),
    eps_m(e),
    rsize_m(r),
    usequad_m(u)
{

    setType(TreeSolver);

    //Initial values for BH solver
    nbody = beam->getTotalNum();
    theta = theta_m;
    eps = eps_m;

    //Defaults
    usequad = usequad_m;
    rsize = rsize_m;

    INFOMSG("Parameters for BH Code are: \nTheta= " << theta << " Eps= " << eps
            << endl << "Usequad= " << usequad << " Rsize= " << rsize << endl
           );


    //Allocate memory for bodies
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
}


////////////////////////////////////////////////////////////////////////////
// destructor
template <class T, unsigned int Dim>
TreePoissonSolver<T,Dim>::~TreePoissonSolver() {}


template<class T, unsigned int Dim>
void TreePoissonSolver<T,Dim>::computeElectricField(ChargedParticles<T,Dim> *beam, Vektor<T,Dim> hr)
{

    //Fill ACTUAL body data into BH body array, at each time, Electric Field is to be calculated!

    int counter=0;
    for (bodyptr p = bodytab; p < bodytab+nbody; p++) {
        Type(p) = BODY;
        Mass(p) = beam->qp_m[counter];
        Pos(p)[0] = beam->R[counter][0];
        Pos(p)[1] = beam->R[counter][1];
        Pos(p)[2] = beam->R[counter][2];
        Vel(p)[0] = beam->P[counter][0];
        Vel(p)[1] = beam->P[counter][1];
        Vel(p)[2] = beam->P[counter][2];
        counter++;
    }

    //Calculates BH Acceleration
    treeforce();

    // Convert BH Acceleration to Electric Field
    counter=0;

    for (bodyptr p = bodytab; p < bodytab+nbody; p++) {

        beam->Ef[counter][0] =  - beam->getCoupling() * Acc(p)[0];
        beam->Ef[counter][1] =  - beam->getCoupling() * Acc(p)[1];
        beam->Ef[counter][2] =  - beam->getCoupling() * Acc(p)[2];


        /*
           cout << "Pos(Pooma)= " << beam->R[counter][1] <<endl;
           cout << "Vel(Pooma)= " << beam->P[counter][1] <<endl;
           cout << "Ef(Pooma)= " << beam->Ef[counter][1] <<endl;
           */

        counter++;
    }

    // INFOMSG("Sum Ef= " << sum(beam->Ef) << endl);
}

