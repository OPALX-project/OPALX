/*

Test program showing the general state of functionality

*/

#include <stdio.h>
#include <vector>

#include "boost/function.hpp"

#include "homotopy.hpp"

using namespace Homotopy;

int main(int argc, char** argv){

	//Choose a Problem
	Problem::Interface * P = Problem::Factory("FON", 2, 3);
	//Problem::Interface * P = Problem::Factory("WFG5", 3, 30);
	//Problem::Interface * P = Problem::Factory("DTLZ2", 3, 30);
    //Problem::Interface * P = Problem::Factory("SURROGATE", 2, 2);

    //Set up the communication framework        
    Homotopy::Communication::AdHoc< Homotopy::Communication::CommImpl::Simulator< functionSet_t > > Comm;
    Comm.comm_.P = P->Objectives;
    Comm.Dispatcher = boost::bind( &Homotopy::Communication::CommImpl::Simulator< functionSet_t >::dispatcher,
        &Comm.comm_, _1, _2 );
    Comm.Handler = boost::bind( &Homotopy::Communication::CommImpl::Simulator< functionSet_t >::handler,
        &Comm.comm_, _1, _2, _3, _4 );
    Comm.Exchanger = boost::bind( &Homotopy::Communication::CommImpl::Simulator< functionSet_t >::exchange,
        &Comm.comm_, _1 );

	//Instantiate the homotopy
	Pareto::homotopy h( P, 1e-3, 1e-6, Comm );

	//Homotopically deform the ansatz
	h.GetFront(19, 1, 0, 1);

	return 0;
}