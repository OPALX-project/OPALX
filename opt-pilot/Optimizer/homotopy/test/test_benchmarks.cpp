//Basic marching homotopy test to validate additional
//problem constraints

#include "gtest/gtest.h"
#include "homotopy.hpp"

#include <algorithm>

namespace {

    //Problem::Interface * P = Problem::Factory("SURROGATE", 3, 2);

    TEST(HOMOTOPY, DISABLED_FON){
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(15, 1, 0, 1);
	}

	TEST(HOMOTOPY, DISABLED_MOTTA1){
        Problem::Interface * P = Problem::Factory("MOTTA1", 3, 3);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(15, 5, 0, 1);
	}
    TEST(HOMOTOPY, DISABLED_MOTTA2){
        Problem::Interface * P = Problem::Factory("MOTTA2", 3, 3);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(15, 1, 0, 1);
	}
    TEST(HOMOTOPY, DISABLED_MOTTA3){
        Problem::Interface * P = Problem::Factory("MOTTA3", 3, 2);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::CENTRAL;

		//Homotopically deform the ansatz
		h.GetFront(15, 1, 0, 1);
	}
    TEST(HOMOTOPY, DISABLED_MOTTA4){
        Problem::Interface * P = Problem::Factory("MOTTA4", 4, 4);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(10, 5, 0, 1);
	}

    TEST(HOMOTOPY, DISABLED_WFG5){
        Problem::Interface * P = Problem::Factory("WFG5", 3, 3);
        /*
        std::vector< std::vector< double > > Design;
        std::vector< std::vector< double > > Objective;
        std::vector< std::vector< double > > Lambda;

        std::vector< double > d(10, 0.5);
        d[0] = 0.001; d[1] = 0.001165 ; d[2] = 0.001;
        Design.push_back( d );
        d[0] = 0.0; d[1] = 1.0 ; d[2] = 0.5;
        Design.push_back( d );
        d[0] = 0.0; d[1] = 0.0 ; d[2] = 0.5;
        Design.push_back( d );

        std::vector< double > o(3, 0.0);
        o[0] = 0.0633; o[1] = 0.37178 ; o[2] = 6.031;
        Objective.push_back( o );
        o[0] = 0.0; o[1] = 1.0 ; o[2] = 0.0;
        Objective.push_back( o );
        o[0] = 1.0; o[1] = 0.0 ; o[2] = 0.0;
        Objective.push_back( o );

        o[0] = 1.0; o[1] = 1.0 ; o[2] = 0.0;
        Lambda.push_back( o );
        o[0] = 1.0; o[1] = 0.0 ; o[2] = 1.0;
        Lambda.push_back( o );
        o[0] = 0.0; o[1] = 1.0 ; o[2] = 1.0;
        Lambda.push_back( o );
        */
        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-6, 1e-6, Comm );
		//Pareto::homotopy h( P, 1e-3, 1e-6, Comm, Design, Objective, Lambda );
        //h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(15, 5, 0, 1);
	}

    TEST(HOMOTOPY, DISABLED_SURROGATE){
        Problem::Interface * P = Problem::Factory("SURROGATE", 3, 2);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(15, 5, 0, 1);
	}

    TEST(HOMOTOPY, DTLZ2){
        Problem::Interface * P = Problem::Factory("DTLZ2", 3, 10);

        std::vector< std::vector< double > > Design;
        std::vector< std::vector< double > > Objective;
        std::vector< std::vector< double > > Lambda;

        std::vector< double > d(10, 0.5);
        d[0] = 1.0; d[1] = 1.0 ; d[2] = 0.5;
        Design.push_back( d );
        d[0] = 0.0; d[1] = 1.0 ; d[2] = 0.5;
        Design.push_back( d );
        d[0] = 0.0; d[1] = 0.0 ; d[2] = 0.5;
        Design.push_back( d );

        std::vector< double > o(3, 0.0);
        o[0] = 0.0; o[1] = 0.0 ; o[2] = 1.0;
        Objective.push_back( o );
        o[0] = 0.0; o[1] = 1.0 ; o[2] = 0.0;
        Objective.push_back( o );
        o[0] = 1.0; o[1] = 0.0 ; o[2] = 0.0;
        Objective.push_back( o );

        o[0] = 1.0; o[1] = 1.0 ; o[2] = 0.0;
        Lambda.push_back( o );
        o[0] = 1.0; o[1] = 0.0 ; o[2] = 1.0;
        Lambda.push_back( o );
        o[0] = 0.0; o[1] = 1.0 ; o[2] = 1.0;
        Lambda.push_back( o );

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
		Pareto::homotopy h( P, 1e-3, 1e-6, Comm, Design, Objective, Lambda );
        //Pareto::homotopy h( P, 1e-3, 1e-6, Comm );
        h.fd_type = FiniteDifferences::FORWARD;

		//Homotopically deform the ansatz
		h.GetFront(15, 5, 0, 1);
	}
}