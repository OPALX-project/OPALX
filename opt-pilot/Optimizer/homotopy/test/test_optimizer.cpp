//Simple test to verify Optimizer class functionality

#include "gtest/gtest.h"

#include <stdio.h>

#include "homotopy.hpp"

using namespace Homotopy;

namespace {

    // Tests that the optimizer works
    TEST(OPTIMIZERTest, Alive) {
	
	    int DesignVars = 3;

        //Choose a Problem to optimize
	    //Problem::Interface * P = Problem::Factory("BASIN", 1, DesignVars);
	    //Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);
	    //Problem::Interface * P = Problem::Factory("CONST_TEST", 1, 2);
        Problem::Interface * P = Problem::Factory("SURROGATE", 2, 2);

        //Scalarize the problem
        FixedScalarization< Evaluator<EvaluationStrategy::Local< functionSet_t > > > S(P, P->Objectives);

        //Establish parameters for finite differences
        FiniteDifferences::Params_t FDpar = { 1e-8, FiniteDifferences::CENTRAL };

        //Create optimizer
	    optimizer * op = new OptNlopt(&S, 1e-7, FDpar);
	
	    printf("starting at: ");
	    std::vector<double> x(S.dimDesign);
	    int i;
	    for(i=0; i<x.size(); i++) { 
		    x[i] = (P->upperBounds[i] - P->lowerBounds[i])/2.0 + P->lowerBounds[i];
		    printf("%lf ", x[i]);
	    }
	    printf("\n");

        //Establish a weighting scheme
	    std::vector<double> w(P->Objectives.size(), 0.0);
        w[1] = 1.0;
	    S.SetWeights(&w);

        //Run the optimization
	    op->RunFrom(x);

        //Get & print results
        double result = S(x);
	    printf("minimum at: ");
	    for(i=0; i<x.size(); i++) { printf("%lf ", x[i]); }
	    printf(" =  %lf\n", result);
    }

}