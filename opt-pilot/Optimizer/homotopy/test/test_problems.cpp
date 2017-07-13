//Simple test to verify Problem class functionality

#include "gtest/gtest.h"

#include <math.h>
#include <cstdio>

#include "homotopy.hpp"

using namespace Homotopy;

namespace {

    TEST(PROBLEMSTest, SingleObjective) {

	    int DesignVars = 3;

	    Problem::Interface * P = Problem::Factory("BASIN", 1, DesignVars);

	    nlopt::opt opt(nlopt::LD_SLSQP, DesignVars);

	    opt.set_lower_bounds(P->lowerBounds);
	    opt.set_upper_bounds(P->upperBounds);

        std::vector <typename Problem::FUNCTION> empty;
        bool valid= true;
        FiniteDifferences::Params_t FDpar = { .000001, FiniteDifferences::CENTRAL };

	    NloptAdapt< typename Problem::FUNCTION > NA(P->Objectives[0], &empty, &empty, valid, FDpar );

	    opt.set_min_objective(&NloptAdapt<typename Problem::FUNCTION>::ObjIface, (void*)&NA);

	    opt.set_xtol_rel(1e-4);

	    std::vector<double> x(DesignVars);
	    int i;
	    for(i=0; i<x.size(); i++) { x[i] = 0.3; }
	    double minf;
	    nlopt::result result = opt.optimize(x, minf);
	    printf("found minimum at f(%lf,%lf) = %lf\n", x[0], x[1], minf);
    }

    TEST(PROBLEMSTest, MultiObjective) {

	    int DesignVars = 3;

	    Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);

	    nlopt::opt opt(nlopt::LD_SLSQP, DesignVars);
	    opt.set_lower_bounds(P->lowerBounds);
	    opt.set_upper_bounds(P->upperBounds);

	    //Scalarize the problem
        FixedScalarization< Evaluator<EvaluationStrategy::Local< functionSet_t > > > S(P, P->Objectives);
	    std::vector<double> weights;
	    weights.push_back(1.0);
	    weights.push_back(1.0);
	    S.SetWeights(&weights);

	    std::vector <typename Problem::FUNCTION> empty;
        bool valid= true;
        FiniteDifferences::Params_t FDpar = { .000001, FiniteDifferences::CENTRAL };

	    NloptAdapt< typename Problem::FUNCTION > NA(S.f, &empty, &empty, valid, FDpar );

	    opt.set_min_objective(&NloptAdapt< typename Problem::FUNCTION >::ObjIface, (void*)&NA);

	    opt.set_xtol_rel(1e-4);

	    std::vector<double> x(S.dimDesign);
	    int i;
	    for(i=0; i<x.size(); i++) { x[i] = 0.3; }
	    printf("starting at (%lf,%lf,%lf) \n", x[0], x[1], x[2]);
	    double minf=0;
	    nlopt::result result = opt.optimize(x, minf);
	    printf("found minimum at f(%lf,%lf,%lf) = %lf\n", x[0], x[1], x[2], minf);

	
	}
}  //namespace
