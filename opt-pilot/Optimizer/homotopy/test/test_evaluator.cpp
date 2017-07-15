//Validate the evaluator

#include "gtest/gtest.h"

#include <vector>
#include "boost/function.hpp"

#include "HomotopyTypes.h"
#include "Problems.h"
#include "evaluator.hpp"

#include <stdio.h>

using namespace Homotopy;

namespace {

	TEST(EVALUATOR, Local){
		//Set up the problem
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the evaluator framework
        Evaluator<EvaluationStrategy::Local< functionSet_t > > E( P->Objectives );

        //independent variables
        std::vector< double > v(3, 0.5);

        //get the raw objectives
        std::vector< double > obj;
        obj.push_back( (P->Objectives[0])( v ) );
        obj.push_back( (P->Objectives[1])( v ) );

        //expect invariance
        EXPECT_EQ( E.eval( v ), obj );
	}

    TEST(EVALUATOR, CachedLocal){
		//Set up the problem
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the evaluator framework
        Evaluator< EvaluationStrategy::Cached < EvaluationStrategy::Local < functionSet_t > > > E( P->Objectives );

        //independent variables
        std::vector< double > v(3, 0.5);

        //get the raw objectives
        std::vector< double > obj1;
        obj1.push_back( (P->Objectives[0])( v ) );
        obj1.push_back( (P->Objectives[1])( v ) );

        //expect invariance
        EXPECT_EQ( E.eval( v ), obj1 );

        //new independent variables
        std::vector< double > w(3, 0.9);
        //get the raw objectives
        std::vector< double > obj2;
        obj2.push_back( (P->Objectives[0])( w ) );
        obj2.push_back( (P->Objectives[1])( w ) );

        //eval at another point
        EXPECT_EQ( E.eval( w ), obj2 );

        //old independent variables
        std::vector< double > x(3, 0.5);
        //expect a cached result
        EXPECT_EQ( E.eval( x ), obj1 );
	}
    
    /*
    //  Decided not to support this until I can think of a better way to
    //  distribute thread-capable code from source.  

    TEST(EVALUATOR, Async){
		//Set up the problem
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the evaluator framework
        Evaluator<EvaluationStrategy::AsyncLocal> E( P );

        //independent variables
        std::vector< double > v(3, 0.5);

        //get the raw objectives
        std::vector< double > obj;
        obj.push_back( (P->Objectives[0])( v ) );
        obj.push_back( (P->Objectives[1])( v ) );

        //expect invariance
        EXPECT_EQ( E.eval( v ), obj );
	}
    */
}