//Validate the evaluator

#include "gtest/gtest.h"

#include "boost/bind.hpp"
#include "boost/function.hpp"

#include "Problems.h"
#include "HomotopyTypes.h"

#include "JobQueue.hpp"
#include "CommInterface.hpp"

#include <stdio.h>

using namespace Homotopy;

namespace {

    void Print(const std::vector< double > & x ){
        int i;
	    for(i=0; i<x.size(); i++) { 
		    printf("%lf ", x[i]);
	    }
	    printf("\n");
    }

	TEST(JOBQUEUE, Alive){
        //Set up the problem
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the evaluator framework
        Communication::SimulatedRemote< functionSet_t > C( P->Objectives );
        JobQueue< Communication::SimulatedRemote< functionSet_t > > E( C );

        //Select some variables
        std::vector< double > v(3,0.5);

        //expect empty (still being computed)
        std::vector<double> empty;
        EXPECT_EQ( E.eval( v ), empty );

        E.Poll();

        //expect results
        std::vector<double> val;
        val.push_back( (P->Objectives[0])(v) );
        val.push_back( (P->Objectives[1])(v) );
        EXPECT_EQ( E.eval( v ), val );
	}
    
    TEST(JOBQUEUE, Multi){
        //Set up the problem
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the evaluator framework
        Communication::SimulatedRemote< functionSet_t > C( P->Objectives );
        JobQueue< Communication::SimulatedRemote< functionSet_t > > E( C );

        //Select sets of variables
        std::vector< std::vector< double > > v;
        int NUM=100;
        int i;
        for(i=0; i<NUM; i++){
            std::vector< double > tmp(3,(double)i/(double)NUM);
            v.push_back( tmp );
        }

        //expect empty (still being computed)
        std::vector<double> empty;
        for(i=0; i<NUM; i++){
            EXPECT_EQ( E.eval( v[i] ), empty );
        }
        
        E.Poll();

        //expect results
        for(i=0; i<NUM; i++){
            std::vector<double> val;
            val.push_back( (P->Objectives[0])(v[i]) );
            val.push_back( (P->Objectives[1])(v[i]) );
            EXPECT_EQ( E.eval( v[i] ), val );
        }
        
	}

    TEST(JOBQUEUE, Groups){
        //Set up the problem
        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        //Set up the evaluator framework
        Communication::SimulatedRemote< functionSet_t > C( P->Objectives );
        JobQueue< Communication::SimulatedRemote< functionSet_t > > E( C );

        //How many evals to request
        int NUM=57;

        //Select sets of variables
        std::vector< std::vector< double > > v;
        int i;
        for(i=0; i<2*NUM; i++){
            std::vector< double > tmp(3,(double)i/(double)(2*NUM));
            v.push_back( tmp );
        }

        //Set a group
        E.NewGroup(5);

        //expect empty (still being computed)
        std::vector<double> empty;
        for(i=0; i<NUM; i++){
            EXPECT_EQ( E.eval( v[i] ), empty );
        }

        //Set a new group
        E.NewGroup(1);

        //expect empty (still being computed)
        for(i=NUM; i<(2*NUM)-1; i++){
            EXPECT_EQ( E.eval( v[i] ), empty );
        }
        
        int ind = E.Poll();
        EXPECT_EQ( ind, 5 );

        //expect results
        for(i=0; i<NUM; i++){
            std::vector<double> val;
            val.push_back( (P->Objectives[0])(v[i]) );
            val.push_back( (P->Objectives[1])(v[i]) );
            EXPECT_EQ( E.eval( v[i] ), val );
        }
	}
}