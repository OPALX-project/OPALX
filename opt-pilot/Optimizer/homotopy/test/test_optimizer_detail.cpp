//Permits detailed optimizer testing with specific cases
//as well as facilities to determine basins of attraction

#include "gtest/gtest.h"

#include "homotopy.hpp"

#include <math.h>

typedef Evaluator< EvaluationStrategy::Local< functionSet_t > > eval_t;

namespace {

    double Dist( std::vector< double > l, std::vector< double > r ){
        double sum = 0.0;
        for( int i=0; i<l.size(); i++){
            sum += (l[i]-r[i])*(l[i]-r[i]);
        }
        return sqrt( sum );
    }
    void Print( std::vector< double > v ){
        for( int i=0; i<v.size(); i++) { printf("%lf ", v[i]); }
        printf("\n");
    }
    void Print2( std::vector< double > v ){
        for( int i=0; i<v.size(); i++) { printf("%lf ", v[i]); }
        
    }

	TEST(OPTIMIZERDETAIL, DISABLED_FON){
		//Set up the problem
		Problem::Interface * P = Problem::Factory("FON", 2, 3);
		//Problem::Interface * P = Problem::Factory("WFG5", 3, 30);
		//Problem::Interface * P = Problem::Factory("DTLZ2", 3, 30);
        //Problem::Interface * P = Problem::Factory("SURROGATE", 2, 2);

        //std::vector< double > v(3, 0.6);
        //P->upperBounds = v;

        //double tolerance = 1e-8;
        double tolerance = 1e-6;
        double fd_step = 1e-8;

        //Establish finite difference parameters
        FiniteDifferences::Params_t FDpar = { fd_step, FiniteDifferences::FORWARD };

        DynamicScalarization< eval_t > Scal( P, P->Objectives );

		//Sum constraints on the Lambda parameters
		BoundSumConstraint LT( BoundSumConstraint::LESS_THAN, 1.0, P->dimDesign, Scal.dimDesign);
		BoundSumConstraint GT( BoundSumConstraint::GREATER_THAN, 0.0, P->dimDesign, Scal.dimDesign);
		Scal.InequalityConstraints.push_back( LT.function );
		Scal.InequalityConstraints.push_back( GT.function );

		//Equidistance constraints
        boost::function<objVars_t (const designVars_t&)> f = boost::bind( &eval_t::eval, &(Scal.e), _1);
		FEqDistanceConstraint< boost::function <objVars_t (const designVars_t&)> > eqd(f, P->dimDesign, NULL, NULL) ;

        //Build the neighbors
        //Nearest neighbor
        double l[2] = { 0.020034786412537 , 0.961522724540055 };
        //Next-nearest jneighbor
        //double l[2] = { 0.000003128154552 , 0.981554382798084  };
        
        std::vector< double > left;
        left.assign( l, l+2 );
        //nearest
        double r[2] = { 0.060098102928507 , 0.921459408023997 };
        //next-nearest
        //double r[2] = { 0.080129761186492 , 0.901427749765968 };
        std::vector< double > right;
        right.assign( r, r+2 );

        //Set the neighbors
        eqd.UpdateFrom( &left, &right );
	    Scal.EqualityConstraints.push_back( eqd.function );

        
        /*
        double x[3] = { 0.529281856518440 , 0.529281856518440 , 0.529281856534945 };
        //double x[3] = { 0.4 , 0.4 , 0.4 };
        std::vector< double > start;
        start.assign( x, x+3 );
        //Add the lambda value
        start.push_back( 0.959184 );


        opt.RunFrom( start );
        */

        double x = 0.100;
        double step = 1e-3;

        printf("x: y:\n");
        while (x < 1.0){
            OptNlopt opt( &Scal, tolerance, FDpar);
            
            std::vector< double > v(3, x);
            v.push_back( 0.959184 );
            opt.RunFrom( v );
            printf("%lf %lf %lf\n", x, v[0], (P->Objectives[0])( v ) );
            x += step;
        }

        /*
        printf("Result:\n");
        Print( start );
        printf("Location:\n");
        Print( Scal.e.eval( start ) );

        //Strip the lambda value
        start.pop_back();
        printf("Distances: %lf %lf\n", Dist( left, Scal.e.eval( start ) ) , Dist( right, Scal.e.eval( start ) ) );
        //Keeps going to 0.758593590638663 0.758593591566017 0.758593591566104
        //Should go to near 0.485797761529296 0.485797761524148 0.485797761531280 -> 0.024832068107697 0.966320003908558 
        */
	}

    TEST(OPTIMIZERDETAIL, DISABLED_FONprobe){

        Problem::Interface * P = Problem::Factory("FON", 2, 3);

        double x = 0.0;
        double step = 1e-3;

        printf("x: y:\n");
        while (x < 1.0){
            std::vector< double > v(3, x);
            printf("%lf %lf\n", x, (P->Objectives[0])( v ) );
            x += step;
        }
    }

    TEST(OPTIMIZERDETAIL, WFG5probe){

        Problem::Interface * P = Problem::Factory("WFG5", 3, 3);

        //Establish finite difference parameters
        FiniteDifferences::Params_t FDpar = { 1e-6, FiniteDifferences::FORWARD };
        FixedScalarization< eval_t > Scal( P, P->Objectives );

        std::vector< double > w (3, 0.0);
        w[2] = 1.0;
        Scal.SetWeights( &w );

        double x = 0.001;
        double step = 5e-2;

        while (x < 1.0){
            std::vector< double > v(3, x);
            Print2(v);
            
            //OptNlopt opt( &Scal, 1e-3, FDpar, nlopt::GN_ISRES );
            OptNlopt opt( &Scal, 1e-6, FDpar, nlopt::LN_COBYLA );
            opt.RunFrom(v);

            std::vector<double> f;
            f.push_back( (P->Objectives[0])( v ) );
            f.push_back( (P->Objectives[1])( v ) );
            f.push_back( (P->Objectives[2])( v ) );

            printf(" -> ");
            Print2(f);
            printf("\n");
            x += step;
        }
    }
}