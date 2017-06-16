//Basic marching homotopy test to validate additional
//problem constraints

#include "gtest/gtest.h"

#include "homotopy.hpp"

namespace {

	void Print(const std::vector< double > &x){
		int i;
		for(i=0; i<x.size(); i++) {printf("%lf ", x[i]);}
	}
	std::vector< double > GetF(const std::vector< Problem::FUNCTION > &f, const std::vector< double > &x){
		std::vector< double > result;
		int i;
		for(i=0; i<f.size(); i++) { result.push_back((f[i])(x)); }
		return result;
	}
	double FDist(const std::vector< Problem::FUNCTION > &f, double step, const std::vector< double > &y, const std::vector< double > &x){
		std::vector< double > fy, fx;
		fy = GetF(f, y);
		fx = GetF(f, x);
		int i;
		double dist = 0.0;
		for(i=0; i<f.size(); i++) {
			dist += (fy[i] - fx[i]) * (fy[i] - fx[i]);
		}
		return fabs( sqrt(dist) - step );
	}
	double Dist(double step, const std::vector< double > &y, const std::vector< double > &x){
		int i;
		double dist = 0.0;
		for(i=0; i<y.size(); i++) {
			dist += (y[i] - x[i]) * (y[i] - x[i]);
		}
		return fabs( sqrt(dist) - step );
	}
	void PrintF(std::vector< double > &x, std::vector< Problem::FUNCTION > &f){ 
		Print(x);
		printf(" -> ");
		Print(GetF(f, x));	
		printf("\n");
	}

    typedef Evaluator< EvaluationStrategy::Cached< EvaluationStrategy::Local< functionSet_t > > > eval_t;

	TEST(MARCHINGHOMOTOPY, 2D){

		int Points = 146;
		double step = 0.01;
		

		//Set up the problem
		Problem::Interface * P = Problem::Factory("FON", 2, 3);

		//Points = 20;
		//step = 0.1;
		//Problem::Interface * P = Problem::Factory("WFG2", 2, 5);

        //Scalarize the problem
        FixedScalarization< Evaluator<EvaluationStrategy::Local< functionSet_t > > > S(P, P->Objectives);
        std::vector< double > w(2, 0.0); w[1] = 1.0;
        S.SetWeights(&w);

        //Finite difference params
        FiniteDifferences::Params_t FDpar = { 1e-6, FiniteDifferences::FORWARD };
		
		//Get the starting point
		optimizer * op = new OptNlopt(&S, 1e-4, FDpar);
		std::vector<double> x1(P->dimDesign, 0.5);
		op->RunFrom(x1);
		PrintF(x1, P->Objectives);
        
        //Re-scalarize the problem
        DynamicScalarization< eval_t > D(P, P->Objectives);
		
        //Establish auxilary stepping constraints
        //StepConstraint< typename Problem::FUNCTION > C(NULL, step);
        boost::function<objVars_t (const designVars_t&)> f = boost::bind( &eval_t::eval, &(D.e), _1);
        FStepConstraint< boost::function<objVars_t (const designVars_t&)> > C(f, P->dimDesign, NULL, step);
		D.EqualityConstraints.push_back( C.function );

		op = new OptNlopt(&D, 1e-4, FDpar);
		x1.push_back(0.5);
		std::vector<double> x(x1);

		//Get the rest of the points
		int i;
		for(i=0; i<Points; i++){
			std::vector< double > last_x(x.begin(), x.begin() + P->dimDesign);
			std::vector< double > last_f( GetF(P->Objectives, last_x) );

			//Build the additional constraint
			C.UpdateFrom( &last_f );
			//C.UpdateFrom( &last_x );

			op->RunFrom(x);
			
			PrintF(x, P->Objectives);

			//printf("DistErr: %lf\n", FDist(P->Objectives, 0, last_x, x));
			//printf("DistErr: %lf\n", Dist(0.0, last_x, x));
		}
	}
}