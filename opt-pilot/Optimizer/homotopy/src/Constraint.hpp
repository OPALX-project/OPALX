/*
General constraint interface
*/

#ifndef Constraint_hpp
#define Constraint_hpp

#include "HomotopyTypes.h"

#include <cassert>

class GeneralConstraint {
protected:
	virtual double eval(const std::vector< double > &x) =0;
public:
	Homotopy::function_t function;
	GeneralConstraint() :
		function( boost::bind( &GeneralConstraint::eval, this, _1 ) ) {}
};

#include <numeric>

class BoundSumConstraint : public GeneralConstraint {
protected:
	int sign;
	double value;
	int from;
	int to;

	double eval(const std::vector< double > &x) {
		return this->sign * ( this->value - std::accumulate(x.begin() + this->from, x.begin() + this->to, 0.0) );
	}

public:
	enum SIGN { LESS_THAN = -1, GREATER_THAN = 1 };
	BoundSumConstraint( SIGN s, double value, int from, int to ) : 
	GeneralConstraint (), sign(s), value(value), from(from), to(to) {}
};

#include <math.h>

class StepConstraint : public GeneralConstraint  {
protected:
	const std::vector< double > * x1;
	double step;
	static double L2Dist(const std::vector< double > &x1, const std::vector< double > &x2 ) {
		double sum = 0.0;
		for(unsigned i=0;i<x1.size();i++){
			sum += pow( x1[i] - x2[i], 2.0 );
		}
		//return sqrt(sum);
        return sum;
	}
	virtual double eval(const std::vector< double > &x){
		return L2Dist(*this->x1, x) - this->step;
	}
public:
	StepConstraint( const std::vector< double >  * x1, double step ) :
	  GeneralConstraint (), x1(x1), step(step) { }
	void UpdateFrom( const std::vector< double > * x1) { this->x1 = x1; }
};
template< typename T >
class FStepConstraint : public StepConstraint  {
protected:
	const T & F_trans;
    const int designVars;
	std::vector< double > f;
	void Feval(const std::vector< double > &x){
        //Only evaluate a subset of the vector or else 
        std::vector< double > vars( x.begin(), x.begin() + this->designVars );
		this->f = ( (F_trans)(vars) );
	}
	virtual double eval(const std::vector< double > &x){
		this->Feval(x);
		return this->L2Dist(*this->x1, this->f) - this->step ;
	}
public:
	FStepConstraint( const T &F_trans, int dV, const std::vector< double >  * f1, double step) :
	  F_trans(F_trans), designVars(dV), f(0), StepConstraint( f1, step) {}
};

template< typename T >
class FEqDistanceConstraint : public FStepConstraint< T >  {
protected:
	std::vector< double > const * x2;

	double eval(const std::vector< double > &x){
		//Assign NULL to deactivate constraints (always satisfied)
		if( !(this->x1 && this->x2) ) return 0;

		this->Feval(x);
		return (this->L2Dist(*this->x1, this->f) - this->L2Dist(*this->x2, this->f) ) ;
	}
public:
	FEqDistanceConstraint( const T &F_trans, int dV, const std::vector< double >  * f1, const std::vector< double >  * f2) :
	  x2(f2), FStepConstraint< T >(F_trans, dV, f1, 0.0) {}
	void UpdateFrom( std::vector< double > * const f1, std::vector< double > * const f2) { 
		this->x1 = f1; 
		this->x2 = f2; 
	}
};

#endif
