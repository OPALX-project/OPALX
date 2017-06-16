#ifndef Scalarization_hpp
#define Scalarization_hpp

/*
Scalarizing adapter allowing dynamic weight changes
*/
#include "boost/bind.hpp"
#include <cassert>
#include <stdio.h>
#include <numeric>

#include "HomotopyTypes.h"

class ScalarizationInterface {
public:
    Homotopy::function_t f;
	Homotopy::functionSet_t EqualityConstraints;
	Homotopy::functionSet_t InequalityConstraints;

    int dimDesign;
	int dimObj;

    std::vector< double > lowerBounds;
	std::vector< double > upperBounds;

    virtual double operator() (const std::vector<double> &x) =0;
    virtual double operator() (const std::vector<double> &x, bool &valid) =0;

    ScalarizationInterface(Problem::Interface * P) :    f(NULL), EqualityConstraints(P->EqualityConstraints),
                                                        InequalityConstraints(P->InequalityConstraints),
                                                        lowerBounds(P->lowerBounds), upperBounds(P->upperBounds)
	{
        this->dimObj = 1;
		this->dimDesign = P->dimDesign;
    }
};

//The public interface is quite similar to the Problem::Interface
//and could probably be refactored.
template< typename T >
class Scalarization : public ScalarizationInterface {
protected:
	virtual double eval(const std::vector<double> &x) =0;
    virtual double eval(const std::vector<double> &x, bool &valid) =0;
	Problem::Interface * P;
public:
    T e;
    Scalarization(Problem::Interface * p, typename T::BaseType e) : ScalarizationInterface(p),
                                            P(p), e(e)
    {        
		this->dimObj = 1;
		this->dimDesign = P->dimDesign;
	}

	//operator () overload
	double operator() (const std::vector<double> &x){
		return (this->eval(x));
	}
    double operator() (const std::vector<double> &x, bool &valid){
        return (this->eval(x, valid));
    }
};

template< typename T >
class FixedScalarization : public Scalarization< T > {
protected:
	std::vector<double> * weights;

	//Could probably get rid of this when 
	//the interface is stabilized
	std::vector<double> default_weights;

    double eval(const std::vector<double> &at){	
        bool valid;
        return this->eval( at, valid );
	}

	double eval(const std::vector<double> &at, bool &valid){	
        Homotopy::objVars_t result( this->e.eval( at, valid ) );
        if(!valid) return 0.0;
        else return std::inner_product( this->weights->begin(), this->weights->end(), result.begin(), 0.0 );
	}

public:
	FixedScalarization(Problem::Interface * p, typename T::BaseType e) : Scalarization<T>(p, e), weights(NULL) {
		int i;
		for(i=0; i<this->P->dimObj; i++) { this->default_weights.push_back(1.0); }
		this->weights = &(this->default_weights);
	}
	
	//Careful: No validation of pointer
	void SetWeights(std::vector<double> * w) { 
		assert ( w->size() == this->P->dimObj );
		this->weights = w; 
	}
};


template< typename T >
class DynamicScalarization : public Scalarization < T > {
protected:
    double eval(const std::vector<double> &at){	
        bool valid;
        return this->eval( at, valid );
	}
    double eval(const std::vector<double> &at, bool &valid){
        //Allow for the weights to be part of the optimization
		//by assuming theyre tacked on to the end of the at vector
		//NOTE: the last weight is the one implied by the 1-sum(others)

		//Determine where the weights start
		int offset = this->P->dimDesign;

		//Separate the set of "Problem" variables
		const Homotopy::designVars_t vars(at.begin(), at.begin() + offset);
        Homotopy::designVars_t weights( at.begin() + offset, at.end() );

        //Evaluate the objectives with the designVars
        Homotopy::objVars_t result( this->e.eval( vars, valid ) );
    
        //Return immediately if result is unavailable
        if(!valid) return 0.0;
        else {
            //Calculate the remaining weight and add it
            weights.push_back( 1.0 - std::accumulate(weights.begin(), weights.end(), 0.0) );
            //Return the weighted sum
            return std::inner_product( weights.begin(), weights.end(), result.begin(), 0.0 );
        }
	}

public:
	DynamicScalarization(Problem::Interface * p, typename T::BaseType e) : Scalarization<T>(p, e) {
		this->dimDesign = this->P->dimDesign + this->P->dimObj - 1;

		//Add constraints on the individual lambda values
		int i;
		for(i=0; i<this->P->dimObj - 1; i++) {
			this->lowerBounds.push_back(0.0);
			this->upperBounds.push_back(1.0);
		}
	}
};

#endif