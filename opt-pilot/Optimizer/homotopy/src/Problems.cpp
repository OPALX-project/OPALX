#include <vector>

#include "boost/bind.hpp"
#include "boost/function.hpp"

#include <string>
#include <cassert>
#include <stdio.h>

#include <iostream>

#include "Problems.h"

#include <numeric>
namespace helpers {
	double addsquare( double l, double r )	{ return l + (r) * (r); }
	double square( double x )	{ return (x) * (x); }
}

//Trivial quadratic problem
class Basin : public Problem::Interface{
private:
	double obj(const std::vector< double > &x){
		return std::accumulate(x.begin(), x.end(), 0.0, helpers::addsquare);
	}
public:
	Basin(int DimObj, int DimDesign)	{

		typename Problem::FUNCTION f( boost::bind(&Basin::obj, this, _1) );
		this->Objectives.push_back(	f );

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( 0.0 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};

#include <math.h>

//Trivial constrainted problem (from nlopt tutorial)
class ConstrainedSR : public Problem::Interface{
private:
	int a1, a2, b1, b2;

	double obj(const std::vector< double > &x){
		return sqrt(x[1]);
	}
	double constr1(const std::vector< double > &x){
		return -1 * x[1];
	}
	double constr2(const std::vector< double > &x){
		return pow(a1*x[0] + b1, 3) - x[1];
	}
	double constr3(const std::vector< double > &x){
		return pow(a2*x[0] + b2, 3) - x[1];
	}
public:
	ConstrainedSR(int DimObj, int DimDesign){

		a1 = 2;
		a2 = -1;

		b1 = 0;
		b2 = 1;

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		typename Problem::FUNCTION f( boost::bind(&ConstrainedSR::obj, this, _1) );
		this->Objectives.push_back(	f );

		typename Problem::FUNCTION c1( boost::bind(&ConstrainedSR::constr1, this, _1) );
		this->InequalityConstraints.push_back(	c1 );

		typename Problem::FUNCTION c2( boost::bind(&ConstrainedSR::constr2, this, _1) );
		this->InequalityConstraints.push_back(	c2 );

		typename Problem::FUNCTION c3( boost::bind(&ConstrainedSR::constr3, this, _1) );
		this->InequalityConstraints.push_back(	c3 );
	}
};

//Problem FON from Pereyra 2009
class FON : public Problem::Interface{
private:
	double obj1(const std::vector< double > &x){
		return 1.0 - exp( -1.0 * ( helpers::square(x[0] - 1.0/sqrt(3.0)) + helpers::square(x[1] - 1.0/sqrt(3.0)) + helpers::square(x[2] - 1.0/sqrt(3.0)) ));
	}
	double obj2(const std::vector< double > &x){
		return 1.0 - exp( -1.0 * ( helpers::square(x[0] + 1.0/sqrt(3.0)) + helpers::square(x[1] + 1.0/sqrt(3.0)) + helpers::square(x[2] + 1.0/sqrt(3.0)) ));
	}
public:
	FON(int DimObj, int DimDesign)	{

		typename Problem::FUNCTION f( boost::bind(&FON::obj1, this, _1) );
		this->Objectives.push_back(	f );
		typename Problem::FUNCTION g( boost::bind(&FON::obj2, this, _1) );
		this->Objectives.push_back(	g );

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( -1.0 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};



//// Toolkit includes. //////////////////////////////////////////////////////

#include "wfgProblems/Toolkit/ExampleProblems.h"
#include "wfgProblems/Toolkit/TransFunctions.h"


//// Used namespaces. ///////////////////////////////////////////////////////

using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;

//WFG2
class WFG5 : public Problem::Interface{
private:
	double obj(const std::vector< double > &x, const int i, const int k, const int M){
		//We're seeing an assertion in the WFGProblems code fail 
		//as the input vectors are "out of bounds" (<0), which is
		//a result of the finite differences scheme.
		
		//return (Problems::WFG2(x, k, M))[i];

		//This should be done with a decorator anyway
		std::vector< double > var_x(x.begin(), x.begin() + this->dimDesign);
		return (Problems::WFG5(var_x, k, M))[i];
	}
public:
	WFG5(int DimObj, int DimDesign)	{
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int k = 1*(DimObj - 1);

		//assert ( k == 1);

		int i;
		for(i=0; i<dimObj; i++){
			typename Problem::FUNCTION f( boost::bind(&WFG5::obj, this, _1, i, k, DimObj) );
			this->Objectives.push_back(	f );
		}
		for(i=0; i<dimDesign; i++){
			//The zero bound is enforced in the WFG code
			//so we perturb it to allow the central difference
			//scheme to work.
			//this->lowerBounds.push_back( 0.0 );
			this->lowerBounds.push_back( 2e-6 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};

//DTLZ2
class DTLZ2 : public Problem::Interface{
private:
	static const double PI = 3.14159;

	double g( const std::vector< double > &x ){
		int i;
		double result = 0.0;
		for(i=this->dimObj-1;i<x.size();i++){
			result += ( x[i] - 0.5 ) * ( x[i] - 0.5 );
		}
		return result;
	}
	double obj(const std::vector< double > &x, const int objNum){
		if( objNum == 0 ) return (1 + this->g( x )) * cos( x[0] * PI/2.0 ) * cos( x[1] * PI/2.0 );
		if( objNum == 1 ) return (1 + this->g( x )) * cos( x[0] * PI/2.0 ) * sin( x[1] * PI/2.0 );
		if( objNum == 2 ) return (1 + this->g( x )) * sin( x[0] * PI/2.0 ) ;
		return 0.0;
	}
    double obj2(const std::vector< double > &x, const int objNum){
		if( objNum == 0 ) return (1 + this->g( x )) * cos( x[0] * PI/2.0 ) * cos( x[1] * PI/2.0 ) * cos( x[2] * PI/2.0 );
		if( objNum == 1 ) return (1 + this->g( x )) * cos( x[0] * PI/2.0 ) * cos( x[1] * PI/2.0 ) * sin( x[2] * PI/2.0 );
		if( objNum == 2 ) return (1 + this->g( x )) * cos( x[0] * PI/2.0 ) * sin( x[1] * PI/2.0 );
        if( objNum == 3 ) return (1 + this->g( x )) * sin( x[0] * PI/2.0 ) ;
		return 0.0;
	}
public:
	DTLZ2(int DimObj, int DimDesign) {
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimObj; i++){
			typename Problem::FUNCTION f( boost::bind(&DTLZ2::obj, this, _1, i) );
            if( dimObj == 4 ) f = boost::bind(&DTLZ2::obj2, this, _1, i) ;
			this->Objectives.push_back(	f );
		}
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( 0.0 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};

//Simulation/Surrogate-based problem
#include "interpolator.hpp"

class Surrogate : public Problem::Interface{
private:
    std::vector< Interpolation::RBF > rbfs;

	double obj(const std::vector< double > &x, const int objNum){
		return this->rbfs[objNum].evaluate( x );
	}
public:
	Surrogate(int DimObj, int DimDesign) {
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

        std::string s[] = {"./data/rms_s.dat", "./data/rms_x.dat", "./data/emit_x.dat"};

		int i;
		for(i=0; i<dimObj; i++){
            //Load the required data
            std::vector< std::vector< double > > data = Interpolation::GetDataFromFile(s[i].c_str()) ;
            //Instatntiate the rbf
            rbfs.push_back( Interpolation::RBF(Interpolation::RBF::RESCALED, data) );
            //bind the function
			typename Problem::FUNCTION f( boost::bind(&Surrogate::obj, this, _1, i) );
            //export to function list
			this->Objectives.push_back(	f );
		}

        //param 1 limits
		//this->lowerBounds.push_back( 0.00026 );
		//this->upperBounds.push_back( 0.00032 );
        this->lowerBounds.push_back( 0.0 );
		this->upperBounds.push_back( 1.0 );
        //param 2 limits
        //this->lowerBounds.push_back( 15.0 );
		//this->upperBounds.push_back( 40.0 );
        this->lowerBounds.push_back( 0.0 );
		this->upperBounds.push_back( 1.0 );
	}
};

class MottaEx1 : public Problem::Interface{
private:
	double obj(const std::vector< double > &x, const int objNum){
		return x[objNum];
	}
    double constr( const std::vector< double > &x, const int constrNum ){
        int other1 = (constrNum + 1)%3;
        int other2 = (constrNum + 2)%3;

        return 1.0/x[other1] + 1.0 / x[other2] - x[constrNum];
    }
public:
	MottaEx1(int DimObj, int DimDesign) {
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimObj; i++){
            //bind the function
			typename Problem::FUNCTION f( boost::bind(&MottaEx1::obj, this, _1, i) );
            //export to function list
			this->Objectives.push_back(	f );
            //bind the function
			typename Problem::FUNCTION g( boost::bind(&MottaEx1::constr, this, _1, i) );
            //export to function list
			this->InequalityConstraints.push_back(	g );
		}

        //param 1 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
        //param 2 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
        //param 3 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
	}
};

class MottaEx2 : public Problem::Interface{
private:
	double obj(const std::vector< double > &x, const int objNum){
		if(objNum==0) return x[0]*x[0]*x[0] + x[1] + 2*x[2];
        if(objNum==1) return x[0] + x[1]*x[1]*x[1] + 2*x[2];
        if(objNum==2) return -x[0]*x[1]*x[2];
	}
    double constr( const std::vector< double > &x, const int constrNum ){
        if(constrNum==0) return x[0]*x[0] + x[1]*x[1] - x[2] - 5;
        if(constrNum==1) return -5*(x[0] + x[1]) + x[2];
    }
public:
	MottaEx2(int DimObj, int DimDesign) {
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimObj; i++){
            //bind the function
			typename Problem::FUNCTION f( boost::bind(&MottaEx2::obj, this, _1, i) );
            //export to function list
			this->Objectives.push_back(	f );
		}
        
		for(i=0; i<2; i++){
            //bind the function
			typename Problem::FUNCTION g( boost::bind(&MottaEx2::constr, this, _1, i) );
            //export to function list
			this->InequalityConstraints.push_back(	g );
		}

        //param 1 limits
        this->lowerBounds.push_back( 0.0 );
		this->upperBounds.push_back( 100.0 );
        //param 2 limits
        this->lowerBounds.push_back( 0.0 );
		this->upperBounds.push_back( 100.0 );
        //param 3 limits
        this->lowerBounds.push_back( 0.0 );
		this->upperBounds.push_back( 100.0 );
	}
};

class MottaEx3 : public Problem::Interface{
private:
    double a[2];
    double b[2];
    double c[2];
    double d[2];
    
    double dist(std::vector<double> x, double* p) {
        return sqrt( (x[0] - p[0])*(x[0] - p[0]) + (x[1] - p[1])*(x[1] - p[1]) );
    }

	double obj(const std::vector< double > &x, const int objNum){
		if(objNum==0) return dist(x, a);
        if(objNum==1) return dist(x, b);
        if(objNum==2) return dist(x, c);
	}
    double constr( const std::vector< double > &x, const int constrNum ){
        return -dist(x, d) + .5;
    }
public:
	MottaEx3(int DimObj, int DimDesign) {
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

        a[0]=0;a[1]=0;
        b[0]=0;b[1]=1;
        c[0]=1;c[1]=0;
        d[0]=1;d[1]=0.4;


		int i;
		for(i=0; i<dimObj; i++){
            //bind the function
			typename Problem::FUNCTION f( boost::bind(&MottaEx3::obj, this, _1, i) );
            //export to function list
			this->Objectives.push_back(	f );
		}
        
		for(i=0; i<1; i++){
            //bind the function
			typename Problem::FUNCTION g( boost::bind(&MottaEx3::constr, this, _1, i) );
            //export to function list
			this->InequalityConstraints.push_back(	g );
		}

        //param 1 limits
        this->lowerBounds.push_back( 0.0 );
		this->upperBounds.push_back( 1.0 );
        //param 2 limits
        this->lowerBounds.push_back( 0 );
		this->upperBounds.push_back( 1.0 );
	}
};

class MottaEx4 : public Problem::Interface{
private:
	double obj(const std::vector< double > &x, const int objNum){
		return x[objNum];
	}
    double constr( const std::vector< double > &x, const int constrNum ){
        int other1 = (constrNum + 1)%4;
        int other2 = (constrNum + 2)%4;
        int other3 = (constrNum + 3)%4;

        return 1.0/x[other1] + 1.0 / x[other2] + 1.0 / x[other3] - x[constrNum];
    }
public:
	MottaEx4(int DimObj, int DimDesign) {
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimObj; i++){
            //bind the function
			typename Problem::FUNCTION f( boost::bind(&MottaEx4::obj, this, _1, i) );
            //export to function list
			this->Objectives.push_back(	f );
            //bind the function
			typename Problem::FUNCTION g( boost::bind(&MottaEx4::constr, this, _1, i) );
            //export to function list
			this->InequalityConstraints.push_back(	g );
		}

        //param 1 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
        //param 2 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
        //param 3 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
        //param 4 limits
        this->lowerBounds.push_back( 0.2 );
		this->upperBounds.push_back( 10.0 );
	}
};

Problem::Interface * Problem::Factory( std::string s, int DimObj, int DimDesign){
	
	Interface * toReturn = NULL;
	
	//Instantiate a simple Basin problem
	if ( s.compare(std::string("BASIN")) == 0 ){
		//This only allows one objective
		assert (DimObj == 1);
		toReturn = new Basin(DimObj, DimDesign);
	}

	//Instantiate a problem to test constraints
	if ( s.compare(std::string("CONST_TEST")) == 0 ){
		//This only allows one objective
		assert (DimObj == 1);
		assert (DimDesign == 2);
		toReturn = new ConstrainedSR(DimObj, DimDesign);
	}

	//Instantiate a FON problem
	if ( s.compare(std::string("FON")) == 0 ){
		assert (DimDesign == 3);
		toReturn = new FON(DimObj, DimDesign);
	}
	
	//Instantiate a WFG5 Problem
	if ( s.compare(std::string("WFG5")) == 0 ){
		//assert (DimDesign % 2 != 0 );
		toReturn = new WFG5(DimObj, DimDesign);
	}

	//Instantiate a DTLZ2 Problem
	if ( s.compare(std::string("DTLZ2")) == 0 ){
		assert (DimObj == 3 || DimObj == 4);
		assert (DimDesign >= 2);
		toReturn = new DTLZ2(DimObj, DimDesign);
	}

    //Instantiate a Simulation Problem
	if ( s.compare(std::string("SURROGATE")) == 0 ){
		assert (DimObj <= 3);
		assert (DimDesign == 2);
		toReturn = new Surrogate(DimObj, DimDesign);
	}
    if ( s.compare(std::string("MOTTA1")) == 0 ){
		assert (DimObj == 3);
		assert (DimDesign == 3);
		toReturn = new MottaEx1(DimObj, DimDesign);
	}
    if ( s.compare(std::string("MOTTA2")) == 0 ){
		assert (DimObj == 3);
		assert (DimDesign == 3);
		toReturn = new MottaEx2(DimObj, DimDesign);
	}
    if ( s.compare(std::string("MOTTA3")) == 0 ){
		assert (DimObj == 3);
		assert (DimDesign == 2);
		toReturn = new MottaEx3(DimObj, DimDesign);
	}
    if ( s.compare(std::string("MOTTA4")) == 0 ){
		assert (DimObj == 4);
		assert (DimDesign == 4);
		toReturn = new MottaEx4(DimObj, DimDesign);
	}
	return toReturn;
}
