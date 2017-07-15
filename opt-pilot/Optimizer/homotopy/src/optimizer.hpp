#ifndef optimizer_hpp
#define optimizer_hpp

#include "HomotopyTypes.h"
namespace Homotopy {

  class optimizer {
  public:
    //virtual double RunFrom(std::vector< double > &) = 0;
    virtual int RunFrom(designVars_t &) = 0;
    virtual void RefreshConstraints() = 0;
    enum EXIT_COND { SUCCESS, RERUN, ERROR };
  };

#include "boost/bind.hpp"
#include "boost/function.hpp"

  //Maintains validity state for groups of evaluations
  //Can be extended to interface with Comms?
  class EvaluationController {
    typedef boost::function<double (const designVars_t &, bool &)> functionWithStatus_t;
    functionWithStatus_t f;
    double evalAdapter( const designVars_t &x ) {
      bool tmp = false;
      double result( (this->f)( x, tmp ) );
      this->valid &= tmp;
      return result;
    }
  public:
    bool valid;
    function_t objFunc;
    EvaluationController ( functionWithStatus_t F ) :  f(F), valid(false) {
      this->objFunc = ( boost::bind( &EvaluationController::evalAdapter, this, _1 ) );
    }
    void Reset() { this->valid = true; }
    bool isValid() { return this->valid; }
  };

#include "Scalarization.hpp"
#include "nlopt.hpp"
#include "NloptAdapt.hpp"

#include <string>
#include <stdio.h>

  //NloptBased optimizer
  class OptNlopt : public optimizer {
  private:    
    double tolerance;
    ScalarizationInterface *S;
    EvaluationController E;
    NloptAdapt< function_t > NA;
    nlopt::opt opt;
    std::vector< double > EqTolerances;
    std::vector< double > InEqTolerances;

    //NOTE: In addition to being easy to change optimizers,
    //we can try more traditional methods using the augmented
    //lagrangian approach

    nlopt::opt local;

    void initialize(){
      //Set the state of the optimizer:
	
      //Boundary values
      if (!S->lowerBounds.empty()) opt.set_lower_bounds(S->lowerBounds);
      if (!S->upperBounds.empty()) opt.set_upper_bounds(S->upperBounds);
        

      //Set the stop conditions
      //this requires some attention

      //Local optimizer settings
      //local.set_xtol_abs(1e-3);
      //local.set_ftol_abs(1e-3);
      //opt.set_local_optimizer( local );

      opt.set_xtol_abs(tolerance);
      opt.set_ftol_abs(tolerance);
      //opt.set_xtol_rel(tolerance);
      //opt.set_ftol_rel(tolerance);

      //Pass a scalarized function through the 
      //Nlopt Adapter to the Nlopt object
      opt.set_min_objective(&NloptAdapt< typename Problem::FUNCTION >::ObjIface, (void*)(&(this->NA)));

      //Likewise pass the constraints to the optimizer
      if ( !S->EqualityConstraints.empty() ) {
	opt.add_equality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::EqConstrIface, (void*)(&(this->NA)), this->EqTolerances);
      }
      if ( !S->InequalityConstraints.empty() ) {
	opt.add_inequality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::InEqConstrIface, (void*)(&(this->NA)), this->InEqTolerances);
      }
    }

  public:
    //Using SLSQP as the default optimization algorithm
    OptNlopt(ScalarizationInterface *s, double tolerance, FiniteDifferences::Params_t fd_par) : 
      optimizer(), S(s), 
      E( boost::bind( &ScalarizationInterface::operator(), S, _1, _2 ) ), 
      NA(E.objFunc, &S->EqualityConstraints, &S->InequalityConstraints, E.valid, fd_par ), tolerance(tolerance),
      opt(nlopt::LD_SLSQP, S->dimDesign), 
      EqTolerances(S->EqualityConstraints.size(), tolerance),
      InEqTolerances(S->InequalityConstraints.size(), tolerance){

      this->initialize();
    }

    OptNlopt(ScalarizationInterface *s, double tolerance, FiniteDifferences::Params_t fd_par, nlopt::algorithm alg) : 
      optimizer(), S(s), 
      E( boost::bind( &ScalarizationInterface::operator(), S, _1, _2 ) ), 
      NA(E.objFunc, &S->EqualityConstraints, &S->InequalityConstraints, E.valid, fd_par ), tolerance(tolerance),
      opt(alg, S->dimDesign), 
      //local(nlopt::LN_COBYLA, S->dimDesign), opt(nlopt::LD_AUGLAG, S->dimDesign),
      EqTolerances(S->EqualityConstraints.size(), tolerance),
      InEqTolerances(S->InequalityConstraints.size(), tolerance){

      this->initialize();
    }

    ~OptNlopt() { }

    int RunFrom(std::vector< double > &x) {
      int result = SUCCESS;
      double minf;
      try {
	opt.optimize(x, minf);
      }
      //catch (const std::exception& ex) {
      catch (const nlopt::forced_stop& ex) {
	//Needs to run again
	result = RERUN;
      }
      catch (const nlopt::roundoff_limited& ex) {
	printf("WARNING: threw nlopt::roundoff_limited\n");
	result = ERROR;
      }
      catch (const std::invalid_argument& ex) {
	printf("WARNING: threw std::invalid_argument with %s\n input: ", ex.what());
	for(unsigned i=0;i<x.size();i++) printf("%lf ", x[i]);
	printf("\n");
	result = ERROR;
      }
      catch (const std::runtime_error& ex) {
	printf("WARNING: threw std::runtime with %s\n", ex.what());
	result = ERROR;
      }

      return result;
    }

    void RefreshConstraints(){
      //printf("Called right func\n");
      opt.remove_equality_constraints();
      opt.remove_inequality_constraints();
      //printf("Called right func\n");
      std::vector< double > E( this->S->EqualityConstraints.size(), this->tolerance) ;
      this->EqTolerances.assign(E.begin(), E.end());
      std::vector< double > I( this->S->InequalityConstraints.size(), this->tolerance);
      this->InEqTolerances.assign(I.begin(), I.end() );
      //printf("Called right func\n");
      if ( !S->EqualityConstraints.empty() ) {
	this->opt.add_equality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::EqConstrIface, (void*)(&(this->NA)), this->EqTolerances);
      }
      if ( !S->InequalityConstraints.empty() ) {
	this->opt.add_inequality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::InEqConstrIface, (void*)(&(this->NA)), this->InEqTolerances);
      }
      //printf("Called right func\n");
    }
  };

}

#endif
