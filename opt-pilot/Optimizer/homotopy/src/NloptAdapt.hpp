#ifndef NloptAdapt_hpp
#define NloptAdapt_hpp

/*
Adapter wrapping a function object in the Nlopt library interface
*/

#include "FiniteDifferences.hpp"

template<typename T>
class NloptAdapt{
private:
	//Templated function object:  could be boost::function, lambda, functor, etc.
	T & obj;
	std::vector< T > * EqConstr;
	std::vector< T > * InEqConstr;
    
    bool &Validity;

    FiniteDifferences::Params_t fd_params;

	static double eval( T& func, const std::vector<double> &at){	return (func)(at);	}

	static void ConstrIface(unsigned m, double* result, unsigned n, const double* x, double* grad, std::vector< T > &Constr, NloptAdapt *nl){
		//Copy array data to vector
		std::vector< double > at(n);
		int i;
		for(i=0; i<n; i++){ at[i] = x[i]; }

        //Could abstract this into a universal OptController interface
        nl->Validity = true;
		
        //Evaluate constraints and gradients
		for(i=0; i<m; i++){
			result[i] = eval(Constr[i], at);
            if(grad) FiniteDifferences::GradEval<double *, T>(grad+i*n, Constr[i], at, nl->fd_params);
		}

        //If all the evaluations weren't valid
        //force stop the evaluator
        if(! nl->Validity) throw nlopt::forced_stop();
	}

public:
    NloptAdapt(T &Obj, std::vector<T> *EqConstr, std::vector<T> *InEqConstr, bool& validity, FiniteDifferences::Params_t fd_params) : 
                obj(Obj), EqConstr(EqConstr), InEqConstr(InEqConstr), Validity( validity ), fd_params( fd_params ) {}

	//C-style callback interfaces expected by Nlopt
	static double ObjIface(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
		assert ( my_func_data != NULL );
		NloptAdapt *nl = reinterpret_cast<NloptAdapt *>(my_func_data);

        //Could abstract this into a universal OptController interface
        nl->Validity = true;

        //Evaluate objectives and gradient
		if(!grad.empty()){
            FiniteDifferences::GradEval<std::vector<double> &, T>(grad, nl->obj,  x,  nl->fd_params );
		}
		double result = eval(nl->obj, x);

        //If all the evaluations weren't valid
        //force stop the evaluator
        if(! nl->Validity) throw nlopt::forced_stop();
        else return result;
	}

	//C-style interface for multiple constriaints
	static void EqConstrIface(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data){
		assert ( data != NULL );
		NloptAdapt *nl = reinterpret_cast<NloptAdapt *>(data);
		ConstrIface( m, result, n, x, grad, *nl->EqConstr, nl);
	}
	static void InEqConstrIface(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data){
		assert ( data != NULL );
		NloptAdapt *nl = reinterpret_cast<NloptAdapt *>(data);
		ConstrIface( m, result, n, x, grad, *nl->InEqConstr, nl);
	}

    //C-style interface for single constraints
	/*static double ConstrIface(const std::vector<double> &x, std::vector<double> &grad, void *data){	}*/
};

#endif