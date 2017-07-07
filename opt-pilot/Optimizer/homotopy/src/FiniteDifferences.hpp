namespace FiniteDifferences {

    //Differencing options
    enum FD_TYPE { CENTRAL, FORWARD };
    //Finite difference parameters
    struct Params_t {
        double step;
        FD_TYPE type;
    };

	template< typename T >
	static double Central(const int dim, std::vector<double> &at, T& func, const double fd_step){
		double tmp = at[dim];

		at[dim] -= fd_step; 
		double left = (func)(at);

		at[dim] = tmp + fd_step; 
		double right = (func)(at);

		at[dim] = tmp;

		return (right - left) / (2*fd_step);
	}
	template< typename T >
	static double Forward(const int dim, std::vector<double> &at, T& func, const double fd_step){
		double tmp = at[dim];

		double center = (func)(at);

		at[dim] += fd_step; 
		double forward = (func)(at);

		at[dim] = tmp;

		return (forward - center) / (fd_step);
	}
	template< typename T >
	static double ForwardOpt(const int dim, std::vector<double> &at, double center, T& func, const double fd_step){
		double tmp = at[dim];

		at[dim] += fd_step; 
		double forward = (func)(at);

		at[dim] = tmp;

		return (forward - center) / (fd_step);
	}
	
    //TODO: Make this a template specialization on par.type
    template<typename Grad, typename T>
	void GradEval(Grad grad, T& func, const std::vector< double > & at, Params_t &par){
		std::vector<double> local_at(at);
        if(par.type == FORWARD) {
            double f_at = (func) (local_at);
            int i;
		    for(i=0; i<at.size(); i++)	grad[i] = ForwardOpt<T>(i, local_at, f_at, func, par.step);
        }
        else{
		    int i;
		    for(i=0; i<at.size(); i++)  grad[i] = Central<T>(i, local_at, func, par.step);   
        }
	}
}