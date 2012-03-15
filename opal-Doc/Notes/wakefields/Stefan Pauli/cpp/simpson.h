#ifndef SIMPSON_H
#define SIMPSON_H
#include <cassert>

/**
* @brief	Simpson-Integration from the function f from a to b with N steps
*
*
* @param[in] 	f the function to integrate
* @param[in] 	a integrate from a
* @param[in] 	b integrate to b
* @param[in] 	N Number of integration points
* @return 	function value of the integration
*
*/
template<class F> double simpson(F &f, double a, double b, unsigned int N) 
	{
	
	assert(b>a);
	assert(N>0);

	double    result=0;
	double      h=(b-a)/N;

	// boundary values
	result += ( f(a) + 4*f(a+h/2) + f(b) ) / 2.0;

	// values between boundaries
	for ( unsigned int i = 1; i <= N-1; ++i ) {
		result += f(a+i*h) + 2*f(a+(i+0.5)*h);
	}

	result *= h/3.0;

	return result;

	}

#endif //SIMPSON_H
