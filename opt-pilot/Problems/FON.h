#ifndef __FON_H__
#define __FON_H__

#include <string>
#include <map>
#include <set>
#include <iterator>

#include <math.h>

#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"

#include "Util/Types.h"
#include "Util/OptPilotException.h"
#include "Expression/Parser/function.hpp"


struct FON {

    inline static double square( double x )	{ return (x) * (x); }

    static double obj1( double x1, double x2, double x3 ) {
        return 1.0 - exp( -1.0 * ( square(x1 - 1.0/sqrt(3.0)) +
            square(x2 - 1.0/sqrt(3.0)) + square(x3 - 1.0/sqrt(3.0)) ));
    }
    static double obj2( double x1, double x2, double x3 ) {
        return 1.0 - exp( -1.0 * ( square(x1 + 1.0/sqrt(3.0)) +
            square(x2 + 1.0/sqrt(3.0)) + square(x3 + 1.0/sqrt(3.0)) ));
    }

    Expressions::Result_t operator()(client::function::arguments_t args) {

        // Parse the args
        double x1  = boost::get<double>(args[0]);
        double x2  = boost::get<double>(args[1]);
        double x3  = boost::get<double>(args[2]);
        double obj = boost::get<double>(args[3]);

        // evaluate result depending on obj argument
        double result = 0.0;
        if( obj == 1.0 ) result = obj1( x1, x2, x3 );
        else result = obj2( x1, x2, x3 );

        return boost::make_tuple(result, true);
    }
};

#endif
