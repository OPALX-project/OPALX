#ifndef __PROBLEMINTERFACE_H__
#define __PROBLEMINTERFACE_H__

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

#include "Problems.h"
//#include "homotopy.hpp"

class ProblemInterface {
private:
    typename Problem::Interface * P;

public:
    ProblemInterface( std::string name, int Obj, int Design ) {
        this->P = Problem::Factory( name, Obj, Design );
    }

    Expressions::Result_t operator()(client::function::arguments_t args) {

        //First arg is the objective number
        int obj = boost::get<double>(args[0]);

        // Parse the args into the design vars
        std::vector< double > x( P->dimDesign );
        int i;
        for(i=0;i<P->dimDesign;i++){
            x[i]  = boost::get<double>(args[i + 1]);
        }
        
        // evaluate result depending on obj argument
        double result = (this->P->Objectives[obj])( x );

        return boost::make_tuple(result, true);
    }
};

#endif
