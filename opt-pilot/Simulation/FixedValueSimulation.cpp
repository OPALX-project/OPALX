#include <iostream>
#include <sstream>
#include <string.h>

#include "Simulation/FixedValueSimulation.h"

#include "Util/OptPilotException.h"

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

FixedValueSimulation::FixedValueSimulation(Expressions::Named_t objectives,
                 Expressions::Named_t constraints, Param_t params,
                 std::string name, MPI_Comm comm, CmdArguments_t args)
               : Simulation(args)
               , objectives_(objectives)
               , constraints_(constraints)
               , params_(params)
{
    type_ = NONE;
    position_ = 0.0;

    std::string objective_file  = "";
    std::string derivation_file = "";

    objective_file = args->getArg<double>("obj-file", true);
    derivative_file = args->getArg<double>("deriv-file", true);

    readProblem(objective_file, derivative_file);
}


FixedValueSimulation::~FixedValueSimulation() {

    derivation_table_.clear();
    evaluation_table_.clear();

    dictionary_.clear();

    requested_vars_.clear();
    user_variables_.clear();
}


void FixedValueSimulation::readProblem(std::string obj_file,
                                       std::string deriv_file) {

    //FIXME: read file into tables

}


void FixedValueSimulation::collectResults() {

    requested_vars_.clear();

    foreach(NamedExpression_t &objective, objectives_) {

        double value = 0.0;
        std::string name = objective->first;

        switch(type_) {
            case EVALUATE:
                value = (evaluation_table_[name])[position_];
                break;

            case DERIVATE:
                if (derivation_table_.size() > 0)
                    value = (derivation_table_[name])[position_];
                else {
                    //FIXME: stencil for derivation
                    value = 0.0;
                }
                break;

            default:
                throw OptPilotException("unknown evaluation type");
                break;
        }

        reqVarInfo_t tmps = {type_, value, true};
        requested_vars_.insert(
                std::pair<std::string, reqVarInfo_t>(it->first, tmps));
    }
}

