#include "Simulation/FunctionEvaluator.h"

#include "boost/algorithm/string.hpp"
#include "boost/variant.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

FunctionEvaluator::FunctionEvaluator(Expressions::Named_t objectives,
                 Expressions::Named_t constraints, Param_t params,
                 std::string name, MPI_Comm comm, CmdArguments_t args)
               : Simulation(args)
               , objectives_(objectives)
               , constraints_(constraints)
               , params_(params)
               , comm_(comm)
{}


FunctionEvaluator::~FunctionEvaluator() {
    requestedVars_.clear();
}


void FunctionEvaluator::collectResults() {

    requestedVars_.clear();
    Expressions::Named_t::iterator it;

    // evaluate all objectives
    for(it = objectives_.begin(); it != objectives_.end(); it++) {

        Expressions::Expr_t *objective = it->second;
        Expressions::Result_t result = objective->evaluate(params_);

        std::vector<double> values;
        values.push_back(boost::get<0>(result));
        bool is_valid = boost::get<1>(result);

        reqVarInfo_t tmps = {EVALUATE, values, is_valid};
        requestedVars_.insert(
                std::pair<std::string, reqVarInfo_t>(it->first, tmps));

    }

    // .. and constraints
    for(it = constraints_.begin(); it != constraints_.end(); it++) {

        Expressions::Expr_t *constraint = it->second;
        Expressions::Result_t result = constraint->evaluate(params_);

        std::vector<double> values;
        values.push_back(boost::get<0>(result));
        bool is_valid = boost::get<1>(result);

        //FIXME: hack to evaluate LHS and RHS
        std::string constr_str = constraint->toString();
        std::vector<std::string> split;
        boost::split(split, constr_str, boost::is_any_of("<>!="),
                     boost::token_compress_on);
        std::string lhs_constr_str = split[0];
        std::string rhs_constr_str = split[1];
        boost::trim_left_if(rhs_constr_str, boost::is_any_of("="));

        functionDictionary_t funcs = constraint->getRegFuncs();
        boost::scoped_ptr<Expressions::Expr_t> lhs(
            new Expressions::Expr_t(lhs_constr_str, funcs));
        boost::scoped_ptr<Expressions::Expr_t> rhs(
            new Expressions::Expr_t(rhs_constr_str, funcs));

        Expressions::Result_t lhs_res = lhs->evaluate(params_);
        Expressions::Result_t rhs_res = rhs->evaluate(params_);

        values.push_back(boost::get<0>(lhs_res));
        values.push_back(boost::get<0>(rhs_res));

        reqVarInfo_t tmps = {EVALUATE, values, is_valid};
        requestedVars_.insert(
                std::pair<std::string, reqVarInfo_t>(it->first, tmps));

    }
}

