#ifndef __FIXED_VALUE_SIMULATION_H__
#define __FIXED_VALUE_SIMULATION_H__

#include <string>
#include <map>
#include <vector>

#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Simulation/Simulation.h"

/**
 *  \class FixedValueSimulation
 *  \brief Queries a table of values for function evaluations and derivations.
 *
 *  This helper class can be used to test an optimizer algorithm by specifying
 *  a matrix (table) with values for function evaluations and derivation to be
 *  requested by the optimizer.
 */
class FixedValueSimulation : public Simulation {

public:

    /**
     *  Constructor
     *  @param[in] objectives of the optimization problem
     *  @param[in] costraints of the optimization problem
     *  @param[in] params simulation parameters
     *  @param[in] name base name of the simulation
     *  @param[in] comm communicator used for solving the simulation
     *  @param[in] args for accessing the user supplied command line arguments
     */
    FixedValueSimulation(Expressions::Named_t objectives,
                         Expressions::Named_t constraints,
                         Param_t params, std::string name, MPI_Comm comm,
                         CmdArguments_t args);

    /// Destructor.
    virtual ~FixedValueSimulation();

    /**
     *  Dummy function
     */
    void run()
    {}

    /**
     *  Fill reqVars_t with predefined values
     */
    void collectResults();

    /**
     *  Return predefined values
     *  @return reqVars_t containing all requested variables with results
     */
    reqVars_t getResults() { return requested_vars_; }


    void setType(int type)       { type_ = (InfoType_t)type; }
    void setPosition(double pos) { position_ = pos; }

private:
    /// variable dictionary
    std::vector<std::string> dictionary_;

    std::map<std::string, std::string> user_variables_;

    /// by the optimizer requested variables
    reqVars_t requested_vars_;

    NamedExpression_t objectives_;
    NamedExpression_t constraints_;

    /// mapping objective to (position, value) pair
    std::map<std::string, std::map<double, double> > derivative_table_;
    std::map<std::string, std::map<double, double> > evaluation_table_;


    InfoType_t  type_;
    double      position_;

    void readProblem(std::string obj_file, std::string deriv_file);

};

#endif
