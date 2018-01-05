#ifndef __SDDSVARIABLE_H__
#define __SDDSVARIABLE_H__

#include <string>

#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/SDDSReader.h"
#include "Expression/Parser/function.hpp"


/**
 *  A simple expression to get SDDS (filename = third value) value near a
 *  specific spos (second argument) for a variable (name = first argument).
 */
struct SDDSVariable {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {

        var_name_      = boost::get<std::string>(args[0]);
        spos_          = boost::get<double>(args[1]);
        stat_filename_ = boost::get<std::string>(args[2]);

        bool is_valid = true;

        boost::scoped_ptr<SDDSReader> sim_stats(new SDDSReader(stat_filename_));
        try {
            sim_stats->parseFile();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }

        double sim_value = 0.0;
        try {
            sim_stats->getInterpolatedValue(spos_, var_name_, sim_value);
        } catch(OptPilotException &e) {
            std::cout << "Exception while getting value "
                      << "from SDDS file: " << e.what()
                      << std::endl;
            is_valid = false;
        }

        return boost::make_tuple(sim_value, is_valid);
    }

private:

    std::string var_name_;
    std::string stat_filename_;
    double spos_;

    // define a mapping to arguments in argument vector
    boost::tuple<std::string, double> argument_types;
    enum {
          var_name
        , spos
        , stat_filename
    } argument_type_id;

};

struct sameSDDSVariable {
    sameSDDSVariable(const std::string & base_filename) {
        size_t pos = base_filename.find_last_of("/");
        std::string tmplfile = base_filename;
        if(pos != std::string::npos)
            tmplfile = base_filename.substr(pos+1);
        pos = tmplfile.find(".");
        // std::string simName =
        stat_filename_ = tmplfile.substr(0,pos) + ".stat";
    }

    Expressions::Result_t operator()(client::function::arguments_t args) {
        args.push_back(stat_filename_);
        return var_(args);
    }

private:
    client::function::argument_t stat_filename_;
    SDDSVariable var_;
};

#endif