#ifndef __PROBEVARIABLE_H__
#define __PROBEVARIABLE_H__

#include <string>

#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/ProbeReader.h"
#include "Expression/Parser/function.hpp"

struct ProbeVariable {
    
    static const std::string name;
    
    Expressions::Result_t operator()(client::function::arguments_t args) {
        
        var_name_       = boost::get<std::string>(args[0]);
        id_             = boost::get<double>(args[1]); //FIXME Can't we use integer?
        probe_filename_ = boost::get<std::string>(args[2]);
        
        bool is_valid = true;
        
        boost::scoped_ptr<ProbeReader> sim_probe(new ProbeReader(probe_filename_));
        
        try {
            sim_probe->parseFile();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }
        
        double sim_value = 0.0;
        try {
            sim_probe->getVariableValue(id_, var_name_, sim_value);
        } catch(OptPilotException &e) {
            std::cout << "Exception while getting value "
                      << "from Probe file: " << e.what()
                      << std::endl;
            is_valid = false;
        }
        
        return boost::make_tuple(sim_value, is_valid);
    }
    
private:
    std::string var_name_;
    int id_;
    std::string probe_filename_;

    // define a mapping to arguments in argument vector
    boost::tuple<std::string, int, std::string> argument_types;
    enum {
          var_name
        , id
        , probe_filename
    } argument_type_id;
};

#endif
