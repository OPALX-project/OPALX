#ifndef __SUMERRSQ_H__
#define __SUMERRSQ_H__

#include <map>
#include <string>
#include <fstream>
#include <iterator>

#include "boost/type_traits/remove_cv.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

#include "Util/Types.h"
#include "Util/SDDSReader.h"
#include "Expression/Parser/function.hpp"

class Measurement {
public:
    double spos;
    double measurement;
};

/**
 *  A simple expression computing the sum of all measurement errors (given as
 *  first and third argument) for a variable (second argument) according to
 *
 *  \f[
 *    result = \frac{1}{n} * \sqrt{\sum_{i=0}^n (measurement_i - value_i)^2}
 *  \f]
 *
 */
struct SumErrSq {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {

        std::string measurement_filename = boost::get<std::string>(args[0]);
        var_name_                        = boost::get<std::string>(args[1]);
        stat_filename_                   = boost::get<std::string>(args[2]);

        //FIXME: we could assume measurements don't change
        parseMeasurements(measurement_filename);
        bool is_valid = true;

        boost::scoped_ptr<SDDSReader> sim_stats(new SDDSReader(stat_filename_));
        try {
            sim_stats->parseFile();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }

        double sum = 0;

        foreach(Measurement measurement, measurements_) {
            double sim_value = 0.0;
            try {
                sim_stats->getInterpolatedValue(
                        measurement.spos, var_name_, sim_value);
            } catch(OptPilotException &e) {
                std::cout << "Exception while getting value "
                          << "from SDDS file: " << e.what()
                          << std::endl;
                is_valid = false;
            }
            double val = measurement.measurement - sim_value;
            sum += val * val;
        }

        return boost::make_tuple(sqrt(sum/measurements_.size()), is_valid);
    }

private:

    std::vector<Measurement> measurements_;

    std::string var_name_;
    std::string stat_filename_;

    void parseMeasurements(std::string measurement_filename);

    // define a mapping to arguments in argument vector
    boost::tuple<std::string, std::string, std::string> argument_types;
    enum {
          measurement_filename
        , var_name
        , stat_filename
    } argument_type_id;

};

#endif
