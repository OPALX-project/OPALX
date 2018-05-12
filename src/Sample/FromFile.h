#ifndef OPAL_SEQUENCE_H
#define OPAL_SEQUENCE_H

#include "Sample/SamplingMethod.h"
#include "Utilities/OpalException.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iterator>

#include <vector>

/**
 * Parse file that contains design variable values.
 * Each column belongs to a design variable.
 * The first line is considered as header and consists of the
 * design variable name. The name has to agree with the string
 * in the input file.
 */
class FromFile : public SamplingMethod
{

public:

    FromFile(const std::string &filename, const std::string &dvarName, size_t modulo)
        : n_m(0)
        , counter_m(0)
        , mod_m(modulo)
    {
        std::ifstream in(filename);

        if ( !in.is_open() ) {
            throw OpalException("FromFile()",
                                "Couldn't open file \"" + filename + "\".");
        }

        std::string header;
        std::getline(in, header);
        std::istringstream iss(header);
        std::vector<std::string> dvars({std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}});
        size_t j = 0;
        for (const std::string str: dvars) {
            if (str == dvarName) break;
            ++ j;
        }

        if (j == dvars.size()) {
            throw OpalException("FromFile()",
                                "Couldn't find the dvar '" + dvarName + "' in the file '" + filename + "'");
        }

        std::string line;
        std::getline(in, line);
        while (in.good() && !in.eof()) {
            std::istringstream iss(line);
            std::vector<std::string> numbers({std::istream_iterator<std::string>{iss},
                                              std::istream_iterator<std::string>{}});

            chain_m.push_back(std::stod(numbers[j]));

            std::getline(in, line);
        }
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = getNext();
    }

    double getNext() {
        double sample = chain_m[n_m];
        incrementCounter();

        return sample;
    }

private:
    std::vector<double> chain_m;
    unsigned int n_m;
    size_t counter_m;
    size_t mod_m;

    void incrementCounter() {
        ++ counter_m;
        if (counter_m % mod_m == 0)
            ++ n_m;
        if (n_m >= chain_m.size())
            n_m = 0;
    }
};

#endif