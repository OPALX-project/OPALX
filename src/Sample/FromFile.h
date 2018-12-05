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
        , filename_m(filename)
        , dvarName_m(dvarName)
    {
        // we need to count the number of lines
        std::ifstream in(filename_m);

        if ( !in.is_open() ) {
            throw OpalException("FromFile()",
                                "Couldn't open file \"" + filename_m + "\".");
        }

        int nLines = std::count(std::istreambuf_iterator<char>(in),
                                std::istreambuf_iterator<char>(), '\n');

        // make sure we do not count empty lines at end
        in.seekg(-1, std::ios_base::end);
        std::size_t pos =  in.tellg();

        std::string line;
        std::getline(in, line);

        while ( line.empty() ) {
            --nLines;
            --pos;
            in.seekg(pos, std::ios_base::beg);
            std::getline(in, line);
        }

        if ( nLines < 0 )
            throw OpalException("FromFile()", "Empty file \"" + filename_m + "\".");

        globalSize_m = nLines;

        in.close();
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = getNext();
    }

    void allocate(const CmdArguments_t& args, const Comm::Bundle_t& comm) {
        std::ifstream in(filename_m);

        if ( !in.is_open() ) {
            throw OpalException("FromFile()",
                                "Couldn't open file \"" + filename_m + "\".");
        }

        std::string header;
        std::getline(in, header);
        std::istringstream iss(header);
        std::vector<std::string> dvars({std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}});
        size_t j = 0;
        for (const std::string str: dvars) {
            if (str == dvarName_m) break;
            ++ j;
        }

        if (j == dvars.size()) {
            throw OpalException("FromFile()",
                                "Couldn't find the dvar '" + dvarName_m + "' in the file '" + filename_m + "'");
        }

        int nSamples = args->getArg<int>("nsamples", true);
        int nMasters = args->getArg<int>("num-masters", true);

        int nLocSamples = nSamples / nMasters;
        int rest = nSamples - nMasters * nLocSamples;

        int id = comm.island_id;
        if ( id < rest )
            nLocSamples++;

        int skip = 0;

        if ( rest == 0 )
            skip = nLocSamples * id;
        else {
            if ( id < rest ) {
                skip = nLocSamples * id;
            } else {
                skip = (nLocSamples + 1) * rest + (id - rest) * nLocSamples;
            }
        }

        while ( skip > 0 ) {
            skip--;
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        std::string line;
        std::getline(in, line);
        while (nLocSamples-- > 0) {
            std::istringstream iss(line);
            std::vector<std::string> numbers({std::istream_iterator<std::string>{iss},
                                              std::istream_iterator<std::string>{}});

            chain_m.push_back(std::stod(numbers[j]));

            std::getline(in, line);
        }
        in.close();
    }

    double getNext() {
        double sample = chain_m[n_m];
        incrementCounter();

        return sample;
    }

    unsigned int getSize() const {
        return globalSize_m;
    }

private:
    std::vector<double> chain_m;
    unsigned int n_m;
    size_t counter_m;
    size_t mod_m;
    std::string filename_m;
    std::string dvarName_m;

    unsigned int globalSize_m;

    void incrementCounter() {
        ++ counter_m;
        if (counter_m % mod_m == 0)
            ++ n_m;
        if (n_m >= chain_m.size())
            n_m = 0;
    }
};

#endif