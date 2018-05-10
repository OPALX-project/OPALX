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
    
    FromFile(std::string filename)
        : filename_m(filename)
        , line_m(0)
        
    {
        in_m.open(filename);
        
        if ( !in_m.is_open() ) {
            throw OpalException("FromFile()",
                                "Couldn't open file \"" + filename + "\".");
        }
        
        std::string header;
        std::getline(in_m, header);
        std::istringstream iss(header);
        dvars_m = {std::istream_iterator<std::string>{iss},
                   std::istream_iterator<std::string>{}};
        line_m++;
    }
    
    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        
        std::string line;
        std::getline(in_m, line);
        
        // 9. Mai 2018
        // https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string
        std::istringstream iss(line);
        std::vector<std::string> numbers{std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>{}};
        
        // the file column does not need to agree with gene index i
        size_t j = getColumn_m(ind->getName(i));
        
        ind->genes[i] = std::stod(numbers[j]);
        
        ++line_m;
    }
    
    ~FromFile() {
        in_m.close();
    }
    
private:
    
    size_t getColumn_m(std::string name) {
        auto res = std::find(std::begin(dvars_m), std::end(dvars_m), name);
        
        if (res == std::end(dvars_m)) {
            throw OpalException("FromFile::getColumn()",
                                "Variable '" + name + "' not contained.");
        }
        return std::distance(std::begin(dvars_m), res);
    }
    
private:
    std::ifstream in_m;
    std::string filename_m;
    int line_m;
    
    std::vector<std::string> dvars_m;
};

#endif
