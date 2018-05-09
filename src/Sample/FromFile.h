#ifndef OPAL_SEQUENCE_H
#define OPAL_SEQUENCE_H

#include "Sample/SamplingMethod.h"
#include "Utilities/OpalException.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iterator>

#include <vector>

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
    }
    
    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        
        std::string line;
        std::getline(in_m, line);
        
        // 9. Mai 2018
        // https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string
        std::istringstream iss(line);
        std::vector<std::string> numbers{std::istream_iterator<std::string>{iss},
                                         std::istream_iterator<std::string>{}};
        
        if ( i > numbers.size() ) {
            throw OpalException("FromFile()",
                                "Number of DVARS is greater than file columns.");
        }
        
        ind->genes[i] = std::stod(numbers[i]);
        
        ++line_m;
    }
    
    ~FromFile() {
        in_m.close();
    }
    
private:
    std::ifstream in_m;
    std::string filename_m;
    int line_m;
};

#endif
