#ifndef __SAMPLE_INDIVIDUAL_H__
#define __SAMPLE_INDIVIDUAL_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>


#include "boost/smart_ptr.hpp"

/**
 * \class SampleIndividual
 *  Structure for an individual in the population holding genes
 *  values.
 *
 *  @see Types.h
 */
class SampleIndividual {

public:

    /// representation of genes
    typedef std::vector<double> genes_t;
    /// gene names
    typedef std::vector<std::string> names_t;

    SampleIndividual()
    {}

    SampleIndividual(names_t names)
        : names_m(names)
    {
        genes.resize(names.size(), 0.0);
    }

    /// serialization of structure
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & genes;
        ar & id;
    }

    /// genes of an individual
    genes_t      genes;
    /// id
    unsigned int id = 0;

private:
    /// gene names
    names_t names_m;
};

#endif
