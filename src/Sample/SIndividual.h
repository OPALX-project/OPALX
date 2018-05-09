#ifndef __SAMPLE_INDIVIDUAL_H__
#define __SAMPLE_INDIVIDUAL_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>


#include "boost/smart_ptr.hpp"

/**
 * \class Individual
 *  Structure for an individual in the population holding genes
 *  values.
 *
 *  @see Types.h
 */
class SIndividual {

public:

    /// representation of genes
    typedef std::vector<double> genes_t;
    /// gene names
    typedef std::vector<std::string> names_t;

    SIndividual()
    {}

    /// create a new SIndividual
    SIndividual(names_t names)
        : names_m(names)
    {
        genes.resize(names.size(), 0.0);
    }

    /// copy another SIndividual
    SIndividual(boost::shared_ptr<SIndividual> individual) {
        genes         =       genes_t(individual->genes);
        names_m       =       names_t(individual->names_m);
        id            = individual->id;
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