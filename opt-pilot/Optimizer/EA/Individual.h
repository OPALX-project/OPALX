#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <vector>

#include "boost/smart_ptr.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

/**
 * \class Individual
 *  Structure for an individual in the population holding genes and objective
 *  values.
 *
 *  @see Types.h
 */
class Individual {

public:

    /// representation of genes
    typedef std::vector<double> genes_t;
    /// objectives array
    typedef std::vector<double> objectives_t;
    /// bounds on design variables
    typedef std::vector< std::pair<double, double> > bounds_t;

    Individual()
    {}

    /// create a new individual and initialize with random genes
    Individual(bounds_t gene_bounds) {
        bounds = gene_bounds;
        genes.resize(bounds.size(), 0.0);
        objectives.resize(bounds.size(), 0.0);

        for(size_t i=0; i < bounds.size(); i++) {
            new_gene(i);
        }
    }

    /// copy another individual
    Individual(boost::shared_ptr<Individual> individual) {
        bounds     = bounds_t(individual->bounds);
        genes      = genes_t(individual->genes);
        objectives = objectives_t(individual->objectives);
    }

    /// serialization of structure
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & genes;
        ar & objectives;
        ar & id;
    }

    /// initialize the gene with index gene_idx with a new random value
    /// contained in the specified gene boundaries
    double new_gene(size_t gene_idx) {
        double max = std::max(bounds[gene_idx].first, bounds[gene_idx].second);
        double min = std::min(bounds[gene_idx].first, bounds[gene_idx].second);
        double delta = fabs(max - min);
        genes[gene_idx] = rand() / (RAND_MAX + 1.0) * delta + min;
        return genes[gene_idx];
    }


    /// genes of an individual
    genes_t      genes;
    /// values of objectives of an individual
    objectives_t objectives;
    /// id
    unsigned int id;

    /// bounds on each gene
    bounds_t bounds;
};

#endif