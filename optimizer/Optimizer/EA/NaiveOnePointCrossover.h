#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"

template <class T> struct NaiveOnePointCrossover
{
    void crossover(boost::shared_ptr<T> ind1, boost::shared_ptr<T> ind2,
                   CmdArguments_t args) {

        typedef typename T::genes_t genes_t;
        genes_t genes_ind2;
        genes_ind2 = ind2->genes;

        // determine crossover position u.a.r.
        size_t position = static_cast<size_t>(
            ((double) ind1->genes.size() * (double) rand() / (RAND_MAX + 1.0))
        );

        for(size_t i = position; i < ind1->genes.size(); i++) {
            ind2->genes[i] = ind1->genes[i];
            ind1->genes[i] = genes_ind2[i];
        }
    }
};
