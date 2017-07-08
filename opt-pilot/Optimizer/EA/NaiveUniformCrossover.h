#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"

/// decide for each gene if swapped with other gene
template <class T> struct NaiveUniformCrossover
{
    void crossover(boost::shared_ptr<T> ind1, boost::shared_ptr<T> ind2,
                   CmdArguments_t args) {

        genes_t genes_ind2 = ind2->genes;

        for(int i = 0; i < ind1->genes.size(); i++) {
            int choose = (int) (2.0 * (double) rand() / (RAND_MAX + 1.0));
            if(choose == 1) {
                ind2->genes[i] = ind1->genes[i];
                ind1->genes[i] = genes_ind2[i];
            }
        }
    }
};
