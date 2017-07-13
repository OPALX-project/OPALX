#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"

/// Mutate each gene with probability p
template <class T> struct IndependentBitMutation
{
    void mutate(boost::shared_ptr<T> ind, CmdArguments_t args) {

        const double probability =
            args->getArg<double>("gene-mutation-probability", 0.5);

        for(size_t i = 0; i < ind->genes.size(); i++) {
            double rval = static_cast<double>(rand() / (RAND_MAX + 1.0));
            if(rval < probability) {
                ind->new_gene(i);
            }
        }
    }
};