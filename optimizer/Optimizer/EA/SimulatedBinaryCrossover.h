#include <cmath>

#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"

/**
 *  Deb (1995) Simulated Binary Crossover (SBX)
 *  Respects interval schemata.
 *  Offspring are symmetric around parent solutions.
 */

template <class T> struct SimulatedBinaryCrossover
{
    void crossover(boost::shared_ptr<T> ind1, boost::shared_ptr<T> ind2,
                   CmdArguments_t args) {

        double nu_c = args->getArg<double>("simbin-crossover-nu", 2.0, false);

        for(std::size_t i = 0; i < ind1->genes_m.size(); i++) {

            double ui = (double) rand() / (RAND_MAX + 1.0);
            double beta_qi = 0.0;
            if(ui <= 0.5) {
                beta_qi = pow(2 * ui, 1.0/(nu_c + 1.0));
            } else {
                beta_qi = pow(1.0/(2 * (1.0 - ui)), 1.0/(nu_c + 1.0));
            }

            double ming = std::min(ind1->genes_m[i], ind2->genes_m[i]);
            double maxg = std::max(ind1->genes_m[i], ind2->genes_m[i]);

            ind1->genes_m[i] = 0.5 * ((1 + beta_qi) * ming + (1 - beta_qi) * maxg);
            ind2->genes_m[i] = 0.5 * ((1 - beta_qi) * ming + (1 + beta_qi) * maxg);
        }
    }
};
