#ifndef __VARIATOR_H__
#define __VARIATOR_H__

#include <string>
#include <vector>
#include <map>
#include <utility>

#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Optimizer/EA/Population.h"
#include "Optimizer/Optimizer.h"

template<
      class ind_t
    , template <class> class CrossoverOperator
    , template <class> class MutationOperator
>
class Variator : public CrossoverOperator<ind_t>,
                 public MutationOperator<ind_t>
{

public:

    Variator(size_t sizeInitial, size_t nrSelectedParents,
             size_t nrChilderToProduce, size_t dim,
             Optimizer::bounds_t dVarBounds, CmdArguments_t args) {

        //FIXME: pass population as arg to variator
        //boost::shared_ptr< Population<ind_t> >
        population_m.reset(new Population<ind_t>());
        dVarBounds_m = dVarBounds;

        sizeInitial_m = sizeInitial;
        nrSelectedParents_m = nrSelectedParents;
        nrChilderToProduce_m = nrChilderToProduce;

        mutationProbability_m =
            args->getArg<double>("mutation-probability", 0.5);

        recombinationProbability_m =
            args->getArg<double>("recombination-probability", 0.5);

        args_ = args;
    }

    ~Variator() {
    }

    //FIXME access population from outside
    boost::shared_ptr< Population<ind_t> > population() {
        return population_m;
    }

    /// create an initial population
    void initial_population() {
        for(size_t i = 0; i < sizeInitial_m; i++)
            new_individual();

    }

    /// set an individual as individual: replace with a new individual
    void infeasible(boost::shared_ptr<ind_t> ind) {
        population_m->remove_individual(ind);
        new_individual();
    }

    /// returns false if all individuals have been evaluated
    bool hasMoreIndividualsToEvaluate() {
        return !individualsToEvaluate_m.empty();
    }

    /// return next individual to evaluate
    boost::shared_ptr<ind_t> popIndividualToEvaluate() {
        unsigned int ind = individualsToEvaluate_m.front();
        individualsToEvaluate_m.pop();
        return population_m->get_staging(ind);
    }

    /** Performs variation (recombination and mutation) on a set of parent
     *  individuals.
     *
     *  @param[in] parents
     */
    void variate(std::vector<unsigned int> parents) {

        // copying all individuals from parents
        std::vector<unsigned int>::iterator itr;
        for(itr = parents.begin(); itr != parents.end(); itr++) {
            new_individual( population_m->get_individual(*itr) );
        }

        // only variate new offspring, individuals in staging area have been
        // variated already
        std::queue<unsigned int> tmp(individualsToEvaluate_m);
        while(!tmp.empty()) {

            // pop first individual
            unsigned int idx = tmp.front(); tmp.pop();
            boost::shared_ptr<ind_t> a = population_m->get_staging(idx);

            // handle special case where we have an uneven number of offspring
            if(tmp.empty()) {
                if (drand(1) <= mutationProbability_m)
                    this->mutate(a, args_);
                break;
            }

            // and second if any
            idx = tmp.front(); tmp.pop();
            boost::shared_ptr<ind_t> b = population_m->get_staging(idx);

            // do recombination
            if(drand(1) <= recombinationProbability_m) {
                this->crossover(a, b, args_);
            }

            // do mutation
            if (drand(1) <= mutationProbability_m) {
                this->mutate(a, args_);
                this->mutate(b, args_);
            }
        }
    }


protected:

    /// create a new individual
    void new_individual() {
        boost::shared_ptr<ind_t> ind(new ind_t(dVarBounds_m));
        individualsToEvaluate_m.push( population_m->add_individual(ind) );
    }

    /// copy an individual
    void new_individual(boost::shared_ptr<ind_t> ind) {
        boost::shared_ptr<ind_t> return_ind(new ind_t(ind));
        individualsToEvaluate_m.push(
            population_m->add_individual(return_ind) ) ;
    }


private:

    /// population of individuals
    boost::shared_ptr< Population<ind_t> > population_m;
    /// number of individuals in initial population
    size_t sizeInitial_m;
    ///
    size_t nrSelectedParents_m;
    ///
    size_t nrChilderToProduce_m;
    /// number of genes
    size_t dim_m;

    /// user specified command line arguments
    CmdArguments_t args_;

    /// keep a queue of individuals that have to be evaluated
    std::queue<unsigned int> individualsToEvaluate_m;

    /// bounds on design variables
    Optimizer::bounds_t dVarBounds_m;

    /// probability of applying the mutation operator
    double mutationProbability_m;
    /// probability of applying the recombination operator
    double recombinationProbability_m;

    /**
     *  Get a random double between [0, range]
     *  @param[in] range of random number
     *  @return random double value between [0, range]
     */
    double drand(double range) {
        return (range * (double) rand() / (RAND_MAX + 1.0));
    }

    /**
     *  Get a random integer between [0, range]
     *  @param[in] range of random number
     *  @return random integer value between [0, range]
     */
    int irand(int range) {
        return (int) ((double) range * (double) rand() / (RAND_MAX + 1.0));
    }

};

#endif
