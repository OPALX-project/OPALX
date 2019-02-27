#ifndef OPAL_SAMPLE_WEIGHTED_SEQUENCE_H
#define OPAL_SAMPLE_WEIGHTED_SEQUENCE_H

#include "Sample/Uniform.h"

template <typename T>
class SampleRandomizedSequence : public SamplingMethod
{

public:

    SampleRandomizedSequence(T lower, T upper, double step)
        : unif_m(0, size_t((upper - lower) / step))
        , lower_m(lower)
        , step_m(step)
    { }

    SampleRandomizedSequence(T lower, T upper, double step, size_t seed)
        : unif_m(0, size_t((upper - lower) / step), seed)
        , lower_m(lower)
        , step_m(step)
    { }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        size_t idx = unif_m.getNext();
        ind->genes[i] = static_cast<T>(lower_m + idx * step_m);
    }

    void allocate(const CmdArguments_t& args, const Comm::Bundle_t& comm) {
        unif_m.allocate(args, comm);
    }

private:
    std::vector<T> points_m;
    Uniform<size_t> unif_m;
    T lower_m;
    double step_m;
};

#endif
