#ifndef OPAL_H5PARTWRAPPERFORPT_H
#define OPAL_H5PARTWRAPPERFORPT_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  H5PartWrapperForPT: a class that manages all calls to H5Part for the Parallel-T tracker
*/

#include "Structure/H5PartWrapper.h"

#include "Algorithms/PBunchDefs.h"
#include "Utilities/OpalException.h"

#include "H5hut.h"

class H5PartWrapperForPT: public H5PartWrapper {
public:
    H5PartWrapperForPT(const std::string &fileName, h5_int32_t flags = H5_O_WRONLY);
    H5PartWrapperForPT(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags = H5_O_RDWR);
    virtual ~H5PartWrapperForPT();

    virtual void readHeader();
    virtual void readStep(PartBunch&, h5_ssize_t firstParticle, h5_ssize_t lastParticle);

    virtual void writeHeader();
    virtual void writeStep(PartBunch&, const std::map<std::string, double> &additionalStepAttributes);

    virtual bool predecessorIsSameFlavour() const;

private:
    void readStepHeader(PartBunch&);
    void readStepData(PartBunch&, h5_ssize_t, h5_ssize_t);

    void writeStepHeader(PartBunch&, const std::map<std::string, double> &);
    void writeStepData(PartBunch&);
};

inline
bool H5PartWrapperForPT::predecessorIsSameFlavour() const {
    return (predecessorOPALFlavour_m == "opal-t");
}

#endif //OPAL_H5PARTWRAPPERFORPT_H