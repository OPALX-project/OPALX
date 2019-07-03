#ifndef OPAL_MULTI_BUNCH_DUMP_H
#define OPAL_MULTI_BUNCH_DUMP_H

#include "SDDSWriter.h"

#include "Algorithms/MultiBunchHandler.h"


class MultiBunchDump : public SDDSWriter {

public:
    typedef MultiBunchHandler::beaminfo_t beaminfo_t;

    MultiBunchDump(const std::string& fname, bool restart);

    void fillHeader();

    void write(PartBunchBase<double, 3>* beam, const beaminfo_t& binfo);

private:
    bool isFirst_m;
};

#endif
