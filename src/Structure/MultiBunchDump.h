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

    /** \brief
     *  delete the last 'numberOfLines' lines of the statistics file
     */
    unsigned int rewindToSpos(double maxSpos);

private:
    bool isFirst_m;
};


inline
unsigned int MultiBunchDump::rewindToSpos(double maxSPos) {
    if (Ippl::myNode() == 0) {
        return Util::rewindLinesSDDS(this->fname_m, maxSPos);
    }
    return 0;
}

#endif
