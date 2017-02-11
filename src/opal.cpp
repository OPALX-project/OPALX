#include "opal.h"

Ippl *ippl;
Inform *gmsg;


#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Utilities/OpalException.h"
#include "Fields/Fieldmap.h"

#include "OPALconfig.h"


int run_opal(char *arg[], std::string inputfile, int restartStep, MPI_Comm comm) {

    MPI_Barrier(comm);

    int narg = 5, remove = 1;

    if(!ippl)
        ippl = new Ippl(narg, arg, remove, comm);
    //Ippl *aippl = new Ippl(narg, arg, remove, comm);
    //ippl = aippl;

    gmsg = new Inform("OPAL ");

    OpalData *opal = OpalData::getInstance();
    Configure::configure();
    opal->storeInputFn(inputfile);

    //FIXME
    if(restartStep > 0) throw new OpalException("run_opal", "Restart not implemented yet!");

    // FileStream is a RCObject
    FileStream *is = 0;
    try {
        is = new FileStream(inputfile);
    } catch(...) {
        is = 0;
        throw new OpalException("run_opal", "Could not open inputfile: " + inputfile);
    }

    // run simulation
    OpalParser *parser = new OpalParser();
    if(is) parser->run(is);

    Ippl::Comm->barrier();

    // cleanup
    //OPAL->reset();
    OpalData::deleteInstance();
    Fieldmap::clearDictionary();
    delete parser;
    delete gmsg;

    //FIXME: strange side effects
    //ippl = 0;
    //delete aippl;

    //XXX: seems like Ippl is always returning the same instance after the
    //     initial instantiation.
    //delete ippl;

    return 0;
}