#include "opal.h"

Ippl *ippl;
Inform *gmsg;
Inform *gmsg2all;

int run_opal(char *arg[], std::string inputfile, int restartStep, MPI_Comm comm) {

    int narg = 5, remove = 1;
    if(!ippl)
        ippl = new Ippl(narg, arg, remove, comm);
    //Ippl *aippl = new Ippl(narg, arg, remove, comm);
    //ippl = aippl;

    gmsg = new  Inform("OPAL ");
    gmsg2all = new  Inform("OPAL", INFORM_ALL_NODES);
    *gmsg << "running simulation " << inputfile << " in " << getcwd() << endl;

    OpalData *OPAL = OpalData::getInstance();
    Configure::configure();
    OPAL->storeInputFn(inputfile);

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
    delete parser;
    delete gmsg2all;
    delete gmsg;

    //FIXME: strange side effects
    //ippl = 0;
    //delete aippl;

    //XXX: seems like Ippl is always returning the same instance after the
    //     initial instantiation.
    //delete ippl;

    return 0;
}
