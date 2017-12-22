#include "opal.h"

extern Ippl *ippl;
extern Inform *gmsg;


#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Utilities/OpalException.h"
#include "Fields/Fieldmap.h"
#include "Structure/IpplInfoWrapper.h"

#include "OPALconfig.h"


int run_opal(char *arg[], std::string inputfile, int restartStep, MPI_Comm comm) {

    std::string::size_type startExtension    = inputfile.find_last_of('.');
    std::string::size_type startRelativePath = inputfile.find_last_of('/');
    std::string relativePath("");
    if (startRelativePath != std::string::npos) {
        relativePath = inputfile.substr(0, startRelativePath + 1);
    }
    std::string outputFileName = inputfile.substr(0,startExtension) + ".out";
    std::ofstream output(outputFileName.c_str());

    MPI_Barrier(comm);

    // int narg = 5, remove = 1;
    IpplInfoWrapper *newippl = new IpplInfoWrapper(inputfile, comm);//narg, arg, remove, comm);
    gmsg = new Inform("OPAL ", output);
    IpplInfo::Info->setDestination(output);
    IpplInfo::Error->setDestination(output);
    IpplInfo::Warn->setDestination(output);

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

    IpplInfo::Info->setDestination(std::cout);
    IpplInfo::Error->setDestination(std::cout);
    IpplInfo::Warn->setDestination(std::cout);

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
    delete newippl;

    output.close();
    return 0;
}