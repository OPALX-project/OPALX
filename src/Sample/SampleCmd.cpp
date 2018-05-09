#include "Sample/SampleCmd.h"
#include "Sample/SVar.h"
#include "Sample/Sampler.h"
// #include "Sample/Objective.h"
// #include "Sample/Constraint.h"
#include "Optimize/OpalSimulation.h"

#include "Attributes/Attributes.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"

//#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"
#include "Utility/IpplTimings.h"
#include "Track/Track.h"

#include "Sample/SPilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"
#include "Util/OpalInputFileParser.h"

#include "Comm/CommSplitter.h"
#include "Comm/Topology/NoCommTopology.h"
#include "Comm/Splitter/ManyMasterSplit.h"
#include "Comm/MasterGraph/SocialNetworkGraph.h"

#include "Expression/Parser/function.hpp"

#include <boost/filesystem.hpp>

extern Inform *gmsg;

namespace {
    enum {
        INPUT,
        OUTPUT,
        OUTDIR,
        SVARS,
        NUMMASTERS,
        NUMCOWORKERS,
        N,
        SIMTMPDIR,
        TEMPLATEDIR,
        FIELDMAPDIR,
        SIZE
    };
}

SampleCmd::SampleCmd():
    Action(SIZE, "SAMPLE",
           "The \"SAMPLE\" command initiates sampling.") {
    itsAttr[INPUT] = Attributes::makeString
        ("INPUT", "Path to input file");
    itsAttr[OUTPUT] = Attributes::makeString
        ("OUTPUT", "Name used in output file sample");
    itsAttr[OUTDIR] = Attributes::makeString
        ("OUTDIR", "Name of directory used to store sample output files");
    itsAttr[SVARS] = Attributes::makeStringArray
        ("SVARS", "List of sampling variables to be used");
    itsAttr[NUMMASTERS] = Attributes::makeReal
        ("NUM_MASTERS", "Number of master nodes");
    itsAttr[NUMCOWORKERS] = Attributes::makeReal
        ("NUM_COWORKERS", "Number processors per worker");
    itsAttr[N] = Attributes::makeReal
        ("N", "Number of samples");
    itsAttr[SIMTMPDIR] = Attributes::makeString
        ("SIMTMPDIR", "Directory where simulations are run");
    itsAttr[TEMPLATEDIR] = Attributes::makeString
        ("TEMPLATEDIR", "Directory where templates are stored");
    itsAttr[FIELDMAPDIR] = Attributes::makeString
        ("FIELDMAPDIR", "Directory where field maps are stored");

    registerOwnership(AttributeHandler::COMMAND);
}

SampleCmd::SampleCmd(const std::string &name, SampleCmd *parent):
    Action(name, parent)
{ }

SampleCmd::~SampleCmd()
{ }

SampleCmd *SampleCmd::clone(const std::string &name) {
    return new SampleCmd(name, this);
}

void SampleCmd::execute() {
    
    namespace fs = boost::filesystem;

    auto opal = OpalData::getInstance();
    fs::path inputfile(Attributes::getString(itsAttr[INPUT]));

    std::vector<std::string> dvarsstr = Attributes::getStringArray(itsAttr[SVARS]);
    DVarContainer_t dvars;
    
    
    // OPT-Pilot wants to have at least 1 objective --> create a dummy one
    Expressions::Named_t objectives;
    objectives.insert(Expressions::SingleNamed_t(
        "dummy", new Expressions::Expr_t("dummy") ));
    
    // constrains --> type of sampling
    Expressions::Named_t constraints;
    

    // Setup/Configuration
    //////////////////////////////////////////////////////////////////////////
    typedef OpalInputFileParser Input_t;
    typedef OpalSimulation Sim_t;

    typedef CommSplitter< ManyMasterSplit< NoCommTopology > > Comm_t;
    typedef SocialNetworkGraph< NoCommTopology > SolPropagationGraph_t;

    typedef SPilot<Input_t, Sampler, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;
    
    //////////////////////////////////////////////////////////////////////////

    std::vector<std::string> arguments(opal->getArguments());
    std::vector<char*> argv;
    std::map<unsigned int, std::string> argumentMapper({
            {INPUT, "inputfile"},
            {OUTPUT, "outfile"},
            {OUTDIR, "outdir"},
            {NUMMASTERS, "num-masters"},
            {NUMCOWORKERS, "num-coworkers"},
            {N, "nsamples"}
        });

    auto it = argumentMapper.end();
    for (unsigned int i = 0; i < SIZE; ++ i) {
        if ((it = argumentMapper.find(i)) != argumentMapper.end()) {
            std::string type = itsAttr[i].getType();
            if (type == "string") {
                if (Attributes::getString(itsAttr[i]) != "") {
                    std::string argument = "--" + (*it).second + "=" + Attributes::getString(itsAttr[i]);
                    arguments.push_back(argument);
                }
            } else if (type == "real") {
                if (itsAttr[i]) {
                    std::string val = std::to_string (Attributes::getReal(itsAttr[i]));
                    size_t last = val.find_last_not_of('0');
                    if (val[last] != '.') ++ last;
                    val.erase (last, std::string::npos );
                    std::string argument = "--" + (*it).second + "=" + val;
                    arguments.push_back(argument);
                }
            } else if (type == "logical") {
                if (itsAttr[i]) {
                    std::string argument = "--" + (*it).second + "=" + std::to_string(Attributes::getBool(itsAttr[i]));
                    arguments.push_back(argument);
                }
            }
        }
    }
    if (Attributes::getString(itsAttr[INPUT]) == "") {
        throw OpalException("SampleCmd::execute",
                            "The argument INPUT has to be provided");
    }
    if (Attributes::getReal(itsAttr[N]) <= 0) {
        throw OpalException("SampleCmd::execute",
                            "The argument N has to be provided");
    }

    if (Attributes::getString(itsAttr[SIMTMPDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[SIMTMPDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        *gmsg << dir.native() << endl;
        if (!fs::exists(dir)) {
            fs::create_directory(dir);
        }
        std::string argument = "--simtmpdir=" + dir.native();
        arguments.push_back(argument);
    }

    if (Attributes::getString(itsAttr[TEMPLATEDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[TEMPLATEDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        std::string argument = "--templates=" + dir.native();
        arguments.push_back(argument);
    }

    if (Attributes::getString(itsAttr[FIELDMAPDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[FIELDMAPDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        setenv("FIELDMAPS", dir.c_str(), 1);
    }
    
    std::string argument = "--sample=1";
    arguments.push_back(argument);
    

    *gmsg << endl;
    for (size_t i = 0; i < arguments.size(); ++ i) {
        argv.push_back(const_cast<char*>(arguments[i].c_str()));
        *gmsg << arguments[i] << " ";
    }
    *gmsg << endl;

    for (const std::string &name: dvarsstr) {
        Object *obj = opal->find(name);
        SVar* dvar = static_cast<SVar*>(obj);
        std::string var = dvar->getVariable();
        double lowerbound = dvar->getLowerBound();
        double upperbound = dvar->getUpperBound();
        std::string type  = dvar->getType();

        DVar_t tmp = boost::make_tuple(var, lowerbound, upperbound);
        dvars.insert(namedDVar_t(name, tmp));
        
        // type of sampling
        constraints.insert(Expressions::SingleNamed_t(
            name, new Expressions::Expr_t(type) ));
    }

    Inform *origGmsg = gmsg;
    gmsg = 0;
    stashEnvironment();
    try {
        CmdArguments_t args(new CmdArguments(argv.size(), &argv[0]));

        boost::shared_ptr<Comm_t>  comm(new Comm_t(args, MPI_COMM_WORLD));
        boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm, dvars, objectives, constraints));

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -100);
    }
    popEnvironment();
    gmsg = origGmsg;
}

void SampleCmd::stashEnvironment() {
    Ippl::stash();
    IpplTimings::stash();
    Track::stash();
    OpalData::stashInstance();
}

void SampleCmd::popEnvironment() {
    Ippl::pop();
    IpplTimings::pop();
    OpalData::popInstance();
    Track::pop();
}