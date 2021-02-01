//
// Class Option
//   The OPTION command.
//   The user interface allowing setting of OPAL options.
//   The actual option flags are contained in namespace Options.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "BasicActions/Option.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Parser/FileStream.h"
#include "Utilities/Options.h"
#include "Utilities/OptionTypes.h"
#include "Utilities/ClassicRandom.h"
#include "Utility/IpplInfo.h"

#include "Utilities/OpalException.h"

#include "Utility/IpplMemoryUsage.h"

#include <ctime>
#include <iostream>
#include <limits>
#include <cstddef>

extern Inform *gmsg;

using namespace Options;

std::string DumpFrameToString(DumpFrame df);

namespace {
    // The attributes of class Option.
    enum {
        ECHO,
        INFO,
        TRACE,
        WARN,
        SEED,
        TELL,
        PSDUMPFREQ,
        STATDUMPFREQ,
        PSDUMPEACHTURN,
        PSDUMPFRAME,
        SPTDUMPFREQ,
        REPARTFREQ,
        REBINFREQ,
        SCSOLVEFREQ,
        MTSSUBSTEPS,
        REMOTEPARTDEL,
        RHODUMP,
        EBDUMP,
        CSRDUMP,
        AUTOPHASE,
        SURFDUMPFREQ,
        NUMBLOCKS,
         RECYCLEBLOCKS,
        NLHS,
        CZERO,
        RNGTYPE,
        ENABLEHDF5,
        ASCIIDUMP,
        BOUNDPDESTROYFQ,
        BEAMHALOBOUNDARY,
        CLOTUNEONLY,
        IDEALIZED,
        LOGBENDTRAJECTORY,
        VERSION,
#ifdef ENABLE_AMR
        AMR,
        AMR_YT_DUMP_FREQ,
        AMR_REGRID_FREQ,
#endif
        MEMORYDUMP,
        HALOSHIFT,
        DELPARTFREQ,
        MINBINEMITTED,
        MINSTEPFORREBIN,
        SIZE
    };
}


Option::Option():
    Action(SIZE, "OPTION",
           "The \"OPTION\" statement defines OPAL execution options.")
    {

    itsAttr[ECHO] = Attributes::makeBool
                    ("ECHO", "If true, give echo of input", echo);

    itsAttr[INFO] = Attributes::makeBool
                    ("INFO", "If true, print information messages", info);

    itsAttr[TRACE] = Attributes::makeBool
                     ("TRACE", "If true, print execution trace", mtrace);

    itsAttr[WARN] = Attributes::makeBool
                    ("WARN", "If true, print warning messages", warn);

    itsAttr[SEED] = Attributes::makeReal
                    ("SEED", "The seed for the random generator, -1 will use time(0) as seed ");

    itsAttr[TELL] = Attributes::makeBool
                    ("TELL", "If true, print the current settings. "
                     "Must be the last option in the inputfile in "
                     "order to render correct results", false);

    itsAttr[PSDUMPFREQ] = Attributes::makeReal
                          ("PSDUMPFREQ", "The frequency to dump the phase space, "
                           "i.e.dump data when step%psDumpFreq==0, its default value is 10.",
                           psDumpFreq);

    itsAttr[STATDUMPFREQ] = Attributes::makeReal
                            ("STATDUMPFREQ", "The frequency to dump statistical data "
                             "(e.g. RMS beam quantities), i.e. dump data when step%statDumpFreq == 0, "
                             "its default value is 10.", statDumpFreq);

    itsAttr[PSDUMPEACHTURN] = Attributes::makeBool
                              ("PSDUMPEACHTURN", "If true, dump phase space after each "
                               "turn ,only aviable for OPAL-cycl, its default value is false",
                               psDumpEachTurn);

    itsAttr[SCSOLVEFREQ] = Attributes::makeReal
                           ("SCSOLVEFREQ", "The frequency to solve space charge fields. its default value is 1");

    itsAttr[MTSSUBSTEPS] = Attributes::makeReal("MTSSUBSTEPS", "How many small timesteps "
                                                "are inside the large timestep used in multiple "
                                                "time stepping (MTS) integrator");

    itsAttr[REMOTEPARTDEL] = Attributes::makeReal
      ("REMOTEPARTDEL", "Artifically delete the remote particle if its distance "
       "to the beam mass is larger than REMOTEPARTDEL times of the beam rms size, "
       "its default values is 0 (no delete) ",0.0);

    itsAttr[PSDUMPFRAME] = Attributes::makeUpperCaseString
                                ("PSDUMPFRAME", "Controls the frame of phase space dump in "
                                 "stat file and h5 file. If 'GLOBAL' OPAL will dump in the "
                                 "lab (global) Cartesian frame; if 'BUNCH_MEAN' OPAL will "
                                 "dump in the local Cartesian frame of the beam mean; "
                                 "if 'REFERENCE'  OPAL will dump in the local Cartesian "
                                 "frame of the reference particle 0. Only aviable for "
                                 "OPAL-cycl, its default value is 'GLOBAL'");

    itsAttr[SPTDUMPFREQ] = Attributes::makeReal
                           ("SPTDUMPFREQ", "The frequency to dump single "
                            "particle trajectory of particles with ID = 0 & 1, "
                            "its default value is 1. ");

    itsAttr[REPARTFREQ] = Attributes::makeReal
                          ("REPARTFREQ", "The frequency to do particles repartition "
                           "for better load balance between nodes, its "
                           "default value is " + std::to_string(repartFreq) + ".", repartFreq);

    itsAttr[MINBINEMITTED] = Attributes::makeReal
                             ("MINBINEMITTED", "The number of bins that have to be emitted before the bins are squashed into "
                              "a single bin; the default value is " + std::to_string(minBinEmitted) + ".", minBinEmitted);

    itsAttr[MINSTEPFORREBIN] = Attributes::makeReal
                            ("MINSTEPFORREBIN", "The number of steps into the simulation before the bins are squashed into "
                             "a single bin; the default value is " + std::to_string(minStepForRebin) + ".", minStepForRebin);

    itsAttr[REBINFREQ] = Attributes::makeReal
                         ("REBINFREQ", "The frequency to reset energy bin ID for "
                          "all particles, its default value is 100.", rebinFreq);

    itsAttr[RHODUMP] = Attributes::makeBool
                       ("RHODUMP", "If true, in addition to the phase "
                        "space the scalar rho field is also dumped (H5Block)", rhoDump);

    itsAttr[EBDUMP] = Attributes::makeBool
                       ("EBDUMP", "If true, in addition to the phase space the "
                        "E and B field at each particle is also dumped into the H5 file)", ebDump);

    itsAttr[CSRDUMP] = Attributes::makeBool
                       ("CSRDUMP", "If true, the csr E field, line density "
                        "and the line density derivative is dumped into the "
                        "data directory)", csrDump);

    itsAttr[AUTOPHASE] = Attributes::makeReal
                         ("AUTOPHASE", "If greater than zero OPAL is scanning "
                          "the phases of each rf structure in order to get maximum "
                          "acceleration. Defines the number of refinements of the "
                          "search range", autoPhase);

    itsAttr[SURFDUMPFREQ] =  Attributes::makeReal
                             ("SURFDUMPFREQ", "The frequency to dump surface-particle "
                              "interaction data, its default value is -1 (no dump).",
                              surfDumpFreq);

    itsAttr[CZERO] =  Attributes::makeBool
                      ("CZERO", "If set to true a symmetric distribution is "
                       "created -> centroid == 0.0 ", cZero);

    itsAttr[RNGTYPE] =  Attributes::makeUpperCaseString
                        ("RNGTYPE", "RANDOM (default), Quasi-random number "
                         "generators: HALTON, SOBOL, NIEDERREITER (Gsl ref manual 18.5)", rngtype);


    itsAttr[CLOTUNEONLY] =  Attributes::makeBool
                                   ("CLOTUNEONLY", "If set to true stop after "
                                    "CLO and tune calculation ", cloTuneOnly);

    itsAttr[NUMBLOCKS] = Attributes::makeReal
                          ("NUMBLOCKS", "Maximum number of vectors in the Krylov "
                           "space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");
    itsAttr[RECYCLEBLOCKS] = Attributes::makeReal
                          ("RECYCLEBLOCKS", "Number of vectors in the recycle "
                           "space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");
    itsAttr[NLHS] = Attributes::makeReal
                          ("NLHS", "Number of stored old solutions for extrapolating "
                           "the new starting vector. Default value is 1 and just the last solution is used.");

    itsAttr[ENABLEHDF5] = Attributes::makeBool
        ("ENABLEHDF5", "If true, HDF5 actions are enabled", enableHDF5);

    itsAttr[ASCIIDUMP] = Attributes::makeBool
        ("ASCIIDUMP", "If true, some of the elements dump in ASCII instead of HDF5", false);

    itsAttr[BOUNDPDESTROYFQ] = Attributes::makeReal
      ("BOUNDPDESTROYFQ", "The frequency to do boundp_destroy to delete lost particles. Default 10",10.0);

    itsAttr[BEAMHALOBOUNDARY] = Attributes::makeReal
      ("BEAMHALOBOUNDARY", "Defines in terms of sigma where the halo starts. Default 0.0",0.0);

    itsAttr[IDEALIZED] = Attributes::makeBool
        ("IDEALIZED", "Using the hard edge model for the calculation of path length. Default: false", false);

    itsAttr[LOGBENDTRAJECTORY] = Attributes::makeBool
        ("LOGBENDTRAJECTORY", "Writing the trajectory of every bend to disk. Default: false", false);

    itsAttr[VERSION] = Attributes::makeReal
        ("VERSION", "Version of OPAL for which input file was written", 10000);

#ifdef ENABLE_AMR
    itsAttr[AMR] = Attributes::makeBool
        ("AMR", "Use adaptive mesh refinement.", amr);

    itsAttr[AMR_YT_DUMP_FREQ] = Attributes::makeReal("AMR_YT_DUMP_FREQ",
                                                     "The frequency to dump grid "
                                                     "and particle data "
                                                     "(default: 10)", amrYtDumpFreq);

    itsAttr[AMR_REGRID_FREQ] = Attributes::makeReal("AMR_REGRID_FREQ",
                                                    "The frequency to perform a regrid "
                                                    "in multi-bunch mode (default: 10)",
                                                    amrRegridFreq);
#endif
    itsAttr[MEMORYDUMP] = Attributes::makeBool
        ("MEMORYDUMP", "If true, write memory to SDDS file", memoryDump);

    itsAttr[HALOSHIFT] = Attributes::makeReal
        ("HALOSHIFT", "Constant parameter to shift halo value (default: 0.0)", haloShift);

    itsAttr[DELPARTFREQ] = Attributes::makeReal
        ("DELPARTFREQ", "The frequency to delete particles, i.e. delete when step%delPartFreq == 0. Default: 1", delPartFreq);

    registerOwnership(AttributeHandler::STATEMENT);

    FileStream::setEcho(echo);
}


Option::Option(const std::string &name, Option *parent):
    Action(name, parent) {
    Attributes::setBool(itsAttr[ECHO],       echo);
    Attributes::setBool(itsAttr[INFO],       info);
    Attributes::setBool(itsAttr[TRACE],      mtrace);
    Attributes::setBool(itsAttr[WARN],       warn);
    Attributes::setReal(itsAttr[SEED],       seed);
    Attributes::setReal(itsAttr[PSDUMPFREQ], psDumpFreq);
    Attributes::setReal(itsAttr[STATDUMPFREQ], statDumpFreq);
    Attributes::setBool(itsAttr[PSDUMPEACHTURN], psDumpEachTurn);
    Attributes::setUpperCaseString(itsAttr[PSDUMPFRAME], DumpFrameToString(psDumpFrame));
    Attributes::setReal(itsAttr[SPTDUMPFREQ], sptDumpFreq);
    Attributes::setReal(itsAttr[SCSOLVEFREQ], scSolveFreq);
    Attributes::setReal(itsAttr[MTSSUBSTEPS], mtsSubsteps);
    Attributes::setReal(itsAttr[REMOTEPARTDEL], remotePartDel);
    Attributes::setReal(itsAttr[REPARTFREQ], repartFreq);
    Attributes::setReal(itsAttr[MINBINEMITTED], minBinEmitted);
    Attributes::setReal(itsAttr[MINSTEPFORREBIN], minStepForRebin);
    Attributes::setReal(itsAttr[REBINFREQ], rebinFreq);
    Attributes::setBool(itsAttr[RHODUMP], rhoDump);
    Attributes::setBool(itsAttr[EBDUMP], ebDump);
    Attributes::setBool(itsAttr[CSRDUMP], csrDump);
    Attributes::setReal(itsAttr[AUTOPHASE], autoPhase);
    Attributes::setReal(itsAttr[SURFDUMPFREQ], surfDumpFreq);
    Attributes::setBool(itsAttr[CZERO], cZero);
    Attributes::setBool(itsAttr[CLOTUNEONLY], cloTuneOnly);
    Attributes::setUpperCaseString(itsAttr[RNGTYPE], std::string(rngtype));
    Attributes::setReal(itsAttr[NUMBLOCKS], numBlocks);
    Attributes::setReal(itsAttr[RECYCLEBLOCKS], recycleBlocks);
    Attributes::setReal(itsAttr[NLHS], nLHS);
    Attributes::setBool(itsAttr[ENABLEHDF5], enableHDF5);
    Attributes::setBool(itsAttr[ASCIIDUMP], asciidump);
    Attributes::setReal(itsAttr[BOUNDPDESTROYFQ], boundpDestroyFreq);
    Attributes::setReal(itsAttr[BEAMHALOBOUNDARY], beamHaloBoundary);
    Attributes::setBool(itsAttr[IDEALIZED], idealized);
    Attributes::setBool(itsAttr[LOGBENDTRAJECTORY], writeBendTrajectories);
    Attributes::setReal(itsAttr[VERSION], version);
#ifdef ENABLE_AMR
    Attributes::setBool(itsAttr[AMR], amr);
    Attributes::setReal(itsAttr[AMR_YT_DUMP_FREQ], amrYtDumpFreq);
    Attributes::setReal(itsAttr[AMR_REGRID_FREQ], amrRegridFreq);
#endif
    Attributes::setBool(itsAttr[MEMORYDUMP], memoryDump);
    Attributes::setReal(itsAttr[HALOSHIFT], haloShift);
    Attributes::setReal(itsAttr[DELPARTFREQ], delPartFreq);
}


Option::~Option()
{}


Option *Option::clone(const std::string &name) {
    return new Option(name, this);
}


void Option::execute() {
    // Store the option flags.
    echo      = Attributes::getBool(itsAttr[ECHO]);
    info      = Attributes::getBool(itsAttr[INFO]);
    mtrace     = Attributes::getBool(itsAttr[TRACE]);
    warn      = Attributes::getBool(itsAttr[WARN]);
    psDumpEachTurn =   Attributes::getBool(itsAttr[PSDUMPEACHTURN]);
    rhoDump = Attributes::getBool(itsAttr[RHODUMP]);
    ebDump = Attributes::getBool(itsAttr[EBDUMP]);
    csrDump = Attributes::getBool(itsAttr[CSRDUMP]);
    enableHDF5 = Attributes::getBool(itsAttr[ENABLEHDF5]);
    version = Attributes::getReal(itsAttr[VERSION]);
#ifdef ENABLE_AMR
    amr = Attributes::getBool(itsAttr[AMR]);
    amrYtDumpFreq = int(Attributes::getReal(itsAttr[AMR_YT_DUMP_FREQ]));

    if ( amrYtDumpFreq < 1 ) {
        amrYtDumpFreq = std::numeric_limits<int>::max();
    }

    amrRegridFreq = int(Attributes::getReal(itsAttr[AMR_REGRID_FREQ]));
    amrRegridFreq = ( amrRegridFreq < 1 ) ? 1 : amrRegridFreq;
#endif
    memoryDump  = Attributes::getBool(itsAttr[MEMORYDUMP]);
    haloShift   = Attributes::getReal(itsAttr[HALOSHIFT]);
    delPartFreq = Attributes::getReal(itsAttr[DELPARTFREQ]);

    if ( memoryDump ) {
        IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance(
                IpplMemoryUsage::Unit::GB, false);
        memory->sample();
    }

    seed = Attributes::getReal(itsAttr[SEED]);

    /// note: rangen is used only for the random number generator in the OPAL language
    ///       not for the distributions

    if (Options::seed == -1)
      rangen.init55(time(0));
    else
      rangen.init55(seed);


    IpplInfo::Info->on(info);
    IpplInfo::Warn->on(warn);

    handlePsDumpFrame(Attributes::getString(itsAttr[PSDUMPFRAME]));

    if(itsAttr[ASCIIDUMP]) {
        asciidump = Attributes::getBool(itsAttr[ASCIIDUMP]);
    }

    /// note: rangen is used only for the random number generator in the OPAL language
    ///       not for the distributions
    if(itsAttr[SEED]) {
        seed = int(Attributes::getReal(itsAttr[SEED]));
        if (seed == -1)
            rangen.init55(time(0));
        else
            rangen.init55(seed);
    }

    if(itsAttr[PSDUMPFREQ]) {
        psDumpFreq = int(Attributes::getReal(itsAttr[PSDUMPFREQ]));
        if (psDumpFreq==0)
            psDumpFreq = std::numeric_limits<int>::max();
    }

    if(itsAttr[STATDUMPFREQ]) {
        statDumpFreq = int(Attributes::getReal(itsAttr[STATDUMPFREQ]));
        if (statDumpFreq==0)
            statDumpFreq = std::numeric_limits<int>::max();
    }

    if(itsAttr[SPTDUMPFREQ]) {
        sptDumpFreq = int(Attributes::getReal(itsAttr[SPTDUMPFREQ]));
        if (sptDumpFreq==0)
            sptDumpFreq = std::numeric_limits<int>::max();
    }


    if(itsAttr[SCSOLVEFREQ]) {
        scSolveFreq = int(Attributes::getReal(itsAttr[SCSOLVEFREQ]));
        scSolveFreq = ( scSolveFreq < 1 ) ? 1 : scSolveFreq;
    }


    if(itsAttr[MTSSUBSTEPS]) {
        mtsSubsteps = int(Attributes::getReal(itsAttr[MTSSUBSTEPS]));
    }


    if(itsAttr[REMOTEPARTDEL]) {
        remotePartDel = Attributes::getReal(itsAttr[REMOTEPARTDEL]);
    }

    if(itsAttr[REPARTFREQ]) {
        repartFreq = int(Attributes::getReal(itsAttr[REPARTFREQ]));
    }

    if (itsAttr[MINBINEMITTED]) {
        minBinEmitted = int(Attributes::getReal(itsAttr[MINBINEMITTED]));
    }

    if (itsAttr[MINSTEPFORREBIN]) {
        minStepForRebin = int(Attributes::getReal(itsAttr[MINSTEPFORREBIN]));
    }

    if(itsAttr[REBINFREQ]) {
        rebinFreq = int(Attributes::getReal(itsAttr[REBINFREQ]));
    }

    if(itsAttr[AUTOPHASE]) {
        autoPhase = int(Attributes::getReal(itsAttr[AUTOPHASE]));
    }
    if(itsAttr[SURFDUMPFREQ]) {
        surfDumpFreq = int(Attributes::getReal(itsAttr[SURFDUMPFREQ]));
    }
    if(itsAttr[NUMBLOCKS]) {
        numBlocks = int(Attributes::getReal(itsAttr[NUMBLOCKS]));
    }
    if(itsAttr[RECYCLEBLOCKS]) {
        recycleBlocks = int(Attributes::getReal(itsAttr[RECYCLEBLOCKS]));
    }
    if(itsAttr[NLHS]) {
        nLHS = int(Attributes::getReal(itsAttr[NLHS]));
    }

    if(itsAttr[CZERO]) {
        cZero = bool(Attributes::getBool(itsAttr[CZERO]));
    }

    if(itsAttr[RNGTYPE]) {
        rngtype = std::string(Attributes::getString(itsAttr[RNGTYPE]));
    } else {
        rngtype = std::string("RANDOM");
    }

    if(itsAttr[BEAMHALOBOUNDARY]) {
        beamHaloBoundary  = Attributes::getReal(itsAttr[BEAMHALOBOUNDARY]);
    }
    else {
        beamHaloBoundary = 0;
    }

    idealized = Attributes::getBool(itsAttr[IDEALIZED]);

    writeBendTrajectories = Attributes::getBool(itsAttr[LOGBENDTRAJECTORY]);

    if(itsAttr[CLOTUNEONLY]) {
        cloTuneOnly = bool(Attributes::getBool(itsAttr[CLOTUNEONLY]));
    } else {
        cloTuneOnly = false;
    }

    // Set message flags.
    FileStream::setEcho(echo);

    if(Attributes::getBool(itsAttr[TELL])) {
        *gmsg << "\nCurrent settings of options:\n" << *this << endl;
    }

    Option* main = dynamic_cast<Option*>(OpalData::getInstance()->find("OPTION"));
    if (main) {
        main->update(itsAttr);
    }
}

void Option::handlePsDumpFrame(const std::string &dumpFrame) {
    if (dumpFrame == "GLOBAL") {
        psDumpFrame = GLOBAL;
    } else if (dumpFrame == "BUNCH_MEAN") {
        psDumpFrame = BUNCH_MEAN;
    } else if (dumpFrame == "REFERENCE") {
        psDumpFrame = REFERENCE;
    } else {
        std::string msg = std::string("Did not recognise PSDUMPFRAME '")+\
                    dumpFrame+std::string("'. It should be one of 'GLOBAL',")+\
                    std::string(" 'BUNCH_MEAN' or 'REFERENCE'");
        throw OpalException("Option::handlePsDumpFrame", msg);
    }
}

std::string DumpFrameToString(DumpFrame df) {
    switch (df) {
    case BUNCH_MEAN:
        return std::string("BUNCH_MEAN");
    case REFERENCE:
        return std::string("REFERENCE");
    case GLOBAL:
    default:
        return std::string("GLOBAL");
    }
}

void Option::update(const std::vector<Attribute>& othersAttributes) {
    for (int i = 0; i < SIZE; ++ i) {
        itsAttr[i] = othersAttributes[i];
    }
}