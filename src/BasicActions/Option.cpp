// ------------------------------------------------------------------------
// $RCSfile: Option.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Option
//   The class for the OPAL OPTION command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Option.h"
#include "Attributes/Attributes.h"
#include "Parser/FileStream.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/Random.h"
#include <ctime>
#include <iostream>

extern Inform *gmsg;

using namespace Options;

// Class Option
// ------------------------------------------------------------------------

namespace {
    // The attributes of class Option.
    enum {
        ECHO,
        INFO,
        TRACE,
        VERIFY,
        WARN,
        SEED,
        TELL,
        PSDUMPFREQ,
        STATDUMPFREQ,
        PSDUMPEACHTURN,
        PSDUMPLOCALFRAME,
        SPTDUMPFREQ,
        REPARTFREQ,
        REBINFREQ,
        SCSOLVEFREQ,
        MTSSUBSTEPS,
	REMOTEPARTDEL,
        SCAN,
        RHODUMP,
        EBDUMP,
	CSRDUMP,
        AUTOPHASE,
        PPDEBUG,
        SURFDUMPFREQ,
        NUMBLOCKS,
        RECYCLEBLOCKS,
        NLHS,
        CZERO,
        RNGTYPE,
        SCHOTTKYCORR,
        SCHOTTKYRENO,
        ENABLEHDF5,
        ASCIIDUMP,
        BOUNDPDESTROYFQ,
	BEAMHALOBOUNDARY,
	CLOTUNEONLY,
        SIZE
    };
}


Option::Option():
    Action(SIZE, "OPTION",
           "The \"OPTION\" statement defines OPAL execution options.") {
    itsAttr[ECHO] = Attributes::makeBool
                    ("ECHO", "If true, give echo of input", echo);
    itsAttr[INFO] = Attributes::makeBool
                    ("INFO", "If true, print information messages", info);
    itsAttr[TRACE] = Attributes::makeBool
                     ("TRACE", "If true, print execution trace", mtrace);
    itsAttr[VERIFY] = Attributes::makeBool
                      ("VERIFY", "If true, print warnings about assumptions", verify);
    itsAttr[WARN] = Attributes::makeBool
                    ("WARN", "If true, print warning messages", warn);
    itsAttr[SEED] = Attributes::makeReal
                    ("SEED", "The seed for the random generator, -1 will use time(0) as seed ");
    itsAttr[TELL] = Attributes::makeBool
                    ("TELL", "If true, print the current settings", false);
    itsAttr[PSDUMPFREQ] = Attributes::makeReal
                          ("PSDUMPFREQ", "The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0, its default value is 10.");
    itsAttr[STATDUMPFREQ] = Attributes::makeReal
                            ("STATDUMPFREQ", "The frequency to dump statistical data (e.g. RMS beam quantities), i.e. dump data when step%statDumpFreq == 0, its default value is 10.");
    itsAttr[PSDUMPEACHTURN] = Attributes::makeBool
                              ("PSDUMPEACHTURN", "If true, dump phase space after each turn ,only aviable for OPAL-cycl, its default value is false");
    itsAttr[SCSOLVEFREQ] = Attributes::makeReal
                           ("SCSOLVEFREQ", "The frequency to solve space charge fields. its default value is 1");
    itsAttr[MTSSUBSTEPS] = Attributes::makeReal("MTSSUBSTEPS", "How many small timesteps are inside the large timestep used in multiple time stepping (MTS) integrator");
    itsAttr[REMOTEPARTDEL] = Attributes::makeReal
      ("REMOTEPARTDEL", "Artifically delete the remote particle if its distance to the beam mass is larger than REMOTEPARTDEL times of the beam rms size, its default values is 0 (no delete) ",0.0);
    itsAttr[PSDUMPLOCALFRAME] = Attributes::makeBool
                                ("PSDUMPLOCALFRAME", "If true, in local Cartesian frame, otherwise in global Cartesian frame, only aviable for OPAL-cycl, its default value is false");
    itsAttr[SPTDUMPFREQ] = Attributes::makeReal
                           ("SPTDUMPFREQ", "The frequency to dump single particle trajectory of particles with ID = 0 & 1, its default value is 1. ");
    itsAttr[REPARTFREQ] = Attributes::makeReal
                          ("REPARTFREQ", "The frequency to do particles repartition for better load balance between nodes,its default value is 10. ");
    itsAttr[REBINFREQ] = Attributes::makeReal
                         ("REBINFREQ", "The frequency to reset energy bin ID for all particles, its default value is 100. ");

    itsAttr[SCAN] = Attributes::makeBool
                    ("SCAN", "If true, each new track starts at the begin of the lattics with a new distribution", scan);

    itsAttr[RHODUMP] = Attributes::makeBool
                       ("RHODUMP", "If true, in addition to the phase space the scalar rho field is also dumped (H5Block)", rhoDump);

    itsAttr[EBDUMP] = Attributes::makeBool
                       ("EBDUMP", "If true, in addition to the phase space the E and B field at each particle is also dumped into the H5 file)", ebDump);

    itsAttr[CSRDUMP] = Attributes::makeBool
                       ("CSRDUMP", "If true, the csr E field, line density and the line density derivative is dumped into the data directory)", csrDump);

    itsAttr[AUTOPHASE] = Attributes::makeReal
                         ("AUTOPHASE", "If greater than zero OPAL is scaning the phases of each rf structure in order to get maximum acceleration. Defines the number of refinements of the search range", autoPhase);

    itsAttr[PPDEBUG] = Attributes::makeBool
                       ("PPDEBUG", "If true, use special initial velocity distribution for parallel plate and print special debug output", ppdebug);
    itsAttr[SURFDUMPFREQ] =  Attributes::makeReal
                             ("SURFDUMPFREQ", "The frequency to dump surface-partcle interaction data, its default value is -1 (no dump). ");

    itsAttr[CZERO] =  Attributes::makeBool
                      ("CZERO", "If set to true a symmetric distribution is created -> centroid == 0.0 ", cZero);

    itsAttr[RNGTYPE] =  Attributes::makeString
                        ("RNGTYPE", "RANDOM (default), Quasi-random number gernerators: HALTON, SOBOL, NIEDERREITER (Gsl ref manual 18.5)", rngtype);

    itsAttr[SCHOTTKYCORR] =  Attributes::makeBool
                                   ("SCHOTTKYCORR", "If set to true a Schottky correction to the charge is applied ", schottkyCorrection);

    itsAttr[CLOTUNEONLY] =  Attributes::makeBool
                                   ("CLOTUNEONLY", "If set to true stop after CLO and tune calculation ", cloTuneOnly);

    itsAttr[SCHOTTKYRENO] =  Attributes::makeReal
                                         ("SCHOTTKYRENO", "IF set to a value greater than 0.0 the Schottky correction scan is disabled and the value is used for charge renormalization ", schottkyRennormalization);


    itsAttr[NUMBLOCKS] = Attributes::makeReal
                          ("NUMBLOCKS", "Maximum number of vectors in the Krylov space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");
    itsAttr[RECYCLEBLOCKS] = Attributes::makeReal
                          ("RECYCLEBLOCKS", "Number of vectors in the recycle space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");
    itsAttr[NLHS] = Attributes::makeReal
                          ("NLHS", "Number of stored old solutions for extrapolating the new starting vector. Default value is 1 and just the last solution is used.");

    itsAttr[ENABLEHDF5] = Attributes::makeBool
        ("ENABLEHDF5", "If true, HDF5 actions are enabled", enableHDF5);

    itsAttr[ASCIIDUMP] = Attributes::makeBool
        ("ASCIIDUMP", "If true, some of the elements dump in ASCII instead of HDF5", false);

    itsAttr[BOUNDPDESTROYFQ] = Attributes::makeReal
      ("BOUNDPDESTROYFQ", "The frequency to do boundp_destroy to delete lost particles. Default 10",10.0);

    itsAttr[BEAMHALOBOUNDARY] = Attributes::makeReal
      ("BEAMHALOBOUNDARY", "Defines in therms of sigma where the halo starts Default 0.0",0.0);

    FileStream::setEcho(echo);
    rangen.init55(seed);
}


Option::Option(const std::string &name, Option *parent):
    Action(name, parent) {
    Attributes::setBool(itsAttr[ECHO],       echo);
    Attributes::setBool(itsAttr[INFO],       info);
    Attributes::setBool(itsAttr[TRACE],      mtrace);
    Attributes::setBool(itsAttr[VERIFY],     verify);
    Attributes::setBool(itsAttr[WARN],       warn);
    Attributes::setReal(itsAttr[SEED],       seed);
    Attributes::setReal(itsAttr[PSDUMPFREQ], psDumpFreq);
    Attributes::setReal(itsAttr[STATDUMPFREQ], statDumpFreq);
    Attributes::setBool(itsAttr[PSDUMPEACHTURN], psDumpEachTurn);
    Attributes::setBool(itsAttr[PSDUMPLOCALFRAME], psDumpLocalFrame);
    Attributes::setReal(itsAttr[SPTDUMPFREQ], sptDumpFreq);
    Attributes::setReal(itsAttr[SCSOLVEFREQ], scSolveFreq);
    Attributes::setReal(itsAttr[MTSSUBSTEPS], mtsSubsteps);
    Attributes::setReal(itsAttr[REMOTEPARTDEL], remotePartDel);
    Attributes::setReal(itsAttr[REPARTFREQ], repartFreq);
    Attributes::setReal(itsAttr[REBINFREQ], rebinFreq);
    Attributes::setBool(itsAttr[SCAN], scan);
    Attributes::setBool(itsAttr[RHODUMP], rhoDump);
    Attributes::setBool(itsAttr[EBDUMP], ebDump);
    Attributes::setBool(itsAttr[CSRDUMP], csrDump);
    Attributes::setReal(itsAttr[AUTOPHASE], autoPhase);
    Attributes::setBool(itsAttr[PPDEBUG], ppdebug);
    Attributes::setReal(itsAttr[SURFDUMPFREQ], surfDumpFreq);
    Attributes::setBool(itsAttr[CZERO], cZero);
    Attributes::setBool(itsAttr[SCHOTTKYCORR], schottkyCorrection);
    Attributes::setBool(itsAttr[CLOTUNEONLY], cloTuneOnly);
    Attributes::setString(itsAttr[RNGTYPE], std::string(rngtype));
    Attributes::setReal(itsAttr[SCHOTTKYRENO], schottkyRennormalization);
    Attributes::setReal(itsAttr[NUMBLOCKS], numBlocks);
    Attributes::setReal(itsAttr[RECYCLEBLOCKS], recycleBlocks);
    Attributes::setReal(itsAttr[NLHS], nLHS);
    Attributes::setBool(itsAttr[ENABLEHDF5], enableHDF5);
    Attributes::setBool(itsAttr[ASCIIDUMP], asciidump);
    Attributes::setReal(itsAttr[BOUNDPDESTROYFQ], 10);
    Attributes::setReal(itsAttr[BEAMHALOBOUNDARY], 0);
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
    verify    = Attributes::getBool(itsAttr[VERIFY]);
    warn      = Attributes::getBool(itsAttr[WARN]);
    psDumpEachTurn =   Attributes::getBool(itsAttr[PSDUMPEACHTURN]);
    psDumpLocalFrame = Attributes::getBool(itsAttr[PSDUMPLOCALFRAME]);
    scan = Attributes::getBool(itsAttr[SCAN]);
    rhoDump = Attributes::getBool(itsAttr[RHODUMP]);
    ebDump = Attributes::getBool(itsAttr[EBDUMP]);
    csrDump = Attributes::getBool(itsAttr[CSRDUMP]);
    ppdebug = Attributes::getBool(itsAttr[PPDEBUG]);
    enableHDF5 = Attributes::getBool(itsAttr[ENABLEHDF5]);

    if(itsAttr[ASCIIDUMP]) {
        asciidump = Attributes::getBool(itsAttr[ASCIIDUMP]);
    }

    if(itsAttr[SEED]) {
        seed = int(Attributes::getReal(itsAttr[SEED]));
	if (seed == -1)
            rangen.init55(time(0));
	else
            rangen.init55(seed);
    }

    if(itsAttr[PSDUMPFREQ]) {
        psDumpFreq = int(Attributes::getReal(itsAttr[PSDUMPFREQ]));
    }

    if(itsAttr[STATDUMPFREQ]) {
        statDumpFreq = int(Attributes::getReal(itsAttr[STATDUMPFREQ]));
    }

    if(itsAttr[SPTDUMPFREQ]) {
        sptDumpFreq = int(Attributes::getReal(itsAttr[SPTDUMPFREQ]));
    }


    if(itsAttr[SCSOLVEFREQ]) {
        scSolveFreq = int(Attributes::getReal(itsAttr[SCSOLVEFREQ]));
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

    if(itsAttr[SCHOTTKYCORR]) {
        schottkyCorrection = bool(Attributes::getBool(itsAttr[SCHOTTKYCORR]));
    } else {
        schottkyCorrection = false;
    }

    if(itsAttr[SCHOTTKYRENO]) {
        schottkyRennormalization = double(Attributes::getReal(itsAttr[SCHOTTKYRENO]));
    } else {
        schottkyRennormalization = -1.0;
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
}
