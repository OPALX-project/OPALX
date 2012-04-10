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
#include "Utilities/Random.h"
#include "Ippl.h"

#include <iostream>

using namespace Options;

extern Inform *gmsg;

// Class Option
// ------------------------------------------------------------------------

// The global option flags.
namespace Options {
    // The global program options.
    bool echo = true;
    bool info = true;
    bool mtrace = false;
    bool verify = false;
    bool warn = true;
    bool tfsFormat = false;
    bool psDumpEachTurn = false;
    bool psDumpLocalFrame = false;
    bool scan = false;
    bool rhoDump = false;
    bool efDump = false;
    bool ppdebug = false;

    // The global random generator.
    Random rangen;

    // The current random seed.
    int seed  = 123456789;

    // the number of refinements of the search range for the phase with maximum energy
    // if eq 0 then no autophase
    int autoPhase = 0;

    // The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
    int psDumpFreq = 10;
    // The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
    double rDump = 0.0;

    // The frequency to dump statistical quantities such as beam RMS properties, i.e. dump
    // when step%statDumpFreq == 0.
    int statDumpFreq = 10;

    // The frequency to dump single particle trajectory of particles with ID = 0 & 1
    int sptDumpFreq = 1;

    // The frequency to do particles repartition for better load balance between nodes
    int repartFreq = 10;

    // The frequency to reset energy bin ID for all particles
    int rebinFreq = 100;

    /// The frequency to solve space charge fields.
    int scSolveFreq = 1;

    // How many small timesteps are inside the large timestep used in multiple time stepping (MTS) integrator
    int mtsSubsteps = 1;

    // The frequency to dump the particle-geometry surface interation data, -1 stands for no dump.
    int surfDumpFreq = -1;

    // Options for the Belos solver
    int numBlocks = 0;
    int recycleBlocks = 0;
    int nLHS = 1;

    // If true create symmetric distribution
    bool cZero = false;

    string rngtype = string("RANDOM");

    bool schottkyCorrection = false;

    double schottkyRennormalization = -1;

    // If true change the time step during emission from the cathode so that that one bin
    // of the time histogram that describes the longitudinal beam distribution is emitted
    // during each time step. If false then the time step during emission is set so that one
    // energy bin of the beam is emitted during each time step.
    bool fineEmission = true;
}


namespace {
    // The attributes of class Option.
    enum {
        ECHO,
        INFO,
        TRACE,
        VERIFY,
        WARN,
        TFS_FORMAT,
        SEED,
        TELL,
        PSDUMPFREQ,
        RDUMP,
        STATDUMPFREQ,
        PSDUMPEACHTURN,
        PSDUMPLOCALFRAME,
        SPTDUMPFREQ,
        REPARTFREQ,
        REBINFREQ,
        SCSOLVEFREQ,
        MTSSUBSTEPS,
        SCAN,
        RHODUMP,
        EFDUMP,
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
        FINEEMISSION,
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
    itsAttr[TFS_FORMAT] = Attributes::makeBool
                          ("TFS", "If true, print tables in TFS format", tfsFormat);
    itsAttr[SEED] = Attributes::makeReal
                    ("SEED", "The seed for the random generator");
    itsAttr[TELL] = Attributes::makeBool
                    ("TELL", "If true, print the current settings", false);
    itsAttr[PSDUMPFREQ] = Attributes::makeReal
                          ("PSDUMPFREQ", "The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0, its default value is 10.");
    itsAttr[RDUMP] = Attributes::makeReal
                     ("RDUMP", "Dump central beam when the radius is bigger than RDUMP");
    itsAttr[STATDUMPFREQ] = Attributes::makeReal
                            ("STATDUMPFREQ", "The frequency to dump statistical data (e.g. RMS beam quantities), i.e. dump data when step%statDumpFreq == 0, its default value is 10.");
    itsAttr[PSDUMPEACHTURN] = Attributes::makeBool
                              ("PSDUMPEACHTURN", "If true, dump phase space after each turn ,only aviable for OPAL-cycl, its default value is false");
    itsAttr[SCSOLVEFREQ] = Attributes::makeReal
                           ("SCSOLVEFREQ", "The frequency to solve space charge fields. its default value is 1");
    itsAttr[MTSSUBSTEPS] = Attributes::makeReal("MTSSUBSTEPS", "How many small timesteps are inside the large timestep used in multiple time stepping (MTS) integrator");
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

    itsAttr[EFDUMP] = Attributes::makeBool
                      ("EFDUMP", "If true, in addition to the phase space the E vector field is also dumped (H5Block)", efDump);

    itsAttr[AUTOPHASE] = Attributes::makeReal
                         ("AUTOPHASE", "If greater than zero OPAL is scaning the phases of each rf structure in order to get maximum acceleration. Defines the number of refinements of the search range", autoPhase);

    itsAttr[PPDEBUG] = Attributes::makeBool
                       ("PPDEBUG", "If true, use special initial velocity distribution for parallel plate and print special debug output", ppdebug);
    itsAttr[SURFDUMPFREQ] =  Attributes::makeReal
                             ("SURFDUMPFREQ", "The frequency to dump surface-partcle interaction data, its default value is -1 (no dump). ");

    itsAttr[CZERO] =  Attributes::makeBool
                      ("CZERO", "If set to true a symmetric distribution is created -> centroid == 0.0 ", cZero);

    itsAttr[RNGTYPE] =  Attributes::makeString
                        ("RNGTYPE", " RANDOM (default), Quasi-random number gernerators: HALTON, SOBOL, NIEDERREITER (Gsl ref manual 18.5)", rngtype);

    itsAttr[SCHOTTKYCORR] =  Attributes::makeBool
                                   ("SCHOTTKYCORR", "If set to true a Schottky correction to the charge is applied ", schottkyCorrection);

    itsAttr[SCHOTTKYRENO] =  Attributes::makeReal
                                         ("SCHOTTKYRENO", "IF set to a value greater than 0.0 the Schottky correction scan is disabled and the value is used for charge renormalization ", schottkyRennormalization);

    itsAttr[FINEEMISSION] = Attributes::makeBool
                            ("FINEEMISSION", "If true uses fine time step during particle emission from cathode.", fineEmission);

    itsAttr[NUMBLOCKS] = Attributes::makeReal
                          ("NUMBLOCKS", "Maximum number of vectors in the Krylov space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");
    itsAttr[RECYCLEBLOCKS] = Attributes::makeReal
                          ("RECYCLEBLOCKS", "Number of vectors in the recycle space (for RCGSolMgr). Default value is 0 and BlockCGSolMgr will be used.");
    itsAttr[NLHS] = Attributes::makeReal
                          ("NLHS", "Number of stored old solutions for extrapolating the new starting vector. Default value is 1 and just the last solution is used.");
    FileStream::setEcho(echo);
    rangen.init55(seed);
}


Option::Option(const string &name, Option *parent):
    Action(name, parent) {
    Attributes::setBool(itsAttr[ECHO],       echo);
    Attributes::setBool(itsAttr[INFO],       info);
    Attributes::setBool(itsAttr[TRACE],      mtrace);
    Attributes::setBool(itsAttr[VERIFY],     verify);
    Attributes::setBool(itsAttr[WARN],       warn);
    Attributes::setBool(itsAttr[TFS_FORMAT], tfsFormat);
    Attributes::setReal(itsAttr[SEED],       seed);
    Attributes::setReal(itsAttr[PSDUMPFREQ], psDumpFreq);
    Attributes::setReal(itsAttr[RDUMP], rDump);
    Attributes::setReal(itsAttr[STATDUMPFREQ], statDumpFreq);
    Attributes::setBool(itsAttr[PSDUMPEACHTURN], psDumpEachTurn);
    Attributes::setBool(itsAttr[PSDUMPLOCALFRAME], psDumpLocalFrame);
    Attributes::setReal(itsAttr[SPTDUMPFREQ], sptDumpFreq);
    Attributes::setReal(itsAttr[SCSOLVEFREQ], scSolveFreq);
    Attributes::setReal(itsAttr[MTSSUBSTEPS], mtsSubsteps);
    Attributes::setReal(itsAttr[REPARTFREQ], repartFreq);
    Attributes::setReal(itsAttr[REBINFREQ], rebinFreq);
    Attributes::setBool(itsAttr[SCAN], scan);
    Attributes::setBool(itsAttr[RHODUMP], rhoDump);
    Attributes::setBool(itsAttr[EFDUMP], efDump);
    Attributes::setReal(itsAttr[AUTOPHASE], autoPhase);
    Attributes::setBool(itsAttr[PPDEBUG], ppdebug);
    Attributes::setReal(itsAttr[SURFDUMPFREQ], surfDumpFreq);
    Attributes::setBool(itsAttr[CZERO], cZero);
    Attributes::setBool(itsAttr[SCHOTTKYCORR], schottkyCorrection);
    Attributes::setString(itsAttr[RNGTYPE], string(rngtype));
    Attributes::setReal(itsAttr[SCHOTTKYRENO], schottkyRennormalization);
    Attributes::setBool(itsAttr[FINEEMISSION], fineEmission);
    Attributes::setReal(itsAttr[NUMBLOCKS], numBlocks);
    Attributes::setReal(itsAttr[RECYCLEBLOCKS], recycleBlocks);
    Attributes::setReal(itsAttr[NLHS], nLHS);
}


Option::~Option()
{}


Option *Option::clone(const string &name) {
    return new Option(name, this);
}


void Option::execute() {
    // Store the option flags.
    echo      = Attributes::getBool(itsAttr[ECHO]);
    info      = Attributes::getBool(itsAttr[INFO]);
    mtrace     = Attributes::getBool(itsAttr[TRACE]);
    verify    = Attributes::getBool(itsAttr[VERIFY]);
    warn      = Attributes::getBool(itsAttr[WARN]);
    tfsFormat = Attributes::getBool(itsAttr[TFS_FORMAT]);
    psDumpEachTurn =   Attributes::getBool(itsAttr[PSDUMPEACHTURN]);
    psDumpLocalFrame = Attributes::getBool(itsAttr[PSDUMPLOCALFRAME]);
    scan = Attributes::getBool(itsAttr[SCAN]);
    rhoDump = Attributes::getBool(itsAttr[RHODUMP]);
    efDump = Attributes::getBool(itsAttr[EFDUMP]);
    ppdebug = Attributes::getBool(itsAttr[PPDEBUG]);



    if(itsAttr[SEED]) {
        seed = int(Attributes::getReal(itsAttr[SEED]));
        rangen.init55(seed);
    }

    if(itsAttr[PSDUMPFREQ]) {
        psDumpFreq = int(Attributes::getReal(itsAttr[PSDUMPFREQ]));
    }
    if(itsAttr[RDUMP]) {
        rDump = double(Attributes::getReal(itsAttr[RDUMP]));
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
        rngtype = string(Attributes::getString(itsAttr[RNGTYPE]));
    } else {
        rngtype = string("RANDOM");
    }

    if(itsAttr[FINEEMISSION]) {
        fineEmission = bool(Attributes::getBool(itsAttr[FINEEMISSION]));
    }

    // Set message flags.
    FileStream::setEcho(echo);

    if(Attributes::getBool(itsAttr[TELL])) {
        *gmsg << "\nCurrent settings of options:\n" << *this << endl;
    }
}
