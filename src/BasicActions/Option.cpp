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
#include <iostream>

using namespace Options;


// Class Option
// ------------------------------------------------------------------------

// The global option flags.
namespace Options {
  // The global program options.
  bool echo = true;
  bool info = true;
  bool mad8 = false;
  bool opal8 = false;	
  bool mtrace = false;
  bool verify = false;
  bool warn = true;
  bool tfsFormat = true;
  bool phaseScan = false;

  // The global random generator.
  Random rangen;

  // The current random seed.
  int seed  = 123456789;

  // The frequency to dump the phase space i.e. timestep%psDumpFreq
  int psDumpFreq = 10;

  // The frequency to dump single particle trajectory of particles with ID = 0 & 1 
  int sptDumpFreq = 1;

  // The frequency to do particles repartition for better load balance between nodes 
  int repartFreq = 10;

  /// Defines if an when we do the eake field calculation
  int wakeCalcStep = -1;  

  /// Defines the number of bins for the line density calculation
  int wakeBins = 64;

}


namespace {
  // The attributes of class Option.
  enum {
    ECHO,
    INFO,
    MAD8,
    OPAL8,
    TRACE,
    VERIFY,
    WARN,
    TFS_FORMAT,
    SEED,
    TELL,
    PHASESCAN,
    PSDUMPFREQ,
    SPTDUMPFREQ,
    REPARTFREQ,
    WAKEBINS,
    WAKECALC,
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
  itsAttr[MAD8] = Attributes::makeBool
    ("MAD8", "If true, write output in MAD-8 format", mad8);
  itsAttr[OPAL8] = Attributes::makeBool
    ("OPAL8", "If true, write output in OPAL-8 format", mad8);
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
  itsAttr[PHASESCAN] = Attributes::makeBool
    ("PHASESCAN", "If true, track only one or a small number of particles with no space charge", false);
  itsAttr[PSDUMPFREQ] = Attributes::makeReal
    ("PSDUMPFREQ", "The frequency to dump the phase space i.e. timestep%psDumpFreq");
  itsAttr[SPTDUMPFREQ] = Attributes::makeReal
    ("SPTDUMPFREQ", "The frequency to dump single particle trajectory of particles with ID = 0 & 1 ");
  itsAttr[REPARTFREQ] = Attributes::makeReal
    ("REPARTFREQ", "The frequency to do particles repartition for better load balance between nodes ");
  itsAttr[WAKECALC] = Attributes::makeReal
    ("WAKECALC", "Defines when to calculate the wake fields.Default is the value  -1 i.e. no wake field calculation."); 
  itsAttr[WAKEBINS] = Attributes::makeReal
    ("WAKEBINS", "Defines the number of bins used for obtaining the line density in the wakefield calculaion. Its default value is 64");
   
  FileStream::setEcho(echo);
  rangen.init55(seed);
}


Option::Option(const string &name, Option *parent):
  Action(name, parent)
{
  Attributes::setBool(itsAttr[ECHO],       echo);
  Attributes::setBool(itsAttr[INFO],       info);
  Attributes::setBool(itsAttr[MAD8],       mad8);
  Attributes::setBool(itsAttr[OPAL8],      opal8);
  Attributes::setBool(itsAttr[TRACE],      mtrace);
  Attributes::setBool(itsAttr[VERIFY],     verify);
  Attributes::setBool(itsAttr[WARN],       warn);
  Attributes::setBool(itsAttr[TFS_FORMAT], tfsFormat);
  Attributes::setReal(itsAttr[SEED],       seed);
  Attributes::setBool(itsAttr[PHASESCAN],  phaseScan);
  Attributes::setReal(itsAttr[PSDUMPFREQ], psDumpFreq);
  Attributes::setReal(itsAttr[SPTDUMPFREQ], sptDumpFreq);
  Attributes::setReal(itsAttr[REPARTFREQ], repartFreq);
  Attributes::setReal(itsAttr[WAKECALC], wakeCalcStep);  
  Attributes::setReal(itsAttr[WAKEBINS], wakeBins);

}


Option::~Option()
{}
  

Option *Option::clone(const string &name)
{
  return new Option(name, this);
}


void Option::execute()
{
  // Store the option flags.
  echo      = Attributes::getBool(itsAttr[ECHO]);
  info      = Attributes::getBool(itsAttr[INFO]);
  mad8      = Attributes::getBool(itsAttr[MAD8]);
  opal8     = Attributes::getBool(itsAttr[OPAL8]);
  mtrace     = Attributes::getBool(itsAttr[TRACE]);
  verify    = Attributes::getBool(itsAttr[VERIFY]);
  warn      = Attributes::getBool(itsAttr[WARN]);
  tfsFormat = Attributes::getBool(itsAttr[TFS_FORMAT]);
  phaseScan = Attributes::getBool(itsAttr[PHASESCAN]);


  if (itsAttr[SEED]) {
    seed = int(Attributes::getReal(itsAttr[SEED]));
    rangen.init55(seed);
  }

  if (itsAttr[PSDUMPFREQ]) {
    psDumpFreq = int(Attributes::getReal(itsAttr[PSDUMPFREQ]));
  }
  
  if (itsAttr[SPTDUMPFREQ]) {
    sptDumpFreq = int(Attributes::getReal(itsAttr[SPTDUMPFREQ]));
  }
  
  if (itsAttr[REPARTFREQ]) {
    repartFreq = int(Attributes::getReal(itsAttr[REPARTFREQ]));
  }

  if (itsAttr[WAKECALC])
    wakeCalcStep = int(Attributes::getReal(itsAttr[WAKECALC]));

  if (itsAttr[WAKEBINS])
    wakeBins = int(Attributes::getReal(itsAttr[WAKEBINS]));

  // Set message flags.
  FileStream::setEcho(echo);
 
  if (Attributes::getBool(itsAttr[TELL])) {
    std::cerr << "\nCurrent settings of options:\n" << *this << std::flush;
  }
}
