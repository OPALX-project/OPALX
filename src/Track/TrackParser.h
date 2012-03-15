#ifndef OPAL_TrackParser_HH
#define OPAL_TrackParser_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackParser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackParser
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/OpalParser.h"
#include "AbstractObjects/Directory.h"


// Class TrackParser
// ------------------------------------------------------------------------
/// The parser class used by the OPAL tracking module.
//  As long as control remains in this class, OPAL recognises only the 
//  commands allowed in tracking mode.  Thus this parser has its own
//  command directory with a find() method which is used to find commands.

class TrackParser: public OpalParser {

public:

  TrackParser();
  virtual ~TrackParser();

protected:

  /// Find object by name in the track command directory.
  virtual Object *find(const string &) const;

private:

  // Not implemented.
  TrackParser(const TrackParser &);
  void operator=(const TrackParser &);

  // The sub-command directory.
  Directory trackDirectory;
};

#endif // OPAL_TrackParser_HH
