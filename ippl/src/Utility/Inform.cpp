// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Utility/Inform.h"

#include "Utility/IpplInfo.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <fstream>
using namespace std;
#else
#include <fstream.h>
#endif

#include <cstring>

/////////////////////////////////////////////////////////////////////
// manipulator functions

// signal we wish to send the message
Inform& endl(Inform& inf)   { inf << '\n'; return inf.outputMessage(); }

// set the current msg level
Inform& level1(Inform& inf) { return inf.setMessageLevel(1); }
Inform& level2(Inform& inf) { return inf.setMessageLevel(2); }
Inform& level3(Inform& inf) { return inf.setMessageLevel(3); }
Inform& level4(Inform& inf) { return inf.setMessageLevel(4); }
Inform& level5(Inform& inf) { return inf.setMessageLevel(5); }


/////////////////////////////////////////////////////////////////////
// perform initialization for this object; called by the constructors.
// arguments = prefix string, print node
void Inform::setup(const char *myname, int pnode) {
  
  On = true;
  OutputLevel = MsgLevel = MIN_INFORM_LEVEL;
  PrintNode = pnode;

  if ( myname != 0 )
    Name = strcpy(new char[strlen(myname) + 1], myname);
  else
    Name = 0;
}


/////////////////////////////////////////////////////////////////////
// class constructor
Inform::Inform(const char *myname, int pnode)
#ifdef IPPL_NO_STRINGSTREAM
    : FormatBuf(ios::out), OpenedSuccessfully(true) {
#else
    : FormatBuf(MsgBuf, MAX_INFORM_MSG_SIZE, ios::out), 
        OpenedSuccessfully(true) {
#endif
  

  // in this case, the default destination stream is cout
  NeedClose = false;
  MsgDest = &cout;

  // perform all other needed initialization
  setup(myname, pnode);
}


/////////////////////////////////////////////////////////////////////
// class constructor specifying a file to open
Inform::Inform(const char *myname, const char *fname, const WriteMode opnmode,
	       int pnode)
#ifdef IPPL_NO_STRINGSTREAM
    : FormatBuf(ios::out), OpenedSuccessfully(true) {
#else
    : FormatBuf(MsgBuf, MAX_INFORM_MSG_SIZE, ios::out), 
        OpenedSuccessfully(true) {
#endif

  // only open a file if we're on the proper node
  MsgDest = 0;
  if (pnode >= 0 && pnode == Ippl::myNode()) {
    if (opnmode == OVERWRITE)
      MsgDest = new ofstream(fname, ios::out);
    else
      MsgDest = new ofstream(fname, ios::app);
  }

  // make sure it was opened properly
  if ( MsgDest == 0 || ! (*MsgDest) ) {
    if (pnode >= 0 && pnode == Ippl::myNode()) {
      cerr << "Inform: Cannot open file '" << fname << "'." << endl;
    }
    NeedClose = false;
    MsgDest = &cout;
    OpenedSuccessfully = false;
  } else {
    NeedClose = true;
  }

  // perform all other needed initialization
  setup(myname, pnode);
}


/////////////////////////////////////////////////////////////////////
// class constructor specifying an output stream to use
Inform::Inform(const char *myname, std::ostream& os, int pnode)
#ifdef IPPL_NO_STRINGSTREAM
  : FormatBuf(ios::out), OpenedSuccessfully(true) {
#else
  : FormatBuf(MsgBuf, MAX_INFORM_MSG_SIZE, ios::out), 
      OpenedSuccessfully(true) {
#endif

  // just store a ref to the provided stream
  NeedClose = false;
  MsgDest = &os;

  // perform all other needed initialization
  setup(myname, pnode);
}


/////////////////////////////////////////////////////////////////////
// class destructor ... frees up space
Inform::~Inform(void) {
  
  delete [] Name;
  if ( NeedClose )
    delete MsgDest;
}


// print out just a single line, from the given buffer
void Inform::display_single_line(char *buf) {

  // output the prefix name if necessary ... if no name was given, do
  // not print any prefix at all
  if ( Name != 0 ) {
    *MsgDest << Name;

    // output the node number if necessary
    if (Ippl::getNodes() > 1)
      *MsgDest << "{" << Ippl::myNode() << "}";

    // output the message level number if necessary
    if ( MsgLevel > 1 )
      *MsgDest << "[" << MsgLevel << "]";

    // output the end of the prefix string if necessary
    if ( Name != 0)
      *MsgDest << "> ";
  }

  // finally, print out the message itself
  *MsgDest << buf << endl;
}


/////////////////////////////////////////////////////////////////////
// Print out the message in the given buffer.
void Inform::display_message(char *buf) {

  // check if we should even print out the message
  if ( On && MsgLevel <= OutputLevel && buf != 0 ) {
    // get location of final string term char
    char *stend = buf + strlen(buf);

    // print blank lines for leading endlines
    while (*buf == '\n') {
      *buf = '\0';
      display_single_line(buf++);
    }

    // print out all lines in the string now
    while ( (buf = strtok(buf, "\n")) != 0 ) {
      display_single_line(buf);
      buf += strlen(buf);
      if (buf < stend)
	buf++;

      // print out contiguous blank lines, if any
      while (*buf == '\n') {
	*buf = '\0';
	display_single_line(buf++);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////
// Set the current output level for this Inform object.
Inform& Inform::setOutputLevel(const int ol) {
  

  if ( ol >= (MIN_INFORM_LEVEL-1) && ol <= MAX_INFORM_LEVEL )
    OutputLevel = ol;
  return *this;
}


/////////////////////////////////////////////////////////////////////
// Set the current message level for the current message in this Inform object.
Inform& Inform::setMessageLevel(const int ol) {
  
  if ( ol >= MIN_INFORM_LEVEL && ol <= MAX_INFORM_LEVEL )
    MsgLevel = ol;
  return *this;
}


/////////////////////////////////////////////////////////////////////
// the signal has been given ... process the message.  Return ref to object.
Inform& Inform::outputMessage(void) {
  
  // print out the message (only if this is the master node)
  if (PrintNode < 0 || PrintNode == Ippl::myNode()) {
    FormatBuf << ends;
#ifdef IPPL_NO_STRINGSTREAM
    // extract C string and display
    MsgBuf = FormatBuf.str();
    char* cstring = const_cast<char*>(MsgBuf.c_str());
    display_message(cstring);
    // clear buffer contents
    // MsgBuf = string("");
    // FormatBuf.str(MsgBuf);
#else
    display_message(MsgBuf);
#endif
  }

  // reset this ostrstream to the start
  FormatBuf.seekp(0, ios::beg);
  return *this;
}


/////////////////////////////////////////////////////////////////////
// test program

#ifdef DEBUG_INFORM_CLASS

int main(int argc, char *argv[]) {
  
  int i;

  // create an Inform instance
  Inform inf("Inform Test");

  // copy in the argv's ... then print them out
  for ( i=0; i < argc ; i++)
    inf << "Argument " << i << " = " << argv[i] << "\n";
  inf << endl << endl;

   // do another one to make sure
  inf.setOutputLevel(3);
  inf << level2 << "This is the second test." << endl;

  return 0;
}

#endif

/***************************************************************************
 * $RCSfile: Inform.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: Inform.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
