// =========================================================================
//
// Main.cc
//    CLASSIC test program
//
// Author: F. C. Iselin, CERN, SL/AP Group
//
// Copyright (c) 1995 (see file Copyright.hh)
//
// =========================================================================

#include "Parser/ClassicParser.h"
#include "Parser/FileStream.h"
#include "Parser/TerminalStream.h"
#include <iostream>


int main(int argc, char *argv[])
{
  std::cerr << "This is the CLASSIC test program, Version 1.0/1 (c) 1996."
       << std::endl;
  ClassicParser parser;

  if (argc < 2) {
    // Run commands from standard input
    TerminalStream is("CLASSIC");
    parser.run(&is);
  } else {
    FileStream is(argv[1]);
    std::cerr << "Reading input stream \"" << argv[1] << "\"." << std::endl;
    FileStream::setEcho(true);
    parser.run(&is);
    std::cerr << "End of input stream \"" << argv[1] << "\"." << std::endl;
  }

  return 0;
}
