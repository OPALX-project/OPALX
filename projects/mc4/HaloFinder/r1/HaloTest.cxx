/*=========================================================================

  Program:   HaloFinder driver program
  Module:    $RCSfile: HaloTest.cxx,v $

  Copyright (c) Chung-Hsing Hsu
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This example shows the driver for the complete HaloFinder regresstion
// test. It will take a .cosmo file containing all particles and output
// two files, one with all particles and the other with the halo catalog.
//

/*
#define  INFILE "../../HaloFinder.files/test.cosmo"
#define OUTFILE "test_serial.out"
#define NP 32
#define RL 100.0
#define BB 0.17
#define PMIN 10
#define PERIODIC true
#define TEXTMODE "ascii"
*/

#define  INFILE "../../HaloFinder.files/test4.cosmo"
#define OUTFILE "test4_serial.out"
#define NP 256
#define RL 90.140846
#define BB 0.20
#define PMIN 10
#define PERIODIC true
#define TEXTMODE "ascii"

#include <iostream>
#include "CosmoHaloFinder.h"

using namespace std;

int main(int argc, char* argv[])

{
  CosmoHaloFinder* haloFinder = new CosmoHaloFinder();
  haloFinder->np = NP;
  haloFinder->pmin = PMIN;
  haloFinder->bb = BB;
  haloFinder->rL = RL;
  haloFinder->periodic = PERIODIC;
  haloFinder->infile   = INFILE;
  haloFinder->outfile  = OUTFILE;
  haloFinder->textmode = TEXTMODE;

  haloFinder->Execute();

  return 0;
}
