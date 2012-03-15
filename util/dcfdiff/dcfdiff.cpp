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
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "dcfdiff_template.h"
#include "dcfdiff.h"

/***************************************************************************
  dcfdiff compares the contents of two fields (which must have the same
  dimensions and total domains), and computes RMSD's, max diffs, min diffs,
  and total diffs.

  Exit codes for dcfdiff:
    0 : OK, files compared within specified tolerance or no tolerance given
    1 : Error in reading data files
    2 : Error in command-line parameters or field configurations
    3 : Files read properly, but failed to compare within specified rmsd
    4 : Unsupported dimension
 ***************************************************************************/

////////////////////////////////////////////////////////////////////////////
// global variables used by diff routine
bool printMaxDiff = false;


////////////////////////////////////////////////////////////////////////////
// print out a simple usage message and exit
void usage(int, char *argv[], int retval) {
  Inform msg("Usage");
  msg << argv[0] << " -file1 <file> -file2 <file> [options]" << endl;
  msg << "  where [options] may be any of the following:";
  msg << "\n   -h or -? : print this help info and quit.";
  msg << "\n   -conf1 <file> : configuration file for file1";
  msg << "\n   -conf2 <file> : configuration file for file2";
  msg << "\n   -field <n> : which field to use from file (0 to numfields-1)";
  msg << "\n   -record <n> : which rec. to use from file (0 to numrecs-1)";
  msg << "\n   -maxrmsd <val> : maximum rmsd value. If the rmsd for a record";
  msg << "\n                    is >= maxrmsd, return error code value of 3";
  msg << "\n   -dim <dim> : dim of field (will be inferred if not given)";
  msg << "\n   -vec or -vector: assume vector instead of scalar data";
  msg << "\n   -tens or -tensor: assume tensor instead of scalar data";
  msg << "\n   -sym or -symtensor: assume symtensor instead of scalar data";
  msg << "\n   -anti or -antisymtensor: assume antisymtensor instead of scalar";
  msg << "\n   -maxloc: print out the location of max diff";
  msg << "\n   -i or -int : assume integer instead of double data";
  msg << "\n   -f or -flt : assume float instead of double data";
  msg << "\n   -d or -dbl : assume double for the data (default)";
  msg << "\n";
  msg << "\n  Possible error codes:";
  msg << "\n   0: OK, files compared within tolerance or no tolerance given";
  msg << "\n   1: Error in reading data files";
  msg << "\n   2: Error in command-line parameters or field configurations";
  msg << "\n   3: Files read properly, but do not compare within given rmsd";
  msg << "\n   4: Unsupported dimension";
  msg << endl;
  exit(retval);
}


////////////////////////////////////////////////////////////////////////////
// main routine
int main(int argc, char *argv[]) {
  unsigned int i;
  int retval = RC_CMDERR;

  // scan early for -h flag
  for (i=1; i < argc; ++i)
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?"))
      usage(argc, argv, 0);

  // initialize ippl
  Ippl ippl(argc,argv);
  Inform errmsg("diff error");

  // we need some kind of arguments
  if (argc == 1)
    usage(argc, argv, retval);

  // parse command-line options to get info
  char *c1 = 0;
  char *c2 = 0;
  char *f1 = 0;
  char *f2 = 0;
  unsigned int dim = 0;
  int datatype = DIFFDOUBLE;
  int fnum = (-1);
  int rnum = (-1);
  double rmsd = -1.0;
  int atype = DIFFSCALAR;

  for (i=1; i < argc; ++i) {
    if (!strcmp(argv[i], "-conf1")) {
      if (i < (argc - 1))
	c1 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-file1")) {
      if (i < (argc - 1))
	f1 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-conf2")) {
      if (i < (argc - 1))
	c2 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-file2")) {
      if (i < (argc - 1))
	f2 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?")) {
      usage(argc, argv, 0);
    } else if (!strcmp(argv[i], "-maxloc")) {
      printMaxDiff = true;
    } else if (!strcmp(argv[i], "-dim")) {
      if (i < (argc - 1))
	dim = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-field")) {
      if (i < (argc - 1))
	fnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-record")) {
      if (i < (argc - 1))
	rnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-maxrmsd")) {
      if (i < (argc - 1))
	rmsd = atof(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-vec") || !strcmp(argv[i], "-vector")) {
      atype = DIFFVECTOR;
    } else if (!strcmp(argv[i], "-tens") || !strcmp(argv[i], "-tensor")) {
      atype = DIFFTENSOR;
    } else if (!strcmp(argv[i], "-sym") || !strcmp(argv[i], "-symtensor")) {
      atype = DIFFSYMTENSOR;
    } else if (!strcmp(argv[i], "-anti")|| !strcmp(argv[i], "-antisymtensor")){
      atype = DIFFANTISYMTENSOR;
    } else if (!strcmp(argv[i], "-int") || !strcmp(argv[i], "-i")) {
      datatype = DIFFINT;
    } else if (!strcmp(argv[i], "-flt") || !strcmp(argv[i], "-f")) {
      datatype = DIFFFLOAT;
    } else if (!strcmp(argv[i], "-dbl") || !strcmp(argv[i], "-d")) {
      datatype = DIFFDOUBLE;
    } else {
      errmsg << "Unknown option '" << argv[i] << "'." << endl;
      usage(argc, argv, retval);
    }
  }

  // sanity checking on args
  if (f1 == 0 || f2 == 0) {
    errmsg << "You must provide two base file names." << endl;
    usage(argc, argv, retval);
  }

  // figure out the dimension of the fields, and make sure they are legal
  retval = RC_OK;
  DiscField<1> io1(f1, c1);
  DiscField<2> io2(f2, c2);
  dim = io1.get_Dimension();
  if (dim != io2.get_Dimension()) {
    errmsg << "Both files must have the same dimension." << endl;
    retval = RC_READERR;
  }

  // figure out which form of the diff routine we need to use, and call it
  if (retval == RC_OK) {
    retval = RC_UNSUPPORTED;
    if (datatype == DIFFDOUBLE) {
      if (dim == 1)
        retval = diff_data_d_1(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 2)
        retval = diff_data_d_2(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 3)
        retval = diff_data_d_3(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 4)
        retval = diff_data_d_4(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else
        errmsg << "Unsupported dimension " << dim << "." << endl;
    } else if (datatype == DIFFFLOAT) {
      if (dim == 1)
        retval = diff_data_f_1(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 2)
        retval = diff_data_f_2(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 3)
        retval = diff_data_f_3(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 4)
        retval = diff_data_f_4(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else
        errmsg << "Unsupported dimension " << dim << "." << endl;
    } else if (datatype == DIFFINT) {
      if (dim == 1)
        retval = diff_data_i_1(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 2)
        retval = diff_data_i_2(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 3)
        retval = diff_data_i_3(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else if (dim == 4)
        retval = diff_data_i_4(atype, f1, c1, f2, c2, fnum, rnum, rmsd, errmsg);
      else
        errmsg << "Unsupported dimension " << dim << "." << endl;
    } else {
      errmsg << "Unsupported data type." << endl;
    }
  }

  return retval;
}


/***************************************************************************
 * $RCSfile: dcfdiff.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcfdiff.cpp,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/
