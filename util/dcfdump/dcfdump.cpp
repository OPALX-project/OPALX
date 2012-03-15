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

#include "Ippl.h"
#include "dcfdump.h"

/***************************************************************************
  dcfdump is a simple utility program which prints out
  the contents of a 1, 2, or 3-dimensional field of doubles
  or Vektors of doubles (can also print fields of floats or ints).

  Exit codes for dcfdump:
    0 : OK, files printed successfully.
    1 : Error in reading or printing.
    2 : Error in input parameters.
    3 : Unsupported dimension or type.
 ***************************************************************************/

// some global variables set by the main routine, used in the dump routine
bool doprint = true;
bool dumpoffset = false;
bool dumplayout = false;

////////////////////////////////////////////////////////////////////////////
// print out a simple usage message and exit
void usage(int, char *argv[], int retval) {
  Inform msg("Usage");
  msg << argv[0] << " -file <file> [options]" << endl;
  msg << "  where [options] may be any of the following:";
  msg << "\n   -h or -? : print this help info and quit.";
  msg << "\n   -conf <file> : configuration file";
  msg << "\n   -field <n> : which field to print from file (0 to numfields-1)";
  msg << "\n   -record <n> : which rec. to print from file (0 to numrecs-1)";
  msg << "\n   -dim <dim> : dim of field (will be inferred if not given)";
  msg << "\n   -vec : assume vector instead of scalar data";
  msg << "\n   -tens or -tensor: assume tensor instead of scalar data";
  msg << "\n   -sym or -symtensor: assume symtensor instead of scalar data";
  msg << "\n   -anti or -antisymtensor: assume antisymtensor instead of scalar data";
  msg << "\n   -i or -int : assume integer instead of double data";
  msg << "\n   -f or -flt : assume float instead of double data";
  msg << "\n   -d or -dbl : assume double for the data (default)";
  msg << "\n   -noprint : just print a summary, not the actual field data";
  msg << "\n";
  msg << "\n  Possible error codes:";
  msg << "\n   0: OK, file printed successfully";
  msg << "\n   1: Error in reading or printing";
  msg << "\n   2: Error in command-line parameters";
  msg << "\n   3: Unsupported dimension or type";
  msg << "\n";
  msg << "\nUse '--connect aclvis' to also visually display the data";
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
  Inform errmsg("dump error");

  // we need some kind of arguments
  if (argc == 1)
    usage(argc, argv, retval);

  // parse command-line options to get info
  char *c1 = 0;
  char *f1 = 0;
  unsigned int dim = 0;
  int datatype = DUMPDOUBLE;
  int fnum = (-1);
  int rnum = (-1);
  int atype = DUMPSCALAR;

  for (i=1; i < argc; ++i) {
    if (!strcmp(argv[i], "-conf")) {
      if (i < (argc - 1))
	c1 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-file")) {
      if (i < (argc - 1))
	f1 = argv[++i];
      else
	usage(argc, argv, retval);
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
    } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?")) {
      usage(argc, argv, 0);

    } else if (!strcmp(argv[i], "-vec") || !strcmp(argv[i], "-vector")) {
      atype = DUMPVECTOR;
    } else if (!strcmp(argv[i], "-tens") || !strcmp(argv[i], "-tensor")) {
      atype = DUMPTENSOR;
    } else if (!strcmp(argv[i], "-sym") || !strcmp(argv[i], "-symtensor")) {
      atype = DUMPSYMTENSOR;
    } else if (!strcmp(argv[i], "-anti")||!strcmp(argv[i], "-antisymtensor")) {
      atype = DUMPANTISYMTENSOR;

    } else if (!strcmp(argv[i], "-int") || !strcmp(argv[i], "-i")) {
      datatype = DUMPINT;
    } else if (!strcmp(argv[i], "-flt") || !strcmp(argv[i], "-f")) {
      datatype = DUMPFLOAT;
    } else if (!strcmp(argv[i], "-dbl") || !strcmp(argv[i], "-d")) {
      datatype = DUMPDOUBLE;

    } else if (!strcmp(argv[i], "-noprint")) {
      doprint = false;
    } else if (!strcmp(argv[i], "-offset")) {
      dumpoffset = true;
    } else if (!strcmp(argv[i], "-layout")) {
      dumplayout = true;

    } else {
      errmsg << "Unknown option '" << argv[i] << "'." << endl;
      usage(argc, argv, retval);
    }
  }

  // sanity checking on args
  if (f1 == 0) {
    errmsg << "You must provide a base file name." << endl;
    usage(argc, argv, retval);
  }

  // figure out the dimension of the fields, and make sure they are legal
  // note that the -dim argument from the user is just ignored.
  DiscField<1> io(f1, c1);
  dim = io.get_Dimension();

  // figure out which form of the print routine we need to use, and call it
  retval = RC_UNSUPPORTED;
  if (datatype == DUMPDOUBLE) {
    if (dim == 1)
      retval = dump_data_d_1(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 2)
      retval = dump_data_d_2(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 3)
      retval = dump_data_d_3(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 4)
      retval = dump_data_d_4(atype, f1, c1, fnum, rnum, errmsg);
    else
      errmsg << "Unsupported dimension " << dim << "." << endl;
  } else if (datatype == DUMPFLOAT) {
    if (dim == 1)
      retval = dump_data_f_1(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 2)
      retval = dump_data_f_2(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 3)
      retval = dump_data_f_3(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 4)
      retval = dump_data_f_4(atype, f1, c1, fnum, rnum, errmsg);
    else
      errmsg << "Unsupported dimension " << dim << "." << endl;
  } else if (datatype == DUMPINT) {
    if (dim == 1)
      retval = dump_data_i_1(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 2)
      retval = dump_data_i_2(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 3)
      retval = dump_data_i_3(atype, f1, c1, fnum, rnum, errmsg);
    else if (dim == 4)
      retval = dump_data_i_4(atype, f1, c1, fnum, rnum, errmsg);
    else
      errmsg << "Unsupported dimension " << dim << "." << endl;
  } else {
    errmsg << "Unsupported data type." << endl;
  }

  return retval;
}


/***************************************************************************
 * $RCSfile: dcfdump.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcfdump.cpp,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/
