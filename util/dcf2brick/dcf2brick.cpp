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
#include "dcf2brick.h"
#include <iostream.h>
#include <stdlib.h>

/***************************************************************************
  dcf2brick converts a given IPPL DiscField file to the format used as
  input to the ACL visualization team's multi-pipe renderer, and other
  programs, where all that is needed is a set of "bricks", where each
  brick is a simple fortran-order block of unsigned char's.  It can convert
  the contents of a 1, 2, or 3-dimensional field of doubles
  or Vektors, Tenzors, SymTenzors of doubles
  (it can also convert fields of floats or ints).

  The output will go to a set of files, named in the format
    outname.NNNN.dat
  where NNNN is a 4-digit number for the record in the file.  Command-line
  options can select whether to do the whole file or just some records.

  Exit codes for dcf2brick:
    0 : OK, files converted successfully.
    1 : Error in reading.
    2 : Error in writing.
    3 : Error in input parameters.
    4 : Unsupported dimension.
 ***************************************************************************/


////////////////////////////////////////////////////////////////////////////
// print out a simple usage message and exit
void usage(int, char *argv[], int retval)
{
  Inform msg("Usage");
  msg << argv[0] << " -file <infile> [-ofile <outfile>] [options]" << endl;
  msg << "  where [options] may be any of the following:";
  msg << "\n   -conf <file> : configuration file";
  msg << "\n   -field <n> : which field to print from file (0 to numfields-1)";
  msg << "\n   -startrecord <n> : starting rec from file (0 to numrecs-1)";
  msg << "\n   -endrecord <n> : ending rec from file (0 to numrecs-1)";
  msg << "\n   -minsize <n> : minimum size of the output field (pow of 2)";
  msg << "\n   -redominmax : recalculate min & max each record, instead of global min/max";
  msg << "\n   -sum : combine multipple-field file into one field";
  msg << "\n   -dim <dim> : dim of field (will be inferred if not given)";
  msg << "\n   -vec or -vector: assume vector instead of scalar data";
  msg << "\n   -tens or -tensor: assume tensor instead of scalar data";
  msg << "\n   -sym or -symtensor: assume symtensor instead of scalar data";
  msg << "\n   -anti or -antisymtensor: assume antisymtensor instead of scalar data";
  msg << "\n   -i or -int : assume integer instead of double data";
  msg << "\n   -f or -flt : assume float instead of double data";
  msg << "\n   -d or -dbl : assume double for the data (default)";
  msg << "\n\n  If -ofile option is not given, then the output filenames ";
  msg << "will be <infile>.NNNN.dat.";
  msg << "\n\n  Possible error codes:";
  msg << "\n   0: OK, file converted successfully";
  msg << "\n   1: Error in reading";
  msg << "\n   2: Error in writing";
  msg << "\n   3: Error in command-line parameters";
  msg << "\n   4: Unsupported dimension";
  msg << "\n\nUse '--connect aclvis' to also visually display the data";
  msg << endl;

  exit(retval);
}


////////////////////////////////////////////////////////////////////////////
// main routine
int main(int argc, char *argv[]) {
  int i, retval = RC_CMDERR;

  // scan early for -h flag
  for (i=1; i < argc; ++i)
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?"))
      usage(argc, argv, RC_OK);

  // initialize ippl
  Ippl ippl(argc,argv);
  Inform errmsg("dump error");

  // we need some kind of arguments
  if (argc == 1)
    usage(argc, argv, retval);

  // initial values for settings set by command-line arguments
  char *config = 0;
  char *fname = 0;
  char *outfname = 0;
  int datatype = CONVDOUBLE;
  int atype = CONVSCALAR;
  unsigned int dim = 0;
  int fieldnum = 0;
  int startrecordnum = (-1);
  int endrecordnum = (-1);
  bool dosum = false;
  bool globalminmax = true;
  int minsize = 0;

  // parse command-line options to get info
  for (i=1; i < argc; ++i) {
    if (!strcmp(argv[i], "-conf")) {
      if (i < (argc - 1))
	config = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-file")) {
      if (i < (argc - 1))
	fname = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-ofile")) {
      if (i < (argc - 1))
	outfname = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-dim")) {
      if (i < (argc - 1))
	dim = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-minsize")) {
      if (i < (argc - 1))
	minsize = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-field")) {
      if (i < (argc - 1))
	fieldnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-startrecord")) {
      if (i < (argc - 1))
	startrecordnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-endrecord")) {
      if (i < (argc - 1))
	endrecordnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-sum")) {
      dosum = true;
    } else if (!strcmp(argv[i], "-redominmax")) {
      globalminmax = false;

    } else if (!strcmp(argv[i], "-vec") || !strcmp(argv[i], "-vector")) {
      atype = CONVVECTOR;
    } else if (!strcmp(argv[i], "-tens") || !strcmp(argv[i], "-tensor")) {
      atype = CONVTENSOR;
    } else if (!strcmp(argv[i], "-sym") || !strcmp(argv[i], "-symtensor")) {
      atype = CONVSYMTENSOR;
    } else if (!strcmp(argv[i], "-anti") || !strcmp(argv[i],"-antisymtensor")){
      atype = CONVANTISYMTENSOR;

    } else if (!strcmp(argv[i], "-int") || !strcmp(argv[i], "-i")) {
      datatype = CONVINT;
    } else if (!strcmp(argv[i], "-flt") || !strcmp(argv[i], "-f")) {
      datatype = CONVFLOAT;
    } else if (!strcmp(argv[i], "-dbl") || !strcmp(argv[i], "-d")) {
      datatype = CONVDOUBLE;

    } else {
      errmsg << "Unknown option '" << argv[i] << "'." << endl;
      usage(argc, argv, retval);
    }
  }

  // sanity checking on args
  if (fname == 0) {
    errmsg << "You must provide a base file name." << endl;
    usage(argc, argv, retval);
  }

  // figure out the dimension of the fields, and make sure they are legal
  // note that the -dim argument from the user is just ignored.
  DiscField<1> io(fname, config);
  dim = io.get_Dimension();

  // figure out which form of the print routine we need to use, and call it
  retval = RC_UNSUPPORTED;
  if (datatype == CONVDOUBLE) {
    if (dim == 1)
      retval = conv_data_d_1(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 2)
      retval = conv_data_d_2(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 3)
      retval = conv_data_d_3(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 4)
      retval = conv_data_d_4(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else
      errmsg << "Unsupported dimension " << dim << "." << endl;

  } else if (datatype == CONVFLOAT) {
    if (dim == 1)
      retval = conv_data_f_1(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 2)
      retval = conv_data_f_2(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 3)
      retval = conv_data_f_3(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 4)
      retval = conv_data_f_4(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else
      errmsg << "Unsupported dimension " << dim << "." << endl;

  } else if (datatype == CONVINT) {
    if (dim == 1)
      retval = conv_data_i_1(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 2)
      retval = conv_data_i_2(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 3)
      retval = conv_data_i_3(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else if (dim == 4)
      retval = conv_data_i_4(atype, fname, config, outfname, fieldnum,
			     startrecordnum, endrecordnum,
			     globalminmax, minsize, dosum, errmsg);
    else
      errmsg << "Unsupported dimension " << dim << "." << endl;

  } else {
    errmsg << "Unsupported data type." << endl;
  }

  return retval;
}


/***************************************************************************
 * $RCSfile: dcf2brick.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:38 $
 * IPPL_VERSION_ID: $Id: dcf2brick.cpp,v 1.1.1.1 2003/01/23 07:41:38 adelmann Exp $ 
 ***************************************************************************/
