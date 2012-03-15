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

#include <vector.h>
#include <stdlib.h>
#include "Ippl.h"

/***************************************************************************
  dcpdump prints the contents of a Particles files, either for a whole set
  of attributes or a single attribute.

  Exit codes for dcpdump:
    0 : OK, all attribs compared within specified tolerance or no tol given
    1 : Error in reading data files
    2 : Error in command-line parameters
    3 : Unsupported dimension or type
 ***************************************************************************/

// enumeration of possible error codes
enum { RC_OK, RC_READERR, RC_CMDERR, RC_UNSUPPORTED };


////////////////////////////////////////////////////////////////////////////
// print out a simple usage message and exit
void usage(int, char *argv[], int retval) {
  Inform msg("Usage");
  msg << argv[0] << " -file <file> [-types <typestring>] [opts]";
  msg << "\n  where [opts] may be any of the following:";
  msg << "\n   -h or -? : print this help info and quit.";
  msg << "\n   -conf <file> : configuration file";
  msg << "\n   -record <n> : which rec. to use from file (0 to numrecs-1)";
  msg << "\n   -attrib <n> : which attrib to print (does all if not given)";
  msg << "\n";
  msg << "\n  <typestring> is a string specifying the type of data for an";
  msg << "\n  attribute.  It should be a single string with a 1-3 char code";
  msg << "\n  for each attribute, with each code separated by dashes.  If";
  msg << "\n  a <typestring> is not included in the DiscParticle file for the";
  msg << "\n  attributes and record requested, you must provide it on the";
  msg << "\n  command-line here.  The type codes are of the format:";
  msg << "\n          [<structure><dim>]<type>";
  msg << "\n  where";
  msg << "\n";
  msg << "\n    <structure> == structure of data; it is a one-char code from:";
  msg << "\n        v ==> vector";
  msg << "\n        t ==> full tensor";
  msg << "\n        s ==> symmetric tensor";
  msg << "\n        a ==> antisymmetric tensor";
  msg << "\n    If <structure> is given, the next character must be a number";
  msg << "\n    from 1 ... 3, giving the dimension of the structure.";
  msg << "\n";
  msg << "\n    <type> == type of data; it is a one-char code from:";
  msg << "\n        c ==> character";
  msg << "\n        s ==> short";
  msg << "\n        i ==> integer";
  msg << "\n        l ==> long";
  msg << "\n        f ==> float";
  msg << "\n        d ==> double";
  msg << "\n        u ==> user or unknown type";
  msg << "\n";
  msg << "\n  Example:";
  msg << "\n    v3d-i-t2f : 3D Vektor (double), integer, 2D Tenzor (float)";
  msg << "\n";
  msg << "\n  If the <typestring> is only one character long, it is assumed";
  msg << "\n  to refer to a scalar of type <type>.  If it is three characters";
  msg << "\n  long, it is assumed to be a vector or tensor.  This program";
  msg << "\n  can only work with known types, if the type is 'u' it will";
  msg << "\n  quit with an error.";
  msg << "\n";
  msg << "\n  NOTE: The first two attributes in the typestring list MUST";
  msg << "\n  specify the type of the particle position and ID attributes,";
  msg << "\n  typically something like 'v3d i' for 3D Vektors and int's,";
  msg << "\n  respectively, if the data was saved as a whole Particles class";
  msg << "\n  instead of a single attribute.";
  msg << "\n";
  msg << "\n  Possible error codes:";
  msg << "\n   0: OK, attribs compared within tolerance or no tolerance given";
  msg << "\n   1: Error in reading data files";
  msg << "\n   2: Error in command-line parameters";
  msg << "\n   3: Unsupported dimension or type";
  msg << endl;
  exit(retval);
}


////////////////////////////////////////////////////////////////////////////
// a simple struct used to print out the values for a given attribute
template<class T>
struct PerformDump {
  static int dump(ParticleAttribBase &ba, int r, int attrib,
                  Inform &msg, Inform &errmsg) {

    // cast the attributes to type T
    ParticleAttrib<T> &a = static_cast<ParticleAttrib<T> &>(ba);

    // get the total number of particles
    int totalnum = a.size();
    reduce(totalnum, totalnum, OpAddAssign());
    
    // print banner
    if (Ippl::myNode() == 0) {
      msg << "\nContents of record " << r << ", attribute " << attrib;
      msg << " with " << totalnum << " total particles:" << endl;
      msg << "-------------------------------------------------" << endl;
    }

    // go through the nodes, and have them print their contents one at a time.
    // use barrier's to control who's doing what.
    int printed = 0;
    for (int n = 0; n < Ippl::getNodes(); ++n) {
      // current node prints its values
      if (n == Ippl::myNode()) {
        for (int p = 0; p < a.size(); ++p, ++printed) {
          msg << "attrib[" << attrib << "][" << printed << "] = ";
          msg << a[p] << endl;
        }
      } else {
        printed = 0;
      }

      // do a reduction to get number of printed values to all nodes
      // printed has been zeroed out on all but the printing node.  The next
      // node will get the new printed valued, and add to it before passing
      // it on.
      reduce(printed, printed, OpAddAssign());
    }

    return RC_OK;
  }

  // print values after reading them for a single attribute
  static int attrdump(DiscParticle &dp, int r, Inform &msg, Inform &errmsg) {

    // create attribute
    ParticleAttrib<T> a;

    // read in from the file
    int retval = RC_OK;
    if (!dp.read(a, r)) {
      errmsg << "Could not read attrib from file for record " << r << endl;
      retval = RC_READERR;
    } else {
      // print out the values
      retval = PerformDump<T>::dump(a, r, 0, msg, errmsg);
    }

    return retval;
  }
};


////////////////////////////////////////////////////////////////////////////
// a simple struct used to create a non-scalar attribute.  Given the dim,
// and the scalar type, this creates the attribute
template<class T, unsigned Dim>
struct CreateAttribTypeDim {
  // create an attribute
  static ParticleAttribBase *create(string &str, Inform &errmsg) {
    int atype = DiscTypeBase::appType(str);
    if (atype == DiscTypeBase::VEKTOR) {
      return new ParticleAttrib< Vektor<T, Dim> >;
    } else if (atype == DiscTypeBase::TENZOR) {
      return new ParticleAttrib< Tenzor<T, Dim> >;
    } else if (atype == DiscTypeBase::SYMTENZOR) {
      return new ParticleAttrib< SymTenzor<T, Dim> >;
    } else {
      errmsg << "In CreateAttribTypeDim::create: with str='" << str.c_str();
      errmsg << "', atype=" << atype << ": ";
      errmsg << "Unsupported attribute apptype." << endl;
      return 0;
    }
  }

  // dump an attribute based on its non-scalar type
  static int dump(ParticleAttribBase &attr, string &str, int r, int a,
		  Inform &msg, Inform &errmsg) {
    int atype = DiscTypeBase::appType(str);
    if (atype == DiscTypeBase::VEKTOR) {
      return PerformDump<Vektor<T, Dim> >::dump(attr, r, a, msg, errmsg);
    } else if (atype == DiscTypeBase::TENZOR) {
      return PerformDump<Tenzor<T, Dim> >::dump(attr, r, a, msg, errmsg);
    } else if (atype == DiscTypeBase::SYMTENZOR) {
      return PerformDump<SymTenzor<T, Dim> >::dump(attr, r, a, msg, errmsg);
    } else {
      errmsg << "In CreateAttribTypeDim::dump: with str='" << str.c_str();
      errmsg << "', atype=" << atype;
      errmsg << ", r=" << r << ", a=" << a << ": ";
      errmsg << "Unsupported attribute apptype." << endl;
      return RC_UNSUPPORTED;
    }
  }

  // dump just a single attribute based on its non-scalar type
  static int attrdump(DiscParticle &dp, vector<string> &str,
                      int r, Inform &msg, Inform &errmsg) {
    // based on the scalar type, invoke the proper dumping routine for
    // a record written using just one attrib
    int atype = DiscTypeBase::appType(str[0]);
    if (atype == DiscTypeBase::VEKTOR) {
      return PerformDump<Vektor<T, Dim> >::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::TENZOR) {
      return PerformDump<Tenzor<T, Dim> >::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::SYMTENZOR) {
      return PerformDump<SymTenzor<T, Dim> >::attrdump(dp, r, msg, errmsg);
    } else {
      errmsg << "In CreateAttribTypeDim::attrdump: with str='";
      errmsg << str[0].c_str();
      errmsg << "', atype=" << atype << ", r=" << r << ": ";
      errmsg << "Unsupported attribute apptype." << endl;
      return RC_UNSUPPORTED;
    }
  }
};


////////////////////////////////////////////////////////////////////////////
// a simple struct used to create a non-scalar attribute.  Given the dim,
// this determines the scalar type.
template<unsigned Dim>
struct CreateAttribDim {
  // create an attribute
  static ParticleAttribBase *create(string &str, Inform &errmsg) {
    int atype = DiscTypeBase::scalarType(str);
    if (atype == DiscTypeBase::CHAR) {
      return CreateAttribTypeDim<char, Dim>::create(str, errmsg);
    } else if (atype == DiscTypeBase::SHORT) {
      return CreateAttribTypeDim<short, Dim>::create(str, errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return CreateAttribTypeDim<int, Dim>::create(str, errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return CreateAttribTypeDim<long, Dim>::create(str, errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return CreateAttribTypeDim<float, Dim>::create(str, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return CreateAttribTypeDim<double, Dim>::create(str, errmsg);
    } else {
      errmsg << "In CreateAttribDim::create: with str='" << str.c_str();
      errmsg << "', atype=" << atype << ": ";
      errmsg << "Unsupported attribute apptype scalar type." << endl;
      return 0;
    }
  }

  // dump a non-scalar attribute 
  static int dump(ParticleAttribBase &attr, string &str, int r, int a,
		  Inform &msg, Inform &errmsg) {
    int atype = DiscTypeBase::scalarType(str);
    if (atype == DiscTypeBase::CHAR) {
      return CreateAttribTypeDim<char, Dim>::dump(attr,str,r,a,msg,errmsg);
    } else if (atype == DiscTypeBase::SHORT) {
      return CreateAttribTypeDim<short, Dim>::dump(attr,str,r,a,msg,errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return CreateAttribTypeDim<int, Dim>::dump(attr,str,r,a,msg,errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return CreateAttribTypeDim<long, Dim>::dump(attr,str,r,a,msg,errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return CreateAttribTypeDim<float, Dim>::dump(attr,str,r,a,msg,errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return CreateAttribTypeDim<double, Dim>::dump(attr,str,r,a,msg,errmsg);
    } else {
      errmsg << "In CreateAttribDim::dump: with str='" << str.c_str();
      errmsg << "', atype=" << atype;
      errmsg << ", r=" << r << ", a=" << a << ": ";
      errmsg << "Unsupported attribute apptype scalar type." << endl;
      return RC_UNSUPPORTED;
    }
  }

  // dump a single non-scalar attribute
  static int attrdump(DiscParticle &dp, vector<string> &str,
                      int r, Inform &msg, Inform &errmsg) {
    // based on the scalar type, invoke the proper dump routine
    // for a record written using just one attrib
    int atype = DiscTypeBase::scalarType(str[0]);
    if (atype == DiscTypeBase::CHAR) {
      return CreateAttribTypeDim<char, Dim>::attrdump(dp, str,
						      r, msg, errmsg);
    } else if (atype == DiscTypeBase::SHORT) {
      return CreateAttribTypeDim<short, Dim>::attrdump(dp, str,
						       r, msg, errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return CreateAttribTypeDim<int, Dim>::attrdump(dp, str,
						       r, msg, errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return CreateAttribTypeDim<long, Dim>::attrdump(dp, str,
						       r, msg, errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return CreateAttribTypeDim<float, Dim>::attrdump(dp, str,
						       r, msg, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return CreateAttribTypeDim<double, Dim>::attrdump(dp, str,
						       r, msg, errmsg);
    } else {
      errmsg << "In CreateAttribDim::attrdump: with str='";
      errmsg << str[0].c_str();
      errmsg << "', atype=" << atype << ", r=" << r << ": ";
      errmsg << "Unsupported attribute apptype scalar type." << endl;
      return RC_UNSUPPORTED;
    }
  }
};


////////////////////////////////////////////////////////////////////////////
// a function to perform the dump of an attribute, either scalar or non-scalar.
int perform_dump(ParticleAttribBase &attr, string &str, int r, int a,
		 Inform &msg, Inform &errmsg) {
  int dim = DiscTypeBase::dim(str);
  if (dim < 0) {
    errmsg << "Unknown attribute type '" << str << "'." << endl;
  } else if (dim == 0) {
    // this is a scalar
    int atype = DiscTypeBase::scalarType(str);
    if (atype == DiscTypeBase::CHAR) {
      return PerformDump<char>::dump(attr, r, a, msg, errmsg);
    } if (atype == DiscTypeBase::SHORT) {
      return PerformDump<short>::dump(attr, r, a, msg, errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return PerformDump<int>::dump(attr, r, a, msg, errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return PerformDump<long>::dump(attr, r, a, msg, errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return PerformDump<float>::dump(attr, r, a, msg, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return PerformDump<double>::dump(attr, r, a, msg, errmsg);
    } else {
      errmsg << "Unsupported attribute scalar type." << endl;
    }
  } else if (dim == 1) {
    return CreateAttribDim<1>::dump(attr, str, r, a, msg, errmsg);
  } else if (dim == 2) {
    return CreateAttribDim<2>::dump(attr, str, r, a, msg, errmsg);
  } else if (dim == 3) {
    return CreateAttribDim<3>::dump(attr, str, r, a, msg, errmsg);
  } else if (dim == 4) {
    return CreateAttribDim<4>::dump(attr, str, r, a, msg, errmsg);
  } else {
    errmsg << "Unsupported attribute type dimension " << dim << "." << endl;
  }

  errmsg << "Encountered in perform_dump: with str='" << str.c_str();
  errmsg << ", r=" << r << ", a=" << a << endl;
  return RC_UNSUPPORTED;
}


////////////////////////////////////////////////////////////////////////////
// a simple function to create the proper type of attribute based on
// the type string.  It determines the dimension and defers the next
// step to the scalar type routine
ParticleAttribBase *create_attrib(string &str, Inform &errmsg) {
  int dim = DiscTypeBase::dim(str);
  if (dim < 0) {
    errmsg << "Unknown attribute type '" << str << "'." << endl;
  } else if (dim == 0) {
    // this is a scalar
    int atype = DiscTypeBase::scalarType(str);
    if (atype == DiscTypeBase::CHAR) {
      return new ParticleAttrib<char>;
    } else if (atype == DiscTypeBase::SHORT) {
      return new ParticleAttrib<short>;
    } else if (atype == DiscTypeBase::INT) {
      return new ParticleAttrib<int>;
    } else if (atype == DiscTypeBase::LONG) {
      return new ParticleAttrib<long>;
    } else if (atype == DiscTypeBase::FLOAT) {
      return new ParticleAttrib<float>;
    } else if (atype == DiscTypeBase::DOUBLE) {
      return new ParticleAttrib<double>;
    } else {
      errmsg << "Unsupported attribute scalar type." << endl;
    }

  } else if (dim == 1) {
    return CreateAttribDim<1>::create(str, errmsg);
  } else if (dim == 2) {
    return CreateAttribDim<2>::create(str, errmsg);
  } else if (dim == 3) {
    return CreateAttribDim<3>::create(str, errmsg);
  } else if (dim == 4) {
    return CreateAttribDim<4>::create(str, errmsg);
  } else {
    errmsg << "Unsupported attribute type dimension " << dim << "." << endl;
  }

  errmsg << "Encountered in create_attrib: with str='" << str.c_str();
  errmsg << "'" << endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////////
// a struct used to do the actual dump work for a series of attributes,
// after the user has determined the data type for the position attribute.
template<class T, unsigned Dim>
struct DumpTypeDim {
  static int dump(DiscParticle &dp, vector<string> &str,
                  int r, int attrib, Inform &msg, Inform &errmsg) {

    // create the Particles object and the position attribute
    typedef ParticleUniformLayout<T, Dim>  Layout_t;
    typedef ParticleBase<Layout_t>         Particles_t;
    Particles_t particles(new Layout_t); // this already has R and ID attribs

    // go through and create the other attributes, and add them in
    int a;
    for (a = 2; a < dp.get_NumAttributes(r); ++a) {
      ParticleAttribBase *pa = create_attrib(str[a], errmsg);
      if (pa == 0)
	return RC_UNSUPPORTED;
      particles.addAttribute(*pa);
    }

    // now read in the particles for this record
    int retval = RC_OK;
    if (!dp.read(particles, r)) {
      errmsg << "Could not read data from file for record " << r << endl;
      retval = RC_READERR;
    }

    // and then go through and do the dumps
    if (retval == RC_OK) {
      for (a = 0; a < dp.get_NumAttributes(r); ++a) {
	if (attrib < 0 || a == attrib) {
	  retval = perform_dump(particles.getAttribute(a),
	 		        str[a], r, a, msg, errmsg);
	  if (retval != RC_OK)
	    break;
	}
      }
    }

    // delete all the attributes when we're done
    for (a = 2; a < dp.get_NumAttributes(r); ++a) {
      delete &(particles.getAttribute(a));
    }

    // and return
    return retval;
  }
};


////////////////////////////////////////////////////////////////////////////
// a simple struct used to determine, after the dimension has been selected,
// the type for the position attribute of a particles object
template<unsigned Dim>
struct DumpDim {
  // dump values for a whole ParticleBase
  static int dump(DiscParticle &dp, vector<string> &str,
                  int r, int attrib, Inform &msg, Inform &errmsg) {
    // based on the scalar type, invoke the proper dumping routine for
    // a record written using ALL attribs
    int atype = DiscTypeBase::scalarType(str[0]);
    if (atype == DiscTypeBase::FLOAT) {
      return DumpTypeDim<float, Dim>::dump(dp, str,
					   r, attrib, msg, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return DumpTypeDim<double, Dim>::dump(dp, str,
					    r, attrib, msg, errmsg);
    } else {
      errmsg << "Unsupported data type for position attribute in record ";
      errmsg << r << endl;
      return RC_UNSUPPORTED;
    }
  }
};


////////////////////////////////////////////////////////////////////////////
// Create and dump a single attribute from a file.  This is used if
// the data was written with just a single attribute instead of written
// as an entire ParticleBase.
int dump_attrib(DiscParticle &dp, vector<string> &str,
                int r, Inform &msg, Inform &errmsg)
{
  // get the dimension of the attrib
  // from the dim info, call the proper dumping routine
  int dim = DiscTypeBase::dim(str[0]);
  if (dim < 0) {
    errmsg << "Unknown attribute type '" << str[0] << "'." << endl;
  } else if (dim == 0) {
    // this is a scalar
    int atype = DiscTypeBase::scalarType(str[0]);
    if (atype == DiscTypeBase::CHAR) {
      return PerformDump<char>::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::SHORT) {
      return PerformDump<short>::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return PerformDump<int>::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return PerformDump<long>::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return PerformDump<float>::attrdump(dp, r, msg, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return PerformDump<double>::attrdump(dp, r, msg, errmsg);
    } else {
      errmsg << "Unsupported attribute scalar type." << endl;
    }
  } else if (dim == 1) {
    return CreateAttribDim<1>::attrdump(dp, str, r, msg, errmsg);
  } else if (dim == 2) {
    return CreateAttribDim<2>::attrdump(dp, str, r, msg, errmsg);
  } else if (dim == 3) {
    return CreateAttribDim<3>::attrdump(dp, str, r, msg, errmsg);
  } else if (dim == 4) {
    return CreateAttribDim<4>::attrdump(dp, str, r, msg, errmsg);
  } else {
    errmsg << "Unsupported attribute type dimension " << dim << "." << endl;
  }

  errmsg << "Encountered in dump_attrib: with str='" << str[0].c_str();
  errmsg << "', r=" << r << endl;
  return RC_UNSUPPORTED;
}


////////////////////////////////////////////////////////////////////////////
// select the proper type of Particles object to create and ask another
// routine to create and dump the attributes.
// This version is for the case where the record stores all the attributes,
// so that the first attribute gives the position attrib type
int dump_all(DiscParticle &dp, vector<string> &str,
             int r, int attrib, Inform &msg, Inform &errmsg)
{
  // there must be at least two attribs, the position and ID attribs
  if (dp.get_NumAttributes(r) < 2) {
    errmsg << "There must be at least two attributes in the records." << endl;
    return RC_READERR;
  }

  // get the dimension and scalar type of the first attrib, which must
  // be a vektor; also, the second attrib must be an integer
  int dim   = DiscTypeBase::dim(str[0]);
  int stype = DiscTypeBase::scalarType(str[0]);
  int atype = DiscTypeBase::appType(str[0]);

  if (atype != DiscTypeBase::VEKTOR || stype == DiscTypeBase::UNKNOWN) {
    errmsg << "The first attribute must be a Vektor of known type." << endl;
    return RC_READERR;
  }
  if (DiscTypeBase::scalarType(str[1]) != DiscTypeBase::INT ||
      DiscTypeBase::dim(str[1]) != 0) {
    errmsg << "The second attribute must be an integer." << endl;
    return RC_READERR;
  }

  // from the dim info, call the proper dumping routine
  if (dim == 1) {
    return DumpDim<1>::dump(dp, str, r, attrib, msg, errmsg);
  } else if (dim == 2) {
    return DumpDim<2>::dump(dp, str, r, attrib, msg, errmsg);
  } else if (dim == 3) {
    return DumpDim<3>::dump(dp, str, r, attrib, msg, errmsg);
  } else {
    errmsg << "Unsupported dimension " << dim << " for position attribute ";
    errmsg << " in record " << r << endl;
    return RC_UNSUPPORTED;
  }
}


////////////////////////////////////////////////////////////////////////////
// main routine
int main(int argc, char *argv[]) {
  unsigned int i;
  int retval = RC_CMDERR;
  vector<string> usertypes;
  vector<string> attribtypes;

  // scan early for -h flag
  for (i=1; i < argc; ++i)
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?"))
      usage(argc, argv, RC_OK);

  // initialize ippl
  Ippl ippl(argc,argv);
  Inform msg("dump", INFORM_ALL_NODES);
  Inform errmsg("dump error");

  // we need some kind of arguments
  if (argc == 1)
    usage(argc, argv, retval);

  // parse command-line options to get info
  char *config = 0;
  char *fname = 0;
  char *typestring = 0;
  int attribnum = (-1);
  int recordnum = (-1);

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
    } else if (!strcmp(argv[i], "-attrib")) {
      if (i < (argc - 1))
	attribnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-record")) {
      if (i < (argc - 1))
	recordnum = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-types")) {
      if (i < (argc - 1))
	typestring = argv[++i];
      else
	usage(argc, argv, retval);
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

  // break out the user-defined type string, if provided
  if (typestring != 0) {
    char *sptr = typestring;
    while (*sptr != 0) {
      // move to next name
      while (*sptr == '-') ++sptr;
      if (*sptr != 0) {
	// find length of current name, and copy into buffer
	char typebuf[4];
	int n = 0;
	while (*sptr != '-' && *sptr != 0) {
	  if (n > 2) {
	    errmsg << "Bad format for user-specified type string." << endl;
	    return RC_CMDERR;
	  }
	  typebuf[n++] = *sptr++;
	}

	// make this into a string
	typebuf[n] = 0;
	string u = typebuf;
	usertypes.push_back(u);
      }
    }
  }

  // open the DiscParticle file
  DiscParticle dp(fname, config, DiscParticle::INPUT);
  int numrec = dp.get_NumRecords();

  // make sure there are sufficient records for the requested dump
  if (numrec < 1 || numrec <= recordnum) {
    errmsg << "Not enough records in file '" << fname << "'." << endl;
    return RC_READERR;
  }

  // process the records, doing a dump for each requested record and
  // attribute.
  retval = RC_OK;
  for (int currrec = 0; currrec < numrec; ++currrec) {

    // if the current record is not the one we want for the first file,
    // skip it; also, if our return value is not OK, skip the record
    if (recordnum >= 0 && currrec != recordnum)
      continue;
    else if (retval != RC_OK)
      break;

    // make sure there are the proper number of attributes
    int numattr = dp.get_NumAttributes(currrec);
    if (numattr < 1 || numattr <= attribnum) {
      errmsg << "Not enough attributes in record " << currrec << " of file '";
      errmsg << fname << "'." << endl;
      return RC_READERR;
    }

    // go through and determine the type strings for the different attributes
    attribtypes.erase(attribtypes.begin(), attribtypes.end());
    for (int currattrib = 0; currattrib < numattr; ++currattrib) {
      bool typefound = false;
      if (dp.get_DiscType(currrec, 0) != 0) {
	string atype = dp.get_DiscType(currrec, currattrib);
        if (DiscTypeBase::scalarType(atype) != DiscTypeBase::UNKNOWN) {
	  attribtypes.push_back(atype);
          typefound = true;
        }
      }

      if (!typefound && currattrib < usertypes.size()) {
	attribtypes.push_back(usertypes[currattrib]);
        typefound = true;
      }

      if (!typefound) {
        errmsg << "You must provide a typestring with the -types argument, ";
        errmsg << "since the input file does not contain known type ";
        errmsg << "info in record " << currrec << " for attribute ";
        errmsg << currattrib << "." << endl;
        return RC_CMDERR;
      }
    }

    // perform the printing for this record and attribute(s)
    if (dp.get_DataMode(currrec) == DiscParticle::ALL) {
      retval = dump_all(dp, attribtypes, currrec, attribnum, msg, errmsg);
    } else {
      retval = dump_attrib(dp, attribtypes, currrec, msg, errmsg);
    }
  }

  return retval;
}


/***************************************************************************
 * $RCSfile: dcpdump.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcpdump.cpp,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/
