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
  dcpdiff compares the contents of two Particles files (that must have the
  same numbers and types of attributes, with the same sizes).  It reports
  the RMD's for each attribute, and the min and max diffs.

  Exit codes for dcpdiff:
    0 : OK, all attribs compared within specified tolerance or no tol given
    1 : Error in reading data files
    2 : Error in command-line parameters
    3 : Files read properly, but failed to compare within specified rmsd
    4 : Unsupported dimension or type
 ***************************************************************************/

// enumeration of possible error codes
enum { RC_OK, RC_READERR, RC_CMDERR, RC_RMSDERR, RC_UNSUPPORTED };

// a small hack: make the maximum allowed RMSD value a global.  If this is
// less than 0, no max RMSD value is used.
double maxrmsd = -1.0;

// if this flag is true, try to sort the data before diffing.  This ill
// not work in parallel.
bool sortdata1 = false;
bool sortascending1 = true;
int sortattrib1 = 1;		// index of attribute to sort
bool sortdata2 = false;
bool sortascending2 = true;
int sortattrib2 = 1;		// index of attribute to sort


////////////////////////////////////////////////////////////////////////////
// print out a simple usage message and exit
void usage(int, char *argv[], int retval) {
  Inform msg("Usage");
  msg << argv[0]<< " -file1 <file> -file2 <file> [-types <typestring>] [opts]";
  msg << "\n  where [opts] may be any of the following:";
  msg << "\n   -h or -? : print this help info and quit.";
  msg << "\n   -conf1 <file> : configuration file for file1";
  msg << "\n   -conf2 <file> : configuration file for file2";
  msg << "\n   -record <n> : which rec. to use from file (0 to numrecs-1)";
  msg << "\n   -attrib <n> : which attrib to compare (does all if not given)";
  msg << "\n   -maxrmsd <val> : maximum rmsd value. If the rmsd for a record";
  msg << "\n                    is >= maxrmsd, return error code value of 3";
  msg << "\n   -sort1 <n> <dir=1|0> : sort the attributes in the first file";
  msg << "\n      based on the values of the Nth attribute (0 to numattr-1)";
  msg << "\n      dir is 1 for ascending, 0 for descending";
  msg << "\n   -sort2 <n> <dir=1|0> : same as -sort1, for the second file";
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
  msg << "\n  typically something like 'v3d-i' for 3D Vektors and int's,";
  msg << "\n  respectively.";
  msg << "\n";
  msg << "\n  Possible error codes:";
  msg << "\n   0: OK, attribs compared within tolerance or no tolerance given";
  msg << "\n   1: Error in reading data files";
  msg << "\n   2: Error in command-line parameters";
  msg << "\n   3: Files read properly, but do not compare within given rmsd";
  msg << "\n   4: Unsupported dimension or type";
  msg << endl;
  exit(retval);
}


////////////////////////////////////////////////////////////////////////////
// simple function to compute field of doubles indicating elementwise diffs
// between two fields
template<class T>
struct CalcDiff {
  static void diff(ParticleAttrib<T> &a,
		   ParticleAttrib<T> &b,
		   ParticleAttrib<double> &c) {
    c = ((a-b) * (a-b));
  }
};

// partial specialization to Vektors
template<class T, unsigned Dim>
struct CalcDiff< Vektor<T,Dim> > {
  static void diff(ParticleAttrib< Vektor<T,Dim> > &a,
		   ParticleAttrib< Vektor<T,Dim> > &b,
		   ParticleAttrib<double> &c) {
    c = dot(a-b, a-b);
  }
};

// partial specialization to Tenzors
template<class T, unsigned Dim>
struct CalcDiff< Tenzor<T,Dim> > {
  static void diff(ParticleAttrib< Tenzor<T,Dim> > &a,
		   ParticleAttrib< Tenzor<T,Dim> > &b,
		   ParticleAttrib<double> &c) {
    c = dotdot(a-b, a-b);
  }
};

// partial specialization to SymTenzors
template<class T, unsigned Dim>
struct CalcDiff< SymTenzor<T,Dim> > {
  static void diff(ParticleAttrib< SymTenzor<T,Dim> > &a,
		   ParticleAttrib< SymTenzor<T,Dim> > &b,
		   ParticleAttrib<double> &c) {
    c = dotdot(a-b, a-b);
  }
};


////////////////////////////////////////////////////////////////////////////
// a simple struct used to invoke the proper CalcDiff function
template<class T>
struct PerformDiffs {
  static int diff(ParticleAttribBase &ba, ParticleAttribBase &bb,
		  int totalnum, int r1, int r2, int attrib, Inform &msg) {

    // cast the attributes to type T
    ParticleAttrib<T> &a = static_cast<ParticleAttrib<T> &>(ba);
    ParticleAttrib<T> &b = static_cast<ParticleAttrib<T> &>(bb);

    // calculate the diff
    ParticleAttrib<double> c;
    c.create(a.size());
    CalcDiff<T>::diff(a, b, c);

    // compute min, max, etc
    double mindiff = sqrt(min(c));
    double maxdiff = sqrt(max(c));
    double rmsdval = sum(c);
    rmsdval = sqrt(rmsdval / (double)totalnum);

    // print out results
    msg << "=" << r1 << ", R2=" << r2 << ", A=" << attrib;
    msg << " : mindiff = " << mindiff << "  maxdiff = " << maxdiff;
    msg << "  rmsd = " << rmsdval << endl;

    // check for tolerance failure ... set return code
    if (maxrmsd > 0.0 && rmsdval > maxrmsd)
      return RC_RMSDERR;
    else
      return RC_OK;
  }

  // compute differences for just two single attributes
  static int attrdiff(DiscParticle &dp1, DiscParticle &dp2,
		      int r1, int r2) {
    Inform msg("diff");
    Inform errmsg("diff error");

    // create two attributes of type T
    ParticleAttrib<T> a;
    ParticleAttrib<T> b;

    // read them in from the files
    int retval = RC_OK;
    if (!dp1.read(a, r1)) {
      errmsg << "Could not read attrib from first particles file for record ";
      errmsg << r1 << endl;
      retval = RC_READERR;
    }

    if (!dp2.read(b, r2)) {
      errmsg << "Could not read attrib from second particles file for record ";
      errmsg << r2 << endl;
      retval = RC_READERR;
    }

    // sort the data now, if necessary.
    ParticleAttribBase::SortList_t sortlist;
    if (sortdata1 || sortdata2) {
      // make sure we're not running parallel, this only works for now
      // on one node.
      if (Ippl::getNodes() > 1) {
	errmsg << "The -sort option cannot be used in parallel." << endl;
	retval = RC_READERR;
      }
    }

    if (sortdata1) {
      a.calcSortList(sortlist, sortascending1);
      a.sort(sortlist);
    }

    if (sortdata2) {
      b.calcSortList(sortlist, sortascending2);
      b.sort(sortlist);
    }

    // compute the total number of particles
    int n = a.size();
    reduce(n, n, OpAddAssign());

    // do the diffs, if necessary
    if (retval == RC_OK)
      retval = PerformDiffs<T>::diff(a, b, n, r1, r2, 0, msg);

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

  // diff two attributes
  static int diff(ParticleAttribBase &attr1, ParticleAttribBase &attr2,
		  string &str, int n, int r1, int r2, int a,
		  Inform &msg, Inform &errmsg) {
    int atype = DiscTypeBase::appType(str);
    if (atype == DiscTypeBase::VEKTOR) {
      return PerformDiffs<Vektor<T, Dim> >::diff(attr1,attr2,n,r1,r2,a,msg);
    } else if (atype == DiscTypeBase::TENZOR) {
      return PerformDiffs<Tenzor<T, Dim> >::diff(attr1,attr2,n,r1,r2,a,msg);
    } else if (atype == DiscTypeBase::SYMTENZOR) {
      return PerformDiffs<SymTenzor<T, Dim> >::diff(attr1,attr2,n,r1,r2,a,msg);
    } else {
      errmsg << "In CreateAttribTypeDim::diff: with str='" << str.c_str();
      errmsg << "', atype=" << atype << ", n=" << n << ", r1=" << r1;
      errmsg << ", r2=" << r2 << ", a=" << a << ": ";
      errmsg << "Unsupported attribute apptype." << endl;
      return RC_UNSUPPORTED;
    }
  }

  // compute differences for just two single attributes
  static int attrdiff(DiscParticle &dp1, DiscParticle &dp2,
		      vector<string> &str, int r1, int r2,
		      Inform &errmsg) {
    // based on the scalar type, invoke the proper diffing routine for
    // computing the difference of a record written using just one attrib
    int atype = DiscTypeBase::appType(str[0]);
    if (atype == DiscTypeBase::VEKTOR) {
      return PerformDiffs<Vektor<T, Dim> >::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::TENZOR) {
      return PerformDiffs<Tenzor<T, Dim> >::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::SYMTENZOR) {
      return PerformDiffs<SymTenzor<T, Dim> >::attrdiff(dp1, dp2, r1, r2);
    } else {
      errmsg << "In CreateAttribTypeDim::attrdiff: with str='";
      errmsg << str[0].c_str();
      errmsg << "', atype=" << atype << ", r1=" << r1;
      errmsg << ", r2=" << r2 << ": ";
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

  // diff two attributes
  static int diff(ParticleAttribBase &attr1, ParticleAttribBase &attr2,
		  string &str, int n, int r1, int r2, int a,
		  Inform &msg, Inform &errmsg) {
    int atype = DiscTypeBase::scalarType(str);
    if (atype == DiscTypeBase::CHAR) {
      return CreateAttribTypeDim<char, Dim>::diff(attr1, attr2, str,
						  n, r1, r2, a, msg, errmsg);
    } else if (atype == DiscTypeBase::SHORT) {
      return CreateAttribTypeDim<short, Dim>::diff(attr1, attr2, str,
						   n, r1, r2, a, msg, errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return CreateAttribTypeDim<int, Dim>::diff(attr1, attr2, str,
						 n, r1, r2, a, msg, errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return CreateAttribTypeDim<long, Dim>::diff(attr1, attr2, str,
						  n, r1, r2, a, msg, errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return CreateAttribTypeDim<float, Dim>::diff(attr1, attr2, str,
						   n, r1, r2, a, msg, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return CreateAttribTypeDim<double, Dim>::diff(attr1, attr2, str,
						    n, r1, r2, a, msg, errmsg);
    } else {
      errmsg << "In CreateAttribDim::diff: with str='" << str.c_str();
      errmsg << "', atype=" << atype << ", n=" << n << ", r1=" << r1;
      errmsg << ", r2=" << r2 << ", a=" << a << ": ";
      errmsg << "Unsupported attribute apptype scalar type." << endl;
      return RC_UNSUPPORTED;
    }
  }

  // compute differences for just two single attributes
  static int attrdiff(DiscParticle &dp1, DiscParticle &dp2,
		      vector<string> &str, int r1, int r2,
		      Inform &errmsg) {
    // based on the scalar type, invoke the proper diffing routine for
    // computing the difference of a record written using just one attrib
    int atype = DiscTypeBase::scalarType(str[0]);
    if (atype == DiscTypeBase::CHAR) {
      return CreateAttribTypeDim<char, Dim>::attrdiff(dp1, dp2, str,
						      r1, r2, errmsg);
    } else if (atype == DiscTypeBase::SHORT) {
      return CreateAttribTypeDim<short, Dim>::attrdiff(dp1, dp2, str,
						       r1, r2, errmsg);
    } else if (atype == DiscTypeBase::INT) {
      return CreateAttribTypeDim<int, Dim>::attrdiff(dp1, dp2, str,
						     r1, r2, errmsg);
    } else if (atype == DiscTypeBase::LONG) {
      return CreateAttribTypeDim<long, Dim>::attrdiff(dp1, dp2, str,
						      r1, r2, errmsg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return CreateAttribTypeDim<float, Dim>::attrdiff(dp1, dp2, str,
						       r1, r2, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return CreateAttribTypeDim<double, Dim>::attrdiff(dp1, dp2, str,
							r1, r2, errmsg);
    } else {
      errmsg << "In CreateAttribDim::attrdiff: with str='";
      errmsg << str[0].c_str();
      errmsg << "', atype=" << atype << ", r1=" << r1;
      errmsg << ", r2=" << r2 << ": ";
      errmsg << "Unsupported attribute apptype scalar type." << endl;
      return RC_UNSUPPORTED;
    }
  }
};


////////////////////////////////////////////////////////////////////////////
// a function to perform the diffs between two attributes that must be of
// the same type.  Returns the error code for the result.
int perform_diff(ParticleAttribBase &attr1, ParticleAttribBase &attr2,
		 string &str, int n, int r1, int r2, int a,
		 Inform &msg, Inform &errmsg) {
  int dim = DiscTypeBase::dim(str);
  if (dim < 0) {
    errmsg << "Unknown attribute type." << endl;
  } else if (dim == 0) {
    // this is a scalar
    int atype = DiscTypeBase::scalarType(str);
    if (atype == DiscTypeBase::CHAR) {
      return PerformDiffs<char>::diff(attr1, attr2, n, r1, r2, a, msg);
    } else if (atype == DiscTypeBase::SHORT) {
      return PerformDiffs<short>::diff(attr1, attr2, n, r1, r2, a, msg);
    } else if (atype == DiscTypeBase::INT) {
      return PerformDiffs<int>::diff(attr1, attr2, n, r1, r2, a, msg);
    } else if (atype == DiscTypeBase::LONG) {
      return PerformDiffs<long>::diff(attr1, attr2, n, r1, r2, a, msg);
    } else if (atype == DiscTypeBase::FLOAT) {
      return PerformDiffs<float>::diff(attr1, attr2, n, r1, r2, a, msg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return PerformDiffs<double>::diff(attr1, attr2, n, r1, r2, a, msg);
    } else {
      errmsg << "Unsupported attribute scalar type." << endl;
    }
  } else if (dim == 1) {
    return CreateAttribDim<1>::diff(attr1, attr2, str, n, r1, r2,a,msg,errmsg);
  } else if (dim == 2) {
    return CreateAttribDim<2>::diff(attr1, attr2, str, n, r1, r2,a,msg,errmsg);
  } else if (dim == 3) {
    return CreateAttribDim<3>::diff(attr1, attr2, str, n, r1, r2,a,msg,errmsg);
  } else if (dim == 4) {
    return CreateAttribDim<4>::diff(attr1, attr2, str, n, r1, r2,a,msg,errmsg);
  } else {
    errmsg << "Unsupported attribute type dimension " << dim << "." << endl;
  }

  errmsg << "Encountered in perform_diff: with str='" << str.c_str();
  errmsg << "', n=" << n << ", r1=" << r1;
  errmsg << ", r2=" << r2 << ", a=" << a << endl;
  return RC_UNSUPPORTED;
}


////////////////////////////////////////////////////////////////////////////
// a simple function to create the proper type of attribute based on
// the type string.  It determines the dimension and defers the next
// step to the scalar type routine
ParticleAttribBase *create_attrib(string &str, Inform &errmsg) {
  int dim = DiscTypeBase::dim(str);
  if (dim < 0) {
    errmsg << "Unknown attribute type." << endl;
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
// a struct used to do the actual diff work for a series of attributes,
// after the user has determined the data type for the position attribute.
template<class T, unsigned Dim>
struct DiffTypeDim {
  static int diff(DiscParticle &dp1, DiscParticle &dp2,
		  vector<string> &str, int r1, int r2, int attrib,
		  Inform &errmsg) {

    Inform msg("diff");

    // create the Particles object and the position attribute
    typedef ParticleUniformLayout<T, Dim>  Layout_t;
    typedef ParticleBase<Layout_t>         Particles_t;
    Particles_t particles1(new Layout_t); // this already has R and ID attribs
    Particles_t particles2(new Layout_t); // this already has R and ID attribs

    // go through and create the other attributes, and add them in
    int a;
    for (a = 2; a < dp1.get_NumAttributes(r1); ++a) {
      ParticleAttribBase *pa1 = create_attrib(str[a], errmsg);
      ParticleAttribBase *pa2 = create_attrib(str[a], errmsg);
      if (pa1 == 0 || pa2 == 0)
	return RC_UNSUPPORTED;
      particles1.addAttribute(*pa1);
      particles2.addAttribute(*pa2);
    }

    // now read in the particles for this record
    int retval = RC_OK;
    if (!dp1.read(particles1, r1)) {
      errmsg << "Could not read data from first particles file for record ";
      errmsg << r1 << endl;
      retval = RC_READERR;
    }

    if (!dp2.read(particles2, r2)) {
      errmsg << "Could not read data from second particles file for record ";
      errmsg << r2 << endl;
      retval = RC_READERR;
    }

    // sort the data now, if necessary.
    ParticleAttribBase::SortList_t sortlist;
    if (sortdata1 || sortdata2) {
      // make sure we're not running parallel, this only works for now
      // on one node.
      if (Ippl::getNodes() > 1) {
	errmsg << "The -sort option cannot be used in parallel." << endl;
	retval = RC_READERR;
      }
    }

    if (sortdata1) {
      particles1.getAttribute(sortattrib1).calcSortList(sortlist,
							sortascending1);
      particles1.sort(sortlist);
    }

    if (sortdata2) {
      particles2.getAttribute(sortattrib2).calcSortList(sortlist,
							sortascending2);
      particles2.sort(sortlist);
    }

    // and then go through and do the diffs
    if (retval == RC_OK) {
      for (a = 0; a < dp1.get_NumAttributes(r1); ++a) {
	if (attrib < 0 || a == attrib) {
	  int r = perform_diff(particles1.getAttribute(a),
			       particles2.getAttribute(a),
			       str[a], dp1.get_NumGlobalParticles(),
			       r1, r2, a, msg, errmsg);
	  if (r == RC_RMSDERR)
	    retval = r;
	  if (r != RC_OK && r != RC_RMSDERR)
	    break;
	}
      }
    }

    // delete all the attributes when we're done
    for (a = 2; a < dp1.get_NumAttributes(r1); ++a) {
      delete &(particles1.getAttribute(a));
      delete &(particles2.getAttribute(a));
    }

    // and return
    return retval;
  }
};


////////////////////////////////////////////////////////////////////////////
// a simple struct used to determine, after the dimension has been selected,
// the type for the position attribute of a particles object
template<unsigned Dim>
struct DiffDim {
  // compute differences for a whole ParticleBase
  static int diff(DiscParticle &dp1, DiscParticle &dp2,
		  vector<string> &str, int r1, int r2, int attrib,
		  Inform &errmsg) {
    // based on the scalar type, invoke the proper diffing routine for
    // computing the difference of a record written using ALL attribs
    int atype = DiscTypeBase::scalarType(str[0]);
    if (atype == DiscTypeBase::FLOAT) {
      return DiffTypeDim<float, Dim>::diff(dp1, dp2, str,
					   r1, r2, attrib, errmsg);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return DiffTypeDim<double, Dim>::diff(dp1, dp2, str,
					    r1, r2, attrib, errmsg);
    } else {
      errmsg << "Unsupported data type for position attribute in records ";
      errmsg << r1 << " and " << r2 << endl;
      return RC_UNSUPPORTED;
    }
  }
};


////////////////////////////////////////////////////////////////////////////
// Create and diff two single attributes from files.  This is used if
// the data was written with just a single attribute instead of written
// as an entire ParticleBase.
int diff_attrib(DiscParticle &dp1, DiscParticle &dp2,
		vector<string> &str, int r1, int r2, Inform &errmsg)
{
  // get the dimension of the attrib
  // from the dim info, call the proper diffing routine
  int dim = DiscTypeBase::dim(str[0]);
  if (dim < 0) {
    errmsg << "Unknown attribute type." << endl;
  } else if (dim == 0) {
    // this is a scalar
    int atype = DiscTypeBase::scalarType(str[0]);
    if (atype == DiscTypeBase::CHAR) {
      return PerformDiffs<char>::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::SHORT) {
      return PerformDiffs<short>::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::INT) {
      return PerformDiffs<int>::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::LONG) {
      return PerformDiffs<long>::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::FLOAT) {
      return PerformDiffs<float>::attrdiff(dp1, dp2, r1, r2);
    } else if (atype == DiscTypeBase::DOUBLE) {
      return PerformDiffs<double>::attrdiff(dp1, dp2, r1, r2);
    } else {
      errmsg << "Unsupported attribute scalar type." << endl;
    }
  } else if (dim == 1) {
    return CreateAttribDim<1>::attrdiff(dp1, dp2, str, r1, r2, errmsg);
  } else if (dim == 2) {
    return CreateAttribDim<2>::attrdiff(dp1, dp2, str, r1, r2, errmsg);
  } else if (dim == 3) {
    return CreateAttribDim<3>::attrdiff(dp1, dp2, str, r1, r2, errmsg);
  } else if (dim == 4) {
    return CreateAttribDim<4>::attrdiff(dp1, dp2, str, r1, r2, errmsg);
  } else {
    errmsg << "Unsupported attribute type dimension " << dim << "." << endl;
  }

  errmsg << "Encountered in diff_attrib: with str='" << str[0].c_str();
  errmsg << "', r1=" << r1 << ", r2=" << r2 << endl;
  return RC_UNSUPPORTED;
}


////////////////////////////////////////////////////////////////////////////
// select the proper type of Particles object to create and ask another
// routine to create and diff the attributes.
// This version is for the case where the record stores all the attributes,
// so that the first attribute gives the position attrib type
int diff_all(DiscParticle &dp1, DiscParticle &dp2,
	     vector<string> &str, int r1, int r2, int attrib,
	     Inform &errmsg)
{
  // there must be at least two attribs, the position and ID attribs
  if (dp1.get_NumAttributes(r1) < 2 ||
      dp2.get_NumAttributes(r2) < 2) {
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

  // from the dim info, call the proper diffing routine
  if (dim == 1) {
    return DiffDim<1>::diff(dp1, dp2, str, r1, r2, attrib, errmsg);
  } else if (dim == 2) {
    return DiffDim<2>::diff(dp1, dp2, str, r1, r2, attrib, errmsg);
  } else if (dim == 3) {
    return DiffDim<3>::diff(dp1, dp2, str, r1, r2, attrib, errmsg);
  } else {
    errmsg << "Unsupported dimension " << dim << " for position attribute ";
    errmsg << " in records " << r1 << " and " << r2 << endl;
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
  Inform errmsg("diff error");

  // we need some kind of arguments
  if (argc == 1)
    usage(argc, argv, retval);

  // parse command-line options to get info
  char *config1 = 0;
  char *config2 = 0;
  char *fname1 = 0;
  char *fname2 = 0;
  char *typestring = 0;
  int attribnum = (-1);
  int recordnum = (-1);
  int rec1 = (-1), rec2 = (-1);

  for (i=1; i < argc; ++i) {
    if (!strcmp(argv[i], "-conf1")) {
      if (i < (argc - 1))
	config1 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-file1")) {
      if (i < (argc - 1))
	fname1 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-conf2")) {
      if (i < (argc - 1))
	config2 = argv[++i];
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-file2")) {
      if (i < (argc - 1))
	fname2 = argv[++i];
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
    } else if (!strcmp(argv[i], "-record1")) {
      if (i < (argc - 1))
	rec1 = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-record2")) {
      if (i < (argc - 1))
	rec2 = atoi(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-maxrmsd")) {
      if (i < (argc - 1))
	maxrmsd = atof(argv[++i]);
      else
	usage(argc, argv, retval);
    } else if (!strcmp(argv[i], "-sort1")) {
      if (i < (argc - 2)) {
	sortdata1 = true;
	sortattrib1 = atoi(argv[++i]);
	sortascending1 = (atoi(argv[++i]) == 1);
      } else {
	usage(argc, argv, retval);
      }
    } else if (!strcmp(argv[i], "-sort2")) {
      if (i < (argc - 2)) {
	sortdata2 = true;
	sortattrib2 = atoi(argv[++i]);
	sortascending2 = (atoi(argv[++i]) == 1);
      } else {
	usage(argc, argv, retval);
      }
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
  if (fname1 == 0 || fname2 == 0) {
    errmsg << "You must provide two base file names." << endl;
    usage(argc, argv, retval);
  }

  // if the user gave both -record and -record1, 2 opts, it is an error
  if (recordnum >= 0 && (rec1 >= 0 || rec2 >= 0)) {
    errmsg << "You cannot specify both -record and -record1 or -record2.";
    errmsg << endl;
    usage(argc, argv, retval);
  }

  // if the user gave a -record1 but not -record (or vice versa), error
  if ((rec1 >= 0 && rec2 < 0) || (rec2 >= 0 && rec1 < 0)) {
    errmsg << "You must specify BOTH -record1 and -record2." << endl;
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

  // open the two DiscParticle files
  DiscParticle dp1(fname1, config1, DiscParticle::INPUT);
  DiscParticle dp2(fname2, config2, DiscParticle::INPUT);
  int numrec1 = dp1.get_NumRecords();
  int numrec2 = dp2.get_NumRecords();

  // make sure there are sufficient records for the requested difference
  if (numrec1 < 1 || numrec1 <= recordnum || numrec1 <= rec1) {
    errmsg << "Not enough records in file '" << fname1 << "'." << endl;
    return RC_READERR;
  }
  if (numrec2 < 1 || numrec2 <= recordnum || numrec2 <= rec2) {
    errmsg << "Not enough records in file '" << fname2 << "'." << endl;
    return RC_READERR;
  }
  if (recordnum < 0 && rec1 < 0 && rec2 < 0 && numrec1 != numrec2) {
    errmsg << "The comparison files do not have the same number of record.";
    errmsg << endl;
    return RC_READERR;
  }

  // process the records, doing a diff for each requested record and
  // attribute.
  retval = RC_OK;
  for (int currrec = 0; currrec < numrec1; ++currrec) {

    // the actual record numbers for the two files that we'll use
    int r1 = currrec;
    int r2 = currrec;
    if (recordnum >= 0) {
      r1 = r2 = recordnum;
    } else if (rec1 >= 0 && rec2 >= 0) {
      r1 = rec1;
      r2 = rec2;
    }

    // if the current record is not the one we want for the first file,
    // skip it; also, if our return value is not OK, skip the record
    if (currrec != r1)
      continue;
    else if (retval != RC_OK)
      break;

    // make sure there are the proper number of attributes
    int numattr1 = dp1.get_NumAttributes(r1);
    int numattr2 = dp2.get_NumAttributes(r2);
    if (numattr1 < 1 || numattr1 <= attribnum) {
      errmsg << "Not enough attributes in record " << r1 << " of file '";
      errmsg << fname1 << "'." << endl;
      return RC_READERR;
    }
    if (numattr2 < 1 || numattr2 <= attribnum) {
      errmsg << "Not enough attributes in record " << r2 << " of file '";
      errmsg << fname2 << "'." << endl;
      return RC_READERR;
    }
    if (attribnum < 0 && numattr1 != numattr2) {
      errmsg << "The comparison files do not have the same number of ";
      errmsg << "attributes in record " << r1 << " and " << r2 << "." << endl;
      return RC_READERR;
    }
    if (dp1.get_DataMode(r1) != dp2.get_DataMode(r2)) {
      errmsg << "Files have different data output modes in record " << currrec;
      errmsg << "." << endl;
      return RC_READERR;
    }

    // go through and determine the type strings for the different attributes
    attribtypes.erase(attribtypes.begin(), attribtypes.end());
    for (int currattrib = 0; currattrib < numattr1; ++currattrib) {
      bool typefound = false;
      if (dp1.get_DiscType(r1, 0) != 0 && dp2.get_DiscType(r2, 0) != 0) {
	string atype1 = dp1.get_DiscType(r1, currattrib);
	string atype2 = dp2.get_DiscType(r2, currattrib);
	if (DiscTypeBase::dim(atype1) != DiscTypeBase::dim(atype2)) {
	  errmsg << "Dimensions for attributes " << currattrib<<" in records ";
	  errmsg << r1 << " and " << r2 << " differ." << endl;
	  return RC_READERR;
	} else if ((DiscTypeBase::scalarType(atype1) !=
		    DiscTypeBase::scalarType(atype2)) ||
		   (DiscTypeBase::appType(atype1) !=
		    DiscTypeBase::appType(atype2))) {
	  errmsg << "Data types for attributes " << currattrib<<" in records ";
	  errmsg << r1 << " and " << r2 << " differ." << endl;
	  return RC_READERR;
	} else if (DiscTypeBase::scalarType(atype1) != DiscTypeBase::UNKNOWN) {
	  attribtypes.push_back(atype1);
          typefound = true;
        }
      }

      if (!typefound && currattrib < usertypes.size()) {
	attribtypes.push_back(usertypes[currattrib]);
        typefound = true;
      }

      if (!typefound) {
        errmsg <<"You must provide a typestring with the -types argument, ";
        errmsg <<"since one or both of the input file does not contain known ";
        errmsg <<"type info in record " << currrec << " for attribute ";
        errmsg <<currattrib << "." << endl;
        return RC_CMDERR;
      }
    }

    // perform the differencing for this record and attribute(s)
    if (dp1.get_DataMode(r1) == DiscParticle::ALL) {
      retval = diff_all(dp1, dp2, attribtypes, r1, r2, attribnum, errmsg);
    } else {
      retval = diff_attrib(dp1, dp2, attribtypes, r1, r2, errmsg);
    }
  }

  return retval;
}


/***************************************************************************
 * $RCSfile: dcpdiff.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcpdiff.cpp,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/
