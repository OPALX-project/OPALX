/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCFDIFF_TEMPLATE_H
#define DCFDIFF_TEMPLATE_H

/***************************************************************************
  dcfdiff compares the contents of two fields (which must have the same
  dimensions and total domains), and computes RMSD's, max diffs, min diffs,
  and total diffs.

  this file contains the templated code for computing diffs.  Other files
  instantiate specific versions of this code, and provide the main function.
 ***************************************************************************/

// include files
#include "Ippl.h"
#include "dcfdiff.h"


////////////////////////////////////////////////////////////////////////////
// simple function to compute field of doubles indicating elementwise diffs
// between two fields
template<class T, unsigned Dim, class M, class C>
struct CalcDiff {
  static void diff(Field<T,Dim,M,C>& a,
		   Field<T,Dim,M,C>& b,
		   Field<double,Dim>& c) {
    c = (a-b) * (a-b);
  }
};

// partial specialization to Vektors
template<class T, unsigned Dim, class M, class C>
struct CalcDiff<Vektor<T,Dim>, Dim, M, C> {
  static void diff(Field<Vektor<T,Dim>,Dim,M,C>& a,
		   Field<Vektor<T,Dim>,Dim,M,C>& b,
		   Field<double,Dim>& c) {
    c = dot(a-b, a-b);
  }
};

// partial specialization to Tenzors
template<class T, unsigned Dim, class M, class C>
struct CalcDiff<Tenzor<T,Dim>, Dim, M, C> {
  static void diff(Field<Tenzor<T,Dim>,Dim,M,C>& a,
		   Field<Tenzor<T,Dim>,Dim,M,C>& b,
		   Field<double,Dim>& c) {
    c = dotdot(a-b, a-b);
  }
};

// partial specialization to SymTenzors
template<class T, unsigned Dim, class M, class C>
struct CalcDiff<SymTenzor<T,Dim>, Dim, M, C> {
  static void diff(Field<SymTenzor<T,Dim>,Dim,M,C>& a,
		   Field<SymTenzor<T,Dim>,Dim,M,C>& b,
		   Field<double,Dim>& c) {
    c = dotdot(a-b, a-b);
  }
};

// partial specialization to AntiSymTenzors
template<class T, unsigned Dim, class M, class C>
struct CalcDiff<AntiSymTenzor<T,Dim>, Dim, M, C> {
  static void diff(Field<AntiSymTenzor<T,Dim>,Dim,M,C>& a,
		   Field<AntiSymTenzor<T,Dim>,Dim,M,C>& b,
		   Field<double,Dim>& c) {
    c = dotdot(a-b, a-b);
  }
};


////////////////////////////////////////////////////////////////////////////
// a function to actually print out the data for a field.
// return error code (listed above).
template<class T, unsigned Dim, class MT>
struct DiffData {
  static int diff(const char *fname1, const char *config1,
		  const char *fname2, const char *config2,
		  int fieldnum, int recordnum, double maxrmsd,
                  Inform &errmsg)
  {
    Inform msg("diff");

    // create DiscFields to read data, and check dimensions
    DiscField<Dim> io1(fname1, config1);
    DiscField<Dim> io2(fname2, config2);
    int numRecs1         = io1.get_NumRecords();
    int numRecs2         = io2.get_NumRecords();
    int numFlds1         = io1.get_NumFields();
    int numFlds2         = io2.get_NumFields();
    NDIndex<Dim> domain1 = io1.get_Domain();
    NDIndex<Dim> domain2 = io2.get_Domain();

    // make sure everything is kosher
    if (numRecs1 < 1 && io1.get_Dimension() != Dim) {
      errmsg << "Bad dimension for file " << fname1 << endl;
      return RC_CMDERR;
    } else if (numRecs2 < 1 && io2.get_Dimension() != Dim) {
      errmsg << "Bad dimension for file " << fname2 << endl;
      return RC_CMDERR;
    } else if (recordnum >= numRecs1 || numRecs1 < 1) {
      errmsg << "Not enough records in file " << fname1 << endl;
      return RC_CMDERR;
    } else if (recordnum >= numRecs2 || numRecs2 < 1) {
      errmsg << "Not enough records in file " << fname2 << endl;
      return RC_CMDERR;
    } else if (fieldnum >= numFlds1 || numFlds1 < 1) {
      errmsg << "Not enough fields in file " << fname1 << endl;
      return RC_CMDERR;
    } else if (fieldnum >= numFlds2 || numFlds2 < 1) {
      errmsg << "Not enough fields in file " << fname2 << endl;
      return RC_CMDERR;
    } else if (!(domain1 == domain2)) {
      errmsg << "Both fields must have the same domain." << endl;
      return RC_CMDERR;
    }

    // initialize the field ... create an initial layout from the domain
    FieldLayout<Dim>                           layout(domain1);
    Field<double,Dim>                          results(layout);
    Field<T,Dim,UniformCartesian<Dim,MT>,Cell> X1(layout);
    Field<T,Dim,UniformCartesian<Dim,MT>,Cell> X2(layout);

    // now read in the data, and compare it
    int numRecords = (numRecs1 < numRecs2 ? numRecs1 : numRecs2);
    int numFields  = (numFlds1 < numFlds2 ? numFlds1 : numFlds2);
    int returncode = 0;
    for (int j=0; j < numRecords; ++j) {
      if ((recordnum < 0 || j == recordnum)) {
	for (int i=0; i < numFields; ++i) {
	  if ((fieldnum < 0 || i==fieldnum)) {
	    // read the next record for this field
	    if (!io1.read(X1, i, j))
	      return RC_READERR;
	    if (!io2.read(X2, i, j))
	      return RC_READERR;

	    // compute the diff^2 between the two; store in X1
	    CalcDiff<T,Dim,UniformCartesian<Dim,MT>,Cell>::diff(X1,X2,results);

	    // compute min, max, RMSD
            NDIndex<Dim> mloc;
	    double mindiff = sqrt(min(results));
	    double maxdiff = sqrt(max(results, mloc));
	    double rmsdval = sqrt(sum(results) / ((double)(domain1.size())));

	    // print out results
	    msg << "F=" << i << ",R=" << j << ": mindiff = " << mindiff;
	    msg << "  maxdiff = " << maxdiff << "  rmsd = " << rmsdval;

            if (printMaxDiff) {
              msg << "  maxdiffloc = {";
	      for (int mdi = 0; mdi < Dim; ++mdi) {
		msg << mloc[mdi].first();
		if (mdi < (Dim - 1))
		  msg << ",";
	      }
	      msg << "}";
	    }

            msg << endl;

	    // check for tolerance failure ... set return code
	    if (maxrmsd > 0.0 && rmsdval > maxrmsd)
	      returncode = RC_RMSDERR;
	  }
	}
      }
    }

    // return with exit status
    return returncode;
  }
};

#endif

/***************************************************************************
 * $RCSfile: dcfdiff_template.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcfdiff_template.h,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/
