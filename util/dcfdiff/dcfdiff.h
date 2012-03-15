/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCFDIFF_H
#define DCFDIFF_H

// forward references
class Inform;


// enumeration of possible error codes
enum { RC_OK, RC_READERR, RC_CMDERR, RC_RMSDERR, RC_UNSUPPORTED };

// enumeration of possible diff types
enum { DIFFSCALAR, DIFFVECTOR, DIFFTENSOR, DIFFSYMTENSOR, DIFFANTISYMTENSOR };

// enumeration of possible diff scalar type
enum { DIFFDOUBLE, DIFFFLOAT, DIFFINT };

// global variables
extern bool printMaxDiff;


//---------------------------------------------------------------------------
// macros used to define prototypes and functions for invoking the templated
// generic diff routine
//---------------------------------------------------------------------------

#define DIFF_PROTOTYPE(FUNCTIONNAME)			\
extern int FUNCTIONNAME(int atype,			\
			const char *f1, const char *c1,	\
                        const char *f2, const char *c2,	\
			int fnum, int rnum,		\
                        double rmsd, Inform &errmsg);

#define DIFF_FUNCTION(FUNCTIONNAME, TYPE, DIM)				     \
int FUNCTIONNAME(int atype,						     \
		 const char *f1, const char *c1,			     \
                 const char *f2, const char *c2,			     \
		 int fnum, int rnum,					     \
                 double rmsd, Inform &errmsg)				     \
{									     \
  int retval = RC_UNSUPPORTED;						     \
									     \
  if (atype == DIFFSCALAR)						     \
    retval = DiffData<TYPE, DIM, TYPE>::diff(f1, c1, f2, c2,		     \
					     fnum, rnum, rmsd, errmsg);	     \
  else if (atype == DIFFVECTOR)						     \
    retval = DiffData<Vektor<TYPE, DIM>, DIM, TYPE>::diff(f1, c1, f2, c2,    \
                                             fnum, rnum, rmsd, errmsg);	     \
  else if (atype == DIFFTENSOR)						     \
    retval = DiffData<Tenzor<TYPE, DIM>, DIM, TYPE>::diff(f1, c1, f2, c2,    \
                                             fnum, rnum, rmsd, errmsg);	     \
  else if (atype == DIFFSYMTENSOR)					     \
    retval = DiffData<SymTenzor<TYPE, DIM>, DIM, TYPE>::diff(f1, c1, f2, c2, \
                                             fnum, rnum, rmsd, errmsg);	     \
  else									     \
    errmsg << "Unsupported user type." << endl;				     \
									     \
  return retval;							     \
}


//---------------------------------------------------------------------------
// prototypes for diff routines; these are defined in the instantiation files
//---------------------------------------------------------------------------

DIFF_PROTOTYPE(diff_data_d_1)
DIFF_PROTOTYPE(diff_data_d_2)
DIFF_PROTOTYPE(diff_data_d_3)
DIFF_PROTOTYPE(diff_data_d_4)

DIFF_PROTOTYPE(diff_data_f_1)
DIFF_PROTOTYPE(diff_data_f_2)
DIFF_PROTOTYPE(diff_data_f_3)
DIFF_PROTOTYPE(diff_data_f_4)

DIFF_PROTOTYPE(diff_data_i_1)
DIFF_PROTOTYPE(diff_data_i_2)
DIFF_PROTOTYPE(diff_data_i_3)
DIFF_PROTOTYPE(diff_data_i_4)

#endif

/***************************************************************************
 * $RCSfile: dcfdiff.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcfdiff.h,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/

