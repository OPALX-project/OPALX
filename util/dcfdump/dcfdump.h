/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCFDUMP_H
#define DCFDUMP_H

// forward references
class Inform;


// enumeration of possible error codes
enum { RC_OK, RC_READERR, RC_CMDERR, RC_UNSUPPORTED };

// enumeration of possible dump types
enum { DUMPSCALAR, DUMPVECTOR, DUMPTENSOR, DUMPSYMTENSOR, DUMPANTISYMTENSOR };

// enumeration of possible dump scalar type
enum { DUMPDOUBLE, DUMPFLOAT, DUMPINT };

// global external variables
extern bool doprint;
extern bool dumpoffset;
extern bool dumplayout;


//---------------------------------------------------------------------------
// macros used to define prototypes and functions for invoking the templated
// generic diff routine
//---------------------------------------------------------------------------

#define DUMP_PROTOTYPE(FUNCTIONNAME)			\
extern int FUNCTIONNAME(int atype,			\
                        const char *f1, const char *c1,	\
                        int fnum, int rnum,		\
                        Inform &errmsg);

#define DUMP_FUNCTION(FUNCTIONNAME, TYPE, DIM)				   \
int FUNCTIONNAME(int atype,						   \
		 const char *f1, const char *c1,			   \
		 int fnum, int rnum,					   \
		 Inform &errmsg)					   \
{									   \
  int retval = RC_UNSUPPORTED;						   \
									   \
  if (atype == DUMPSCALAR)						   \
    retval = DumpData<TYPE, DIM, TYPE>::dump(f1, c1, fnum, rnum, errmsg);  \
									   \
  else if (atype == DUMPVECTOR)						   \
    retval = DumpData<Vektor<TYPE,DIM>, DIM, TYPE>::dump(f1, c1, fnum,	   \
							 rnum, errmsg);	   \
									   \
  else if (atype == DUMPTENSOR)						   \
    retval = DumpData<Tenzor<TYPE,DIM>, DIM, TYPE>::dump(f1, c1, fnum,	   \
							 rnum, errmsg);	   \
									   \
  else if (atype == DUMPSYMTENSOR)					   \
    retval = DumpData<SymTenzor<TYPE,DIM>, DIM, TYPE>::dump(f1, c1, fnum,  \
							    rnum, errmsg); \
									   \
  else if (atype == DUMPANTISYMTENSOR)					   \
    retval = DumpData<AntiSymTenzor<TYPE,DIM>, DIM, TYPE>::dump(f1, c1,	   \
								fnum,	   \
								rnum,	   \
								errmsg);   \
  else									   \
    errmsg << "Unsupported user type." << endl;				   \
									   \
  return retval;							   \
}


//---------------------------------------------------------------------------
// prototypes for dump routines; these are defined in the instantiation files
//---------------------------------------------------------------------------

DUMP_PROTOTYPE(dump_data_d_1)
DUMP_PROTOTYPE(dump_data_d_2)
DUMP_PROTOTYPE(dump_data_d_3)
DUMP_PROTOTYPE(dump_data_d_4)

DUMP_PROTOTYPE(dump_data_f_1)
DUMP_PROTOTYPE(dump_data_f_2)
DUMP_PROTOTYPE(dump_data_f_3)
DUMP_PROTOTYPE(dump_data_f_4)

DUMP_PROTOTYPE(dump_data_i_1)
DUMP_PROTOTYPE(dump_data_i_2)
DUMP_PROTOTYPE(dump_data_i_3)
DUMP_PROTOTYPE(dump_data_i_4)

#endif

/***************************************************************************
 * $RCSfile: dcfdump.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcfdump.h,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/

