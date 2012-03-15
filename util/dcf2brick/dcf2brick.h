/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCF2BRICK_H
#define DCF2BRICK_H

// forward references
class Inform;


// enumeration of possible error codes
enum { RC_OK, RC_READERR, RC_WRITEERR, RC_CMDERR, RC_UNSUPPORTED };

// enumeration of possible data types
enum { CONVSCALAR, CONVVECTOR, CONVTENSOR, CONVSYMTENSOR, CONVANTISYMTENSOR };

// enumeration of possible data scalar type
enum { CONVDOUBLE, CONVFLOAT, CONVINT };


//---------------------------------------------------------------------------
// macros used to define prototypes and functions for invoking the templated
// generic conversion routine
//---------------------------------------------------------------------------

#define CONV_PROTOTYPE(FUNCTIONNAME)				\
extern int FUNCTIONNAME(int atype,				\
                        const char *f1, const char *c1,		\
			const char *f2,				\
                        int fnum,				\
			int startrnum, int endrnum,		\
			bool globalminmax,			\
			int minsize,				\
			bool dosum,				\
                        Inform &errmsg);

#define CONV_FUNCTION(FUNCTIONNAME, TYPE, DIM)				  \
int FUNCTIONNAME(int atype,						  \
		 const char *f1, const char *c1,			  \
		 const char *f2,					  \
		 int fnum,						  \
		 int startrnum, int endrnum,				  \
		 bool globalminmax,					  \
		 int minsize,						  \
		 bool dosum,						  \
		 Inform &errmsg)					  \
{									  \
  int retval = RC_UNSUPPORTED;						  \
									  \
  if (atype == CONVSCALAR)						  \
    retval = ConvData<TYPE, DIM, TYPE>::conv(f1, c1, f2, fnum,		  \
					     startrnum, endrnum,	  \
					     globalminmax, minsize,	  \
					     dosum, errmsg);		  \
									  \
  else if (atype == CONVVECTOR)						  \
    retval = ConvData<Vektor<TYPE,DIM>, DIM, TYPE>::conv(f1, c1, f2,	  \
							 fnum,		  \
							 startrnum,	  \
							 endrnum,	  \
							 globalminmax,	  \
							 minsize,	  \
							 dosum,		  \
							 errmsg);	  \
									  \
  else if (atype == CONVTENSOR)						  \
    retval = ConvData<Tenzor<TYPE,DIM>, DIM, TYPE>::conv(f1, c1, f2,	  \
							 fnum,		  \
							 startrnum,	  \
							 endrnum,	  \
							 globalminmax,	  \
							 minsize,	  \
							 dosum,		  \
							 errmsg);	  \
									  \
  else if (atype == CONVSYMTENSOR)					  \
    retval = ConvData<SymTenzor<TYPE,DIM>, DIM, TYPE>::conv(f1, c1, f2,	  \
							    fnum,	  \
							    startrnum,	  \
							    endrnum,	  \
							    globalminmax, \
							    minsize,	  \
							    dosum,	  \
							    errmsg);	  \
  else if (atype == CONVANTISYMTENSOR)					  \
    retval = ConvData<AntiSymTenzor<TYPE,DIM>, DIM, TYPE>::conv(f1, c1,   \
							    f2,	          \
							    fnum,	  \
							    startrnum,	  \
							    endrnum,	  \
							    globalminmax, \
							    minsize,	  \
							    dosum,	  \
							    errmsg);	  \
  else									  \
    errmsg << "Unsupported user type." << endl;				  \
									  \
  return retval;							  \
}


//---------------------------------------------------------------------------
// prototypes for dump routines; these are defined in the instantiation files
//---------------------------------------------------------------------------

CONV_PROTOTYPE(conv_data_d_1)
CONV_PROTOTYPE(conv_data_d_2)
CONV_PROTOTYPE(conv_data_d_3)
CONV_PROTOTYPE(conv_data_d_4)

CONV_PROTOTYPE(conv_data_f_1)
CONV_PROTOTYPE(conv_data_f_2)
CONV_PROTOTYPE(conv_data_f_3)
CONV_PROTOTYPE(conv_data_f_4)

CONV_PROTOTYPE(conv_data_i_1)
CONV_PROTOTYPE(conv_data_i_2)
CONV_PROTOTYPE(conv_data_i_3)
CONV_PROTOTYPE(conv_data_i_4)

#endif

/***************************************************************************
 * $RCSfile: dcf2brick.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:38 $
 * IPPL_VERSION_ID: $Id: dcf2brick.h,v 1.1.1.1 2003/01/23 07:41:38 adelmann Exp $ 
 ***************************************************************************/

