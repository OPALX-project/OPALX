/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCFDUMP_TEMPLATE_H
#define DCFDUMP_TEMPLATE_H

#include "Ippl.h"
#include "dcfdump.h"

/***************************************************************************
  dcfdump is a simple utility program which prints out
  the contents of a 1, 2, or 3-dimensional field of doubles
  or Vektors of doubles (can also print fields of floats or ints).
 ***************************************************************************/

////////////////////////////////////////////////////////////////////////////
// a function to actually print out the data for a field.
template<class T, unsigned Dim, class MT>
struct DumpData {
  static int dump(const char *fname1, const char *config1,
	  	  int fieldnum, int recordnum, Inform &errmsg)
  {
    Inform msg(fname1);

    // create DiscFields to read data, and check dimensions
    DiscField<Dim> io1(fname1, config1);
    int numRecords         = io1.get_NumRecords();
    int numFields          = io1.get_NumFields();
    NDIndex<Dim> domain    = io1.get_Domain();
    const char *typestring = io1.get_TypeString();

    // make sure everything is kosher
    if (numRecords < 1 && io1.get_Dimension() != Dim) {
      errmsg << "Bad dimension for file " << fname1 << endl;
      return RC_CMDERR;
    } else if (recordnum >= numRecords || numRecords < 1) {
      errmsg << "Not enough records in file " << fname1 << endl;
      return RC_CMDERR;
    } else if (fieldnum >= numFields || numFields < 1) {
      errmsg << "Not enough fields in file " << fname1 << endl;
      return RC_CMDERR;
    }

    // print out simple summary of contents
    msg << "Summary of data in file '" << fname1 << "':" << endl;
    msg << "        numRecords = " << numRecords << endl;
    msg << "  numFields/record = " << numFields << endl;
    msg << "      total domain = " << domain << endl;

    // initialize the field ... create an initial layout from the domain
    FieldLayout<Dim>                           layout(domain);
    Field<T,Dim,UniformCartesian<Dim,MT>,Cell> X1(layout);

    // set up for connection
    bool connected = false;
    DataConnect *conn = 0;

    // now read in the data, and compare it
    for (int j=0; j < numRecords; ++j) {
      if ((recordnum < 0 || j == recordnum)) {
	// get all field data in this record
	for (int i=0; i < numFields; ++i) {
	  if ((fieldnum < 0 || i==fieldnum)) {
	    // read the next record for this field
	    if (!io1.read(X1, i, j))
	      return RC_READERR;

	    // register the Field for displaying data, if this is the
	    // first read
	    if (!connected) {
	      conn = DataConnectCreator::create(fname1);
	      conn->connect(fname1, X1, DataSource::OUTPUT);
	      connected = true;
	    }
	    if (conn != 0) {
	      conn->updateConnections();
	      conn->updateConnections();
	      conn->interact();
	    }

	    // only print data if necessary
	    if (doprint) {
	      msg << "\nData for Field " << i << ", Record " << j << ":"<<endl;
	      msg << "------------------------------------------" << endl;
	      FieldPrint<T,Dim> printField(X1);
	      printField.print(domain);
	    }
	  }
	}
      }
    }

    // return with exit status
    return RC_OK;
  }
};

#endif

/***************************************************************************
 * $RCSfile: dcfdump_template.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcfdump_template.h,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/
