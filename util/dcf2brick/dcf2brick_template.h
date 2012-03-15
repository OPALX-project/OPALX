/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCF2BRICK_TEMPLATE_H
#define DCF2BRICK_TEMPLATE_H

#include "Ippl.h"
#include "dcf2brick.h"
#include <iostream.h>
#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////
// a function to convert from a floating point field of scalars or
// vectors to a field of bytes ... values will be remapped to range
// 0 ... 255.
// This is the general case that will work for scalars.
template<class T, unsigned Dim>
struct ConvertField
{
  // calculate the minimum and maximum for a field of scalars
  template<class F1>
  static void minmax(F1 &a, T &rminval, T &rmaxval, bool firsttime) {
    // find min and max of a
    T minval = min(a);
    T maxval = max(a);

    // change returned values, if necessary
    if (firsttime || minval < rminval)
      rminval = minval;
    if (firsttime || maxval > rmaxval)
      rmaxval = maxval;
  }

  // the conversion function.  F1 should be a field of the floating
  // point values, F2 a field of bytes with the same layout.
  template<class F1, class F2>
  static void convert(F1 &a, F2 &b, T minval, T maxval) {
    // convert to range 0 ... 255, and assign to bytes
    if (minval == maxval) {
      // everything is the same, so just make everything in b be zero
      b = 0;
    } else {
      // rescale the range
      b = (a - minval)/(maxval - minval) * 255.0;
    }
  }

  // accumulate the first field into the second
  template<class F1, class F2>
  static void accum(F1 &X2, F2 &X1, int k, int numFields) {
    X1 = where(X2 > 0.99, k+1, X1);
    X1 = where(X2 <= 0.99 && X2 > 0.0, numFields + 1, X1);
  }
};


////////////////////////////////////////////////////////////////////////////
// a function to convert from a floating point field of scalars or
// vectors to a field of bytes ... values will be remapped to range
// 0 ... 255.
// This is a specialization for Vektors
template<class T, unsigned Dim>
struct ConvertField<Vektor<T,Dim>, Dim>
{
  // calculate the minimum and maximum for a field of vektors
  template<class F1>
  static void minmax(F1 &a, T &rminval, T &rmaxval, bool firsttime) {
    // find min and max of a
    T minval = min(dot(a,a));
    T maxval = max(dot(a,a));

    // change returned values, if necessary
    if (firsttime || minval < rminval)
      rminval = minval;
    if (firsttime || maxval > rmaxval)
      rmaxval = maxval;
  }

  // the conversion function.  F1 should be a field of the vektor
  // values, F2 a field of bytes with the same layout.
  template<class F1, class F2>
  static void convert(F1 &a, F2 &b, T minval, T maxval) {
    // convert to range 0 ... 255, and assign to bytes
    if (minval == maxval) {
      // everything is the same, so just make everything in b be zero
      b = 0;
    } else {
      // rescale the range
      minval = sqrt(minval);
      maxval = sqrt(maxval);
      b = (sqrt(dot(a,a)) - minval)/(maxval - minval) * 255.0;
    }
  }

  // accumulate the first field into the second
  template<class F1, class F2>
  static void accum(F1 &X2, F2 &X1, int, int) {
    X1 = X2;
  }
};


////////////////////////////////////////////////////////////////////////////
// a function to convert from a floating point field of tensors
// to a field of bytes ... values will be remapped to range
// 0 ... 255.
// This is a base struct used for all Tenzors
template<class T, unsigned Dim>
struct ConvertTenzorField
{
  // calculate the minimum and maximum for a field of vektors
  template<class F1>
  static void minmax(F1 &a, T &rminval, T &rmaxval, bool firsttime) {
    // find min and max of a
    T minval = min(dotdot(a,a));
    T maxval = max(dotdot(a,a));

    // change returned values, if necessary
    if (firsttime || minval < rminval)
      rminval = minval;
    if (firsttime || maxval > rmaxval)
      rmaxval = maxval;
  }

  // the conversion function.  F1 should be a field of the
  // values, F2 a field of bytes with the same layout.
  template<class F1, class F2>
  static void convert(F1 &a, F2 &b, T minval, T maxval) {
    // convert to range 0 ... 255, and assign to bytes
    if (minval == maxval) {
      // everything is the same, so just make everything in b be zero
      b = 0;
    } else {
      // rescale the range
      minval = sqrt(minval);
      maxval = sqrt(maxval);
      b = (sqrt(dotdot(a,a)) - minval)/(maxval - minval) * 255.0;
    }
  }

  // accumulate the first field into the second
  template<class F1, class F2>
  static void accum(F1 &X2, F2 &X1, int, int) {
    X1 = X2;
  }
};

// Specializations of ConvertField that really just use the general
// tenzor case defined above by inheritance

template<class T, unsigned Dim>
struct ConvertField<Tenzor<T,Dim>, Dim>
  : public ConvertTenzorField<T,Dim>
{
};

template<class T, unsigned Dim>
struct ConvertField<SymTenzor<T,Dim>, Dim>
  : public ConvertTenzorField<T,Dim>
{
};

template<class T, unsigned Dim>
struct ConvertField<AntiSymTenzor<T,Dim>, Dim>
  : public ConvertTenzorField<T,Dim>
{
};


////////////////////////////////////////////////////////////////////////////
// a function to actually convert and write out the data for a field.
template<class T, unsigned Dim, class MT>
struct ConvData {
  static int conv(const char *fname1, const char *config1,
		  const char *outfname,
		  int fieldnum, int startrecordnum, int endrecordnum,
		  bool globalminmax, int minsize, bool dosum,
		  Inform &errmsg) {
    int i, j, k;
    Inform msg(fname1);

    // create DiscFields to read data, and check dimensions
    DiscField<Dim> io1(fname1, config1);
    int numRecords         = io1.get_NumRecords();
    int numFields          = io1.get_NumFields();
    NDIndex<Dim> domain    = io1.get_Domain();
    const char *typestring = io1.get_TypeString();

    // adjust start/end record nums if necessary
    if (startrecordnum < 0 || startrecordnum >= numRecords)
      startrecordnum = 0;
    if (endrecordnum < 0 || endrecordnum >= numRecords)
      endrecordnum = (numRecords - 1);

    // make sure everything is kosher
    if (numRecords < 1 && io1.get_Dimension() != Dim) {
      errmsg << "Bad dimension for file " << fname1 << endl;
      return RC_CMDERR;
    } else if (numRecords < 1 || startrecordnum > endrecordnum) {
      errmsg << "Not enough records in file " << fname1;
      errmsg << " or bad start/end record range." << endl;
      return RC_CMDERR;
    } else if (fieldnum >= numFields || fieldnum < 0 || numFields < 1) {
      errmsg << "Not enough fields in file " << fname1;
      errmsg << " or bad field number " << fieldnum << endl;
      return RC_CMDERR;
    }

    // determine the base of the output filename.  We'll append
    // an integer code to the end later.
    string outfnamebase = fname1;
    if (outfname != 0)
      outfnamebase = outfname;

    // figure out the domain of the output field.  It should have
    // sizes that are powers of two.  For now, we'll just pad with
    // zeros in the upper region of the data.
    NDIndex<Dim> outdomain;
    int maxlen = 0;
    int numbytes = 1;
    for (i=0; i < Dim; ++i) {
      int length = domain[i].length();
      int twopow = 1;
      while (twopow < length)
	twopow *= 2;
      if (twopow > maxlen)
	maxlen = twopow;
    }

    // if the user requested a minimum output size greater than the one
    // we want so far, use it, but only if its a power of two
    if (minsize > maxlen) {
      // make sure minsize is a power of two
      int testminsize = minsize;
      while (testminsize > 1) {
	if ((testminsize % 2) != 0) {
	  errmsg << "minsize (" << minsize << ") must be a power of two.";
	  errmsg << endl;
	  return RC_CMDERR;
	}
	testminsize /= 2;
      }

      // set our maxlen to the minimum size
      maxlen = minsize;
    }

    // finally, create an N-dimensional domain using the maximum length,
    // so that we get a domain with the same length in all dimensions
    for (i=0; i < Dim; ++i) {
      outdomain[i] = Index(domain[i].first(), domain[i].first() + maxlen - 1);
      numbytes *= maxlen;
    }

    // print out simple summary of contents
    msg << "Summary of data in file '" << fname1 << "':" << endl;
    msg << "        numRecords = " << numRecords << endl;
    msg << "  numFields/record = " << numFields << endl;
    msg << "      total domain = " << domain << endl << endl;

    // print out summary of what we're going to do
    msg << "Converting data to brick file '" << outfnamebase.c_str();
    msg << "':" << endl;
    msg << "      Field number = " << fieldnum << endl;
    msg << "      Start record = " << startrecordnum << endl;
    msg << "        End record = " << endrecordnum << endl;
    msg << "     Output domain = " << outdomain << endl;
    msg << "         Doing sum ? " << (dosum ? "Y" : "N") << endl;

    // initialize the input field ... create an initial layout from the domain
    FieldLayout<Dim>                           layout(domain);
    Field<T,Dim,UniformCartesian<Dim,MT>,Cell> X1(layout);
    Field<T,Dim,UniformCartesian<Dim,MT>,Cell> X2(layout);

    // now create a field of bytes to store the converted data, with the
    // same layout and size
    Field<unsigned char, Dim> bytes(layout);

    // create the output field ... this is always a field of bytes, with
    // just one vnode (that will be on node 0) and of the larger output
    // domain size
    FieldLayout<Dim>          outputLayout(outdomain, 0, 1);
    Field<unsigned char, Dim> outputBytes(outputLayout);
    outputBytes = 0;

    // if we are using a global min/max calculation, must run through
    // the data once and find global min/max
    MT globalmin, globalmax;
    if (globalminmax) {
      msg << "==========================================" << endl;
      msg << "Finding global minimum and maximum ..." << endl;
      msg << "==========================================" << endl;
      for (j=startrecordnum; j <= endrecordnum; ++j) {

        if (j > startrecordnum)
	  msg << "------------------------------------------" << endl;
	msg << "Finding min/max of record " << j;
        if (dosum)
          msg << " by summing all " << numFields << " fields";
        else
          msg << " for just field " << fieldnum;
	msg << " ..." << endl;

        // read in first field to process into X1
        k = (dosum ? 0 : fieldnum);
	msg << "  Reading field number " << k << " ..." << endl;
	if (!io1.read(X1, k, j))
	  return RC_READERR;

	// if we're summing, go through the input file and do the sum into
	// the parallel byte field for all the fields in this record
	if (dosum) {
          // We already read the first field, so read the rest into
          // a temporary and add to X1
	  for (k=1; k < numFields; ++k) {
	    // read just the selected record and field
	    msg << "  Reading field " << k << " ..." << endl;
	    if (!io1.read(X2, k, j))
	      return RC_READERR;

	    // accumulate X2 into X1 in proper places
	    msg << "  Adding in field " << k << " contribution ..." << endl;
            X1 += X2;
	    //ConvertField<T, Dim>::accum(X2, X1, k, numFields);
	  }
	}

	// calculate the minimum and maximum
	msg << "  Calculating min and max ..." << endl;
	ConvertField<T, Dim>::minmax(X1, globalmin, globalmax,
				     (j == startrecordnum));
	msg << "    .. new min=" << globalmin << ", new max=" << globalmax;
	msg << endl;
      }
    }

    // set up for connection
    bool connected = false;
    DataConnect *conn = 0;

    msg << endl;
    msg << "==========================================" << endl;
    msg << "Converting records from data file ..." << endl;
    msg << "==========================================" << endl;

    // now read in the data, convert it, and write it out
    for (j=startrecordnum; j <= endrecordnum; ++j) {

      if (j > startrecordnum)
        msg << "------------------------------------------" << endl;
      msg << "Starting conversion of record " << j;
      if (dosum)
        msg << " by summing all " << numFields << " fields";
      else
        msg << " for just field " << fieldnum;
      msg << " ..." << endl;

      // read in first field to process into X1
      k = (dosum ? 0 : fieldnum);
      msg << "  Reading field number " << k << " ..." << endl;
      if (!io1.read(X1, k, j))
        return RC_READERR;

      // if we're summing, go through the input file and do the sum into
      // the parallel byte field for all the fields in this record
      if (dosum) {
        // we've already read the first field, so just do the rest
	for (k=1; k < numFields; ++k) {
	  // read just the selected record and field
	  msg << "  Reading field " << k << " data ..." << endl;
	  if (!io1.read(X2, k, j))
	    return RC_READERR;

	  // accumulate X2 into X1 in proper places
	  msg << "  Adding in field " << k << " contribution ..." << endl;
          X1 += X2;
	  //ConvertField<T, Dim>::accum(X2, X1, k, numFields);
	}
      }

      // if we are not doing global min/max scaling, we must recompute
      // the min and max here
      if (!globalminmax) {
	msg << "  Recalculating min and max ..." << endl;
	ConvertField<T, Dim>::minmax(X1, globalmin, globalmax, true);
      }
      msg << "  Using min=" << globalmin << ", max=" << globalmax << endl;

      // convert to bytes ... how to do this depends on the data type
      msg << "  Converting data to bytes ..." << endl;
      ConvertField<T, Dim>::convert(X1, bytes, globalmin, globalmax);

      // remap the data to the output field, just a field of bytes
      // with no guard cells or anything.  Just copy over the smaller
      // domain's worth of data; the other values will be zero.
      msg << "  Copying data to node zero ..." << endl;
      outputBytes[domain] = bytes[domain];

      // remaining steps are only done on node 0
      if (Ippl::myNode() == 0) {
	// get pointer to start of bytes
	outputBytes.Uncompress();
	unsigned char *firstnum = (*(outputBytes.begin_if())).second->getP();

	// write out a new file ... first generate output filename
	char numbuf[32];
	sprintf(numbuf, ".%04d.dat", j - startrecordnum);
        string outfnamestr = outfnamebase;
	outfnamestr += numbuf;

	// then open a file
	FILE *f = fopen(outfnamestr.c_str(), "w");
	if (f == 0) {
	  msg << "Cannot open file '"<< outfnamestr << " for writing." << endl;
          return RC_WRITEERR;
	} else {
	  msg << "  Writing " << numbytes << " bytes to file '";
	  msg << outfnamestr << "' ..." << endl;
	  long written = fwrite(firstnum, sizeof(unsigned char), numbytes, f);
	  fclose(f);
          if (written != numbytes) {
            errmsg << "Error writing output file.  Only wrote " << written;
            errmsg << " bytes, trying to write " << numbytes << " bytes.";
            errmsg << endl;
            return RC_WRITEERR; 
	  }
	}
      }

      // register the Field for displaying data, if this is the first read
      if (!connected) {
	conn = DataConnectCreator::create(fname1);
	conn->connect("original data", X1, DataSource::OUTPUT);
	conn->connect("converted brick", outputBytes, DataSource::OUTPUT);
	connected = true;
      }
 
      if (conn != 0) {
	conn->updateConnections();
        if (j == startrecordnum)
	  conn->updateConnections();  // ACLVIS bug workaround
	conn->interact();
      }
    }

    // clean up
    if (conn != 0)
      delete conn;

    // return with exit status
    return RC_OK;
  }
};

#endif

/***************************************************************************
 * $RCSfile: dcf2brick_template.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:41:39 $
 * IPPL_VERSION_ID: $Id: dcf2brick_template.h,v 1.1.1.1 2003/01/23 07:41:39 adelmann Exp $ 
 ***************************************************************************/

