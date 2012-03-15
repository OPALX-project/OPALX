#ifndef OPAL_AttCell_HH
#define OPAL_AttCell_HH 1

// ------------------------------------------------------------------------
// $RCSfile: AttCell.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttCell
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <map>
#include <string>

using std::string;


//  Class AttCell.
// ------------------------------------------------------------------------
/// The abstract base class for attribute cells.
//  Representation of a table cell for ATTLIST command.

class AttCell {

public:

    AttCell();
    virtual ~AttCell();

    /// Clear the value.
    //  Reset the value to undefined.
    virtual void clearValue() = 0;

    /// Print the attribute format.
    //  Print the format string in C style.
    virtual void printFormat(std::ostream &os) const = 0;

    /// Print the attribute value.
    //  According to the format string.
    virtual void printValue(std::ostream &os) const = 0;

    /// Store the value.
    //  Set the cell value to the given double.
    virtual void setReal(double);

    /// Store the value.
    //  Set the cell value to the given string.
    virtual void setString(const std::string &);

private:

    // Not implemented.
    AttCell(const AttCell &);
    void operator=(const AttCell &);
};


//  Class AttReal.
// ------------------------------------------------------------------------
/// The class for attribute cells with a real value.

class AttReal: public AttCell {

public:

    AttReal();
    virtual ~AttReal();

    /// Clear the value.
    virtual void clearValue();

    /// Print the attribute format.
    //  Prints "%le".
    virtual void printFormat(std::ostream &os) const;

    /// Print the attribute value.
    virtual void printValue(std::ostream &os) const;

    /// Store the value.
    virtual void setReal(double);

private:

    // Not implemented.
    AttReal(const AttReal &);
    void operator=(const AttReal &);

    // The attribute value.
    double itsValue;
};


//  Class AttString.
// ------------------------------------------------------------------------
/// The class for attribute cells with a string value.

class AttString: public AttCell {

public:

    AttString();
    virtual ~AttString();

    /// Clear the value.
    virtual void clearValue();

    /// Print the attribute format.
    //  Prints "%s".
    virtual void printFormat(std::ostream &os) const;

    /// Print the attribute value.
    virtual void printValue(std::ostream &os) const;

    /// Store the value.
    virtual void setString(const std::string &);

private:

    // Not implemented.
    AttString(const AttString &);
    void operator=(const AttString &);

    // The attribute value.
    string itsValue;
};

#endif // OPAL_AttCell_HH



