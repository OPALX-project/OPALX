#ifndef OPAL_ARefBas_HH
#define OPAL_ARefAttr_HH

// ------------------------------------------------------------------------
// $RCSfile: ARefAttr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class ARefAttr<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/04/07 12:02:51 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Expressions/ADeferred.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <vector>


namespace Expressions {

    // Class ARefAttr
    // ----------------------------------------------------------------------
    /// An attribute defined as a reference to an array.
    //  The components of the array may be real, logical or string.
    //  When a reference is seen, the pointers to the relevant object and
    //  attribute are left zero.  When the expression value is required,
    //  the object and the attribute are searched, and the pointers cached.
    //  The reference is registered with the object.  If the object referred
    //  to is deleted, it calls the invalidate() method of all reference
    //  expressions referring to it.  This resets the pointers to zero, so
    //  that the next evaluation forces to search for a replacement object.

    template <class T>
    class ARefAttr: public AttributeBase {

    public:

        /// Constructor.
        //  Use object name [b]oName[/b] to identify the object containing
        //  the array, and [b]aName[/b] to identify the array itself.
        ARefAttr(const string &oName, const string &aName);

        ARefAttr(const ARefAttr &);
        virtual ~ARefAttr();

        /// Make clone.
        virtual ARefAttr<T> *clone() const;

        /// Evaluate.
        //  Evaluate the reference and return the value.
        virtual std::vector<T> evaluate() const;

        /// Return real array value.
        //  This function has been added for speed of access.
        virtual vector<double> getRealArray();

        /// Invalidate.
        //  Force re-evaluation of the reference.
        virtual void invalidate();

        /// Print the reference.
        virtual void print(std::ostream &) const;

        /// Store new value.
        //  Evaluate the reference and assign to the array referred to.
        virtual void set(const std::vector<T> &) const;

    private:

        // Not implemented.
        ARefAttr();
        void operator=(const ARefAttr &);

        // Fill in the reference.
        void fill() const;

        // The name of the type referred to.
        static const string typeName;

        // The referred object and attribute.
        const string obj_name;
        const string att_name;

        // The object and attribute referred to.
        mutable Object    *itsObject;
        mutable Attribute *itsAttr;
    };


    template <class T>
    inline std::ostream &operator<<(std::ostream &os, const ARefAttr<T> &a) {
        a.print(os);
        return os;
    }


    // Implementation of class ARefAttr<T>.
    // ------------------------------------------------------------------------

    template <class T>
    ARefAttr<T>::ARefAttr(const string &oName, const string &aName):
        obj_name(oName), att_name(aName), itsObject(0), itsAttr(0)
    {}


    template <class T>
    ARefAttr<T>::ARefAttr(const ARefAttr &rhs):
        obj_name(rhs.obj_name), att_name(rhs.att_name),
        itsObject(rhs.itsObject), itsAttr(rhs.itsAttr)
    {}


    template <class T>
    ARefAttr<T>::~ARefAttr() {
        if(itsObject) itsObject->unregisterReference(this);
    }


    template <class T>
    ARefAttr<T> *ARefAttr<T>::clone() const {
        return new ARefAttr<T>(*this);
    }


    template <class T>
    std::vector<T> ARefAttr<T>::evaluate() const {
        fill();

        if(AttributeBase *base = &itsAttr->getBase()) {
            if(ADeferred<T> *value = dynamic_cast<ADeferred<T> *>(base)) {
                return value->evaluate();
            } else {
                throw OpalException("Real::get()", "Attribute \"" +
                                    itsAttr->getName() + "\" is of the wrong type.");
            }
        } else {
            return 0.0;
        }
    }


    template <class T>
    void ARefAttr<T>::invalidate() {
        itsObject = 0;
        itsAttr = 0;
    }


    template <class T>
    void ARefAttr<T>::print(std::ostream &os) const {
        os << obj_name;
        if(! att_name.empty()) os << "->" << att_name;
        return;
    }


    template <class T>
    void ARefAttr<T>::fill() const {
        if(itsObject == 0) {
            itsObject = OpalData::getInstance()->find(obj_name);
            if(itsObject == 0) {
                throw OpalException("ARefAttr::fill()",
                                    "Object \"" + obj_name + "\" is unknown.");
            }

            // Register the reference with the object, to allow invalidation
            // when the object is deleted.
            itsObject->registerReference(const_cast<ARefAttr<T>*>(this));

            if(att_name.empty()) {
                itsAttr = itsObject->findAttribute("VALUE");
                if(itsAttr == 0) {
                    throw OpalException("ARefAttr::fill()", "Object \"" + obj_name +
                                        "\" is not a variable, constant or vector.");
                }
            } else {
                itsAttr = itsObject->findAttribute(att_name);
                if(itsAttr == 0) {
                    throw OpalException("ARefAttr::fill()", "Attribute \"" + obj_name +
                                        "->" + att_name + "\" is unknown.");
                }
            }
        }
    }


    template <class T>
    vector<double> ARefAttr<T>::getRealArray() {
        throw OpalException("AValue<T>::getRealArray()",
                            "Attribute is not of real array type.");
    }


    template <> inline
    vector<double> ARefAttr<double>::getRealArray() {
        return evaluate();
    }


    template <class T>
    void ARefAttr<T>::set(const std::vector<T> &value) const {
        fill();

        if(AttributeBase *base = &itsAttr->getBase()) {
            if(dynamic_cast<ADeferred<T> *>(base)) {
                return itsAttr->set(new ADeferred<T>(value));
            } else {
                throw OpalException("Real::get()", "Attribute \"" +
                                    itsAttr->getName() + "\" is of the wrong type.");
            }
        }
    }

}

#endif // OPAL_ARefAttr_HH
