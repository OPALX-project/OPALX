#ifndef PyElement_H

#include <Python.h>
#include <structmember.h>

#include <memory>
#include <exception>
#include <iostream>
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>
#include <boost/mpl/front.hpp>

#include "Utilities/OpalException.h"
#include "Elements/OpalElement.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbsBeamline/Component.h"
#include "Attributes/Attributes.h"
#include "OpalConfigure/Configure.h"

namespace PyOpal {

/** PyElementNS namespace contains PyElement, a wrapper for OpalElement objects,
 *  and various supporting objects.
 *
 *  PyElement<C>: basic element wrapper for C, which should be a subclass of
 *                OpalElement.
 *  AttributeType: enumeration of Opal Attribute Types (real, string, etc)
 *  AttributeDef: struct containing all of the things PyElement needs to know
 *                about each attribute that should be exposed to the python api.
 *  PyElementGetProperty: call policy to handle access of an Attribute for a 
 *                python property
 *  PyElementSetProperty: call policy to handle setting of an Attribute for a
 *                python property
 *
 *  To wrap an OpalElement, say MyOpalElementClass, you need to:
 *  1. Define static member data for the PyElement<MyOpalElementClass>. This 
 *     includes defining the attributes and setting a few options for different
 *     methods to expose.
 *  2. In BOOST_PYTHON_MODULE(my_opal_element_module), call 
 *
 *          PyElement<MyOpalElementType>.make_class()
 *
 *  Nb: apologies, this is heavy template stuff so almost everything has to go
 *  in the header file.
 */
namespace PyElementNS {

// forward declarations
template <class C> class PyElement;
template <class C> struct PyElementGetProperty;
template <class C> struct PyElementSetProperty;

/** AttributeType is used to control conversion from python to OpalAttribute 
 *  - Float will convert to RealAttribute
 *  - String will convert to StringAttribute
 *  - Bool will convert to BoolAttribute (tho in python Bool is alias to long)
 *  - Long will convert to RealAttribute
 */
enum AttributeType {DOUBLE, STRING, BOOL, INT};

/** Maps the AttributeType to a string representation for docstrings/etc */
extern std::map<AttributeType, std::string> attributeName; // defined in PyElement.cpp

/** AttributeDef defines an attribute
 *  opalName_m: the name of the opal Attribute
 *  pyName_m: the name that will be visible to the user in python. Properties
 *        should be_lower_case_with_underscores to comply with python API
 *  docString_m: docstring. If left empty (""), PyElement will generate a
 *        docstring like "py_name (type): Attribute help string"
 *  type_m: python type.
 */
struct AttributeDef {
    std::string opalName_m;
    std::string pyName_m;
    std::string docString_m;
    AttributeType type_m;
    // could also add a "read only" flag as some attributes are read only
};

/** PyElement<C> is the basic wrapper class for OpalElements
 *
 *  PyElement<C> is a wrapper for an OpalElement C. Opal Attributes are 
 *  implemented in python as properties (public member data). Some default
 *  method calls can be enabled in the concrete implementation by setting flags.
 *  Docstring can be enabled by setting classDocString.
 *
 *  element_m: pointer to the base element
 *  attributes: vector of attributes. The order of the vector is the order that
 *      attributes will appear in the docstring
 *  classDocString: overall docstring. docstrings on each property and method 
 *      are added automatically by the framework.
 *  hasGetFieldValue: set to true to add a get_field_value method to the class;
 *      this is useful for OpalElements that expose a field object.
 *
 *  property access is set up in make_class using add_property from
 *  boost::python. Normally, the "get" method can't take any data; this would
 *  mean that we can't pass an Attribute name for the OpalElement lookup. In
 *  order to get around this limitation, we have a "dummy" getter dummyGet.
 *  Then there is a "real" getter, which is hidden in the call policy for the
 *  property (a sort of pre- and post- decorator). We use PyElementGetProperty
 *  to hold the Attribute name and then dummyGet just tells PyElementGetProperty
 *  the address of the PyElement. PyElementGetProperty does a postcall action to
 *  overwrite the return value with the "real" return value, based on the stored
 *  Attribute name. 
 *  
 *  Setters work in exactly the same way, except now we have dummySet and
 *  PyElementSetProperty to handle the interface.
 *  
 *  The routine then looks like:
 *     python::boost calls dummyGet;
 *           dummyGet tells PyElementGetProperty pointer to PyElement
 *     python::boost calls PyElementGetProperty.postcall
 *           PyElementGetProperty.postcall calls PyElement->getAttribute; where
 *           getAttribute takes the name of the Attribute as a string.
 *  It's a bit of a faff.
 */
template <class C>
class PyElement {
public:
    typedef PyElement<C> PyC; // for convenience
    /** Default constructor */
    inline PyElement();
    /** Constructor taking the element as argument */
    PyElement(std::shared_ptr<C> element) : element_m(element) {}
    /** Copying is disabled */
    PyElement(const PyElement<C>& rhs);
    /** Default destructor */
    ~PyElement() {}

    /** This is the basic method to make a class. It should normally be called 
     *  when the module is declared, BOOST_PYTHON_MODULE.
     */
    boost::python::class_<PyC> make_class(const char* className);

    /** Add attributes to the python class.
     */
    void addAttributes(boost::python::class_<PyC>& pyclass);


    /** dummyGet sets the element ptr for PyElementGetProperty but doesn't
     *  actually do the get(...)
     */
    template <class ValueType>
    ValueType dummyGet() const {PyElementGetProperty<C>::setElement(this); return 0;}

    /** dummySet sets the element ptr for PyElementSetProperty but doesn't
     *  actually do the set(...)
     */
    template <class ValueType>
    void dummySet(ValueType test) {PyElementSetProperty<C>::setElement(this);}

    /** Get the value of an attribute
     *    - type: the type expected for python
     *    - opalName: the name of the Opal Attribute
     *  Returns a PyObject holding data; INCREF is called on the object i.e. it
     *  is a new object, not borrowed.
     *
     *  It is an error if the Attribute does not exist (throws OpalException)
     */
    PyObject* getAttribute(AttributeType type, std::string opalName) const;

    /** Set the value of an attribute
     *    - type: the type expected for python
     *    - opalName: the name of the Opal Attribute
     *    - value: Python object holding the value to be stored. value refcnt is
     *      not changed - caller is responsible for managing the refcnt.
     *  It is an error if the Attribute does not exist (throws OpalException)
     */
    void setAttribute(AttributeType type, std::string opalName, PyObject* value);

    /** Get the field at a point in position and time, for the OpalElement
     *    - x: position in local coordinate system of the OpalElement
     *    - y: position in local coordinate system of the OpalElement
     *    - z: position in local coordinate system of the OpalElement
     *    - t: time
     *  Returns a tuple like (out_of_bounds, Bx, By, Bz, Ex, Ey, Ez).
     *
     *  Only added to the python class if hasGetFieldValue is true.
     */
    boost::python::object getFieldValue(double x, double y, double z, double t) const;

    /** */
    boost::python::object getPyOpalElement();

    /** Returns the OpalElement from the PyElement */
    OpalElement* getOpalElement();

    C* getElement() {return element_m.get();}

protected:
    static std::vector<AttributeDef> attributes; /** class data (attributes) */
    static std::string classDocstring; /** class docstring */
    static bool converterRegistered; /** set to true if the converter has been registered */
    static bool hasGetFieldValue; /** field value */
    std::shared_ptr<C> element_m; /** pointer to the element */

    /** Generates a docstring from the attribute */
    std::string getDocString(AttributeDef& def);
};

/** Helper class to handle getting Attributes from python */
template <class C>
struct PyElementGetProperty : boost::python::default_call_policies {
public:
    /** Constructor
     *  - type: type of the python object
     *  - opalName: name of the Opal Attribute
     */
    PyElementGetProperty(AttributeType type, std::string opalName): type_m(type), opalName_m(opalName) {}
    /** destructor */
    ~PyElementGetProperty() {;}

    /** postcall action
     *  - ArgumentPackage: this is really a PyObject* holding the argument calls
     *  - result: filled by the postcall operation.
     *  Sets the element_m to NULL. Note that the whole element_m stuff is a bit
     *  hacky. Probably possible to extract the element from ArgumentPackage, I
     *  didn't manage to do it yet.
     */
    template <class ArgumentPackage>
    PyObject* postcall(ArgumentPackage const&, PyObject* result);

    /** Set pointer to the element; should be called before each postcall*/
    static void setElement(const PyElement<C>* element) {element_m = element;}

private:
    AttributeType type_m;
    std::string opalName_m;
    static const PyElement<C>* element_m;
};

/** Helper class to handle setting Attributes from python */
template <class C>
struct PyElementSetProperty : boost::python::default_call_policies {
public:
    /** Constructor
     *  - type: type of the python object
     *  - opalName: name of the Opal Attribute
     */
    PyElementSetProperty(AttributeType type, std::string opalName): type_m(type), opalName_m(opalName) {}
    /** destructor */
    ~PyElementSetProperty() {;}

    /** postcall action
     *  - ArgumentPackage: this is really a PyObject* holding the argument calls
     *  - result: filled by the postcall operation (with PyNone).
     *  Sets the element_m to NULL. Note that the whole element_m stuff is a bit
     *  hacky. Probably possible to extract the element from ArgumentPackage, I
     *  didn't manage to do it yet.
     */
    template <class ArgumentPackage>
    PyObject* postcall(ArgumentPackage const& args, PyObject* result);

    /** Set pointer to the element; should be called before each postcall*/
    static void setElement(PyElement<C>* element) {element_m = element;}

private:
    AttributeType type_m;
    std::string opalName_m;
    static PyElement<C>* element_m;
};

///////// templated implementations for PyElement ////////////////

// defined in PyOpalElement
template <>
PyElement<OpalElement>::PyElement();

template <class C>
PyElement<C>::PyElement() : element_m(new C) {}

template <class C>
bool PyElement<C>::converterRegistered = false;

template <class C>
boost::python::object PyElement<C>::getFieldValue(double x, double y, double z, double t) const {
    Vector_t R(x, y, z);
    Vector_t P(0.0, 0.0, 0.0);
    Vector_t B;
    Vector_t E;
    element_m->update();
    ElementBase* element = element_m->getElement()->removeWrappers();
    Component* component = dynamic_cast<Component*>(element);
    if (component == NULL) {
        throw OpalException("PyElement<C>::getFieldValue",
                            "Failed to deduce Component from ElementBase.");
    }
    bool outOfBounds = component->apply(R, P, t, E, B);
    return boost::python::make_tuple(outOfBounds,
                                     B[0], B[1], B[2],
                                     E[0], E[1], E[2]);
}

template <class C>
PyObject* PyElement<C>::getAttribute(AttributeType type, std::string opalName) const {
    if (!element_m) {
        throw OpalException("PyElement<C>::getRealAttribute",
                            "Element was not initialised");       
    }
    Attribute* attribute = element_m->findAttribute(opalName);
    if (attribute == NULL) {
        throw OpalException("PyElement<C>::getRealAttribute",
                            "Failed to parse attribute "+opalName);
    }
    PyObject* pyvalue;
    // I spent quite a bit of time trying to get this type unwinding to work
    // using templates (so it is done at compile time). In the end I couldn't
    // fight the template syntax and had to do it at runtime using an enum; it's
    // not so bad - the memory footprint is smaller and, tbh, template syntax is
    // horrible so might be easier to use.
    if (type == DOUBLE) {
        double value = Attributes::getReal(*attribute);
        pyvalue = PyFloat_FromDouble(value);
    } else if (type == INT) {
        double value = Attributes::getReal(*attribute);
        pyvalue = PyLong_FromDouble(value);
    }
    Py_INCREF(pyvalue);
    return pyvalue;
}


template <class C>
void PyElement<C>::setAttribute(AttributeType type, std::string opalName, PyObject* pyvalue) {
    if (!element_m) {
        throw OpalException("PyElement<C>::getRealAttribute",
                            "Element was not initialised");       
    }
    Attribute* attribute = element_m->findAttribute(opalName);
    if (attribute == NULL) {
        throw OpalException("PyElement<C>::getRealAttribute",
                            "Failed to parse attribute "+opalName);
    }
    if (type == DOUBLE) {
        double value = PyFloat_AsDouble(pyvalue);
        Attributes::setReal(*attribute, value);
    } else if (type == INT) {
        double value = PyLong_AsDouble(pyvalue);
        Attributes::setReal(*attribute, value);
    }
    return;
}

template <class C>
PyElement<C>::PyElement(const PyElement<C>& rhs) : element_m(rhs.element_m) { 
}

template <class C>
std::string PyElement<C>::getDocString(AttributeDef& def) {
    Attribute* attribute = element_m->findAttribute(def.opalName_m);
    if (attribute == NULL) {
        throw OpalException("PyElement<C>::getRealAttribute",
                            "Failed to parse attribute "+def.opalName_m);
    }
    std::string docString = def.pyName_m+" ("+attributeName[def.type_m]+"): "+attribute->getHelp();
    if (def.docString_m != "") {
        docString = def.pyName_m+" ("+attributeName[def.type_m]+"): "+def.docString_m;
    }
    return docString;
}

template <class C>
boost::python::class_<PyElement<C> > PyElement<C>::make_class(const char* className) {
    boost::python::docstring_options docop(true, true, false); // user_def, py_sig, cpp_sig
    auto pyclass = boost::python::class_<PyC>(className);
    addAttributes(pyclass);
    return pyclass;
}

template <class C>
void PyElement<C>::addAttributes(boost::python::class_<PyElement<C> >& pyclass) {
    for (std::vector<AttributeDef>::iterator iter = attributes.begin(); iter != attributes.end(); ++iter) {
        PyElementGetProperty<C> getProp(iter->type_m, iter->opalName_m);
        PyElementSetProperty<C> setProp(iter->type_m, iter->opalName_m);
        std::string docString = getDocString(*iter);
        std::string pyname = iter->pyName_m.c_str();
        if (iter->type_m == DOUBLE) {
            pyclass.add_property(pyname.c_str(),
                                 boost::python::make_function(&PyC::dummyGet<double>, getProp),
                                 boost::python::make_function(&PyC::dummySet<double>, setProp),
                                 docString.c_str()
            );
        } else if (iter->type_m == INT) {
            pyclass.add_property(pyname.c_str(),
                                 boost::python::make_function(&PyC::dummyGet<int>, getProp),
                                 boost::python::make_function(&PyC::dummySet<int>, setProp),
                                 docString.c_str()
            );
        }
    }
    pyclass.def("get_opal_element", &PyC::getPyOpalElement);
    if (hasGetFieldValue) {
        pyclass.def("get_field_value", &PyC::getFieldValue);
    }
}


template <class C>
OpalElement* PyElement<C>::getOpalElement() {
    OpalElement* element = dynamic_cast<OpalElement*>(element_m.get());
    if (element == NULL) {
        throw OpalException("PyElement<C>::getOpalElement",
                            "Failed to cast to OpalElement");
    }
    return element;
}

template <class C>
boost::python::object PyElement<C>::getPyOpalElement() {
    std::shared_ptr<OpalElement> elementPtr = std::dynamic_pointer_cast<OpalElement, C>(element_m);
    if (elementPtr.get() == NULL) {
        throw OpalException("PyElement<C>::getPyOpalElement", "Wrapped object was not an OpalElement");
    }
    PyElement<OpalElement> element(elementPtr);
    boost::python::object pyelement(element);
    return pyelement;
}

///////// templated implementations for PyElementGetProperty ////////////////

template <class C>
template <class ArgumentPackage>
PyObject* PyElementGetProperty<C>::postcall(ArgumentPackage const&, PyObject* result) {
        Py_DECREF(result);
        result = element_m->getAttribute(DOUBLE, opalName_m);
        return result;
}

template <class C>
const PyElement<C>* PyElementGetProperty<C>::element_m = NULL;

///////// templated implementations for PyElementSetProperty ////////////////

template <class C>
template <class ArgumentPackage>
PyObject* PyElementSetProperty<C>::postcall(ArgumentPackage const& args, PyObject* result) {
    PyObject* value;
    PyObject* pyElement; // this is a direct pointer to the C but I don't know how to unwrap it...
    if (!PyArg_ParseTuple(args, "OO", &pyElement, &value)) {
        return NULL; // ParseTuple sets the error message
    }
    Py_DECREF(result);
    element_m->setAttribute(DOUBLE, opalName_m, value);
    element_m = NULL;
    Py_RETURN_NONE;
}

template <class C>
PyElement<C>* PyElementSetProperty<C>::element_m = NULL;

} // PyElementNS
} // PyOpal

#endif // PyElement_H
