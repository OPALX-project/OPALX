#ifndef PyOpalObject_H

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

/** PyOpalObjectNS namespace contains PyOpalObject, a wrapper for Object objects,
 *  and various supporting objects.
 *
 *  PyOpalObject<C>: basic element wrapper for C, which should be a subclass of
 *                OpalElement.
 *  AttributeType: enumeration of Opal Attribute Types (real, string, etc)
 *  AttributeDef: struct containing all of the things PyOpalObject needs to know
 *                about each attribute that should be exposed to the python api.
 *  PyElementGetProperty: call policy to handle access of an Attribute for a 
 *                python property
 *  PyElementSetProperty: call policy to handle setting of an Attribute for a
 *                python property
 *
 *  To wrap an OpalElement, say MyOpalObjectClass, you need to:
 *  1. Define static member data for the PyOpalObject<MyOpalObjectClass>. This 
 *     includes defining the attributes and setting a few options for different
 *     methods to expose.
 *  2. In BOOST_PYTHON_MODULE(my_opal_element_module), call 
 *
 *          PyOpalObject<MyOpalObjectType>.make_class()
 *
 *  Nb: apologies, this is heavy template stuff so almost everything has to go
 *  in the header file.
 */
namespace PyOpalObjectNS {

// forward declarations
template <class C> class PyOpalObject;
template <class C> struct PyOpalObjectGetProperty;
template <class C> struct PyOpalObjectSetProperty;

/** AttributeType is used to control conversion from python to OpalAttribute 
 *  - Float will convert to RealAttribute
 *  - String will convert to StringAttribute
 *  - Bool will convert to BoolAttribute (tho in python Bool is alias to long)
 *  - Long will convert to RealAttribute
 *  - VectorDouble will convert list to RealArray
 */
enum AttributeType {DOUBLE, STRING, BOOL, INT, FLOATLIST};

/** Maps the AttributeType to a string representation for docstrings/etc */
extern std::map<AttributeType, std::string> attributeName; // defined in PyOpalObject.cpp

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

/** PyOpalObject<C> is the basic wrapper class for Opal Objects
 *
 *  PyOpalObject<C> is a wrapper for an Opal Object C. Opal Attributes are 
 *  implemented in python as properties (public member data). Some default
 *  method calls can be enabled in the concrete implementation by setting flags.
 *  Docstring can be enabled by setting classDocString.
 *
 *  element_m: pointer to the base element
 *  attributes: vector of attributes. The order of the vector is the order that
 *      attributes will appear in the docstring
 *  classDocString: overall docstring. docstrings on each property and method 
 *      are added automatically by the framework.
  *
 *  property access is set up in make_class using add_property from
 *  boost::python. Normally, the "get" method can't take any data; this would
 *  mean that we can't pass an Attribute name for the OpalElement lookup. In
 *  order to get around this limitation, we have a "dummy" getter dummyGet.
 *  Then there is a "real" getter, which is hidden in the call policy for the
 *  property (a sort of pre- and post- decorator). We use PyOpalObjectGetProperty
 *  to hold the Attribute name and then dummyGet just tells PyElementGetProperty
 *  the address of the PyOpalObject. PyElementGetProperty does a postcall action to
 *  overwrite the return value with the "real" return value, based on the stored
 *  Attribute name. 
 *  
 *  Setters work in exactly the same way, except now we have dummySet and
 *  PyElementSetProperty to handle the interface.
 *  
 *  The routine then looks like:
 *     python::boost calls dummyGet;
 *           dummyGet tells PyElementGetProperty pointer to PyOpalObject
 *     python::boost calls PyElementGetProperty.postcall
 *           PyElementGetProperty.postcall calls PyOpalObject->getAttribute; where
 *           getAttribute takes the name of the Attribute as a string.
 *  It's a bit of a faff.
 */
template <class C>
class PyOpalObject {
public:
    typedef PyOpalObject<C> PyC; // for convenience
    /** Default constructor */
    inline PyOpalObject();
    /** Constructor taking the element as argument */
    PyOpalObject(std::shared_ptr<C> object) : object_m(object) {}
    /** Copying is disabled */
    PyOpalObject(const PyOpalObject<C>& rhs);
    /** Default destructor */
    ~PyOpalObject() {}

    /** This is the basic method to make a class. It should normally be called 
     *  when the module is declared, BOOST_PYTHON_MODULE.
     *
     *  Note that while in principle this could be static, in the end we need an
     *  instance of C so that we can access the docString from the attribute.
     */
    inline boost::python::class_<PyC> make_class(const char* className);

    /** Add attributes to the python class.
     */
    template <class PYCLASS>
    void addAttributes(PYCLASS& pyclass);

    /** Add an "execute" method to the python class (e.g. if it is an ACTION)
     */
    template <class PYCLASS>
    void addExecute(PYCLASS& pyclass);

    /** Add a "register" method to the python class (to register against opal global data header)
     */
    template <class PYCLASS>
    void addRegister(PYCLASS& pyclass);

    /** Add a "get_opal_element" method to the python class (to overload as an PyOpalElement)
     */
    template <class PYCLASS>
    void addGetOpalElement(PYCLASS& pyclass);

    /** dummyGet sets the object ptr for PyOpalObjectGetProperty but doesn't
     *  actually do the get(...)
     */
    template <class ValueType>
    ValueType dummyGet() const {PyOpalObjectGetProperty<C>::setObject(this); return ValueType();}

    /** dummySet sets the element ptr for PyOpalObjectSetProperty but doesn't
     *  actually do the set(...)
     */
    template <class ValueType>
    void dummySet(ValueType test) {PyOpalObjectSetProperty<C>::setObject(this);}

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

    /** Returns the Opal Object from the PyOpalObject */
    std::shared_ptr<C> getOpalShared() {return object_m;}

protected:
    static std::vector<AttributeDef> attributes; /** class data (attributes) */
    static std::string classDocstring; /** class docstring */
    static bool converterRegistered; /** set to true if the converter has been registered */
    std::shared_ptr<C> object_m; /** pointer to the element */

    /** Generates a docstring from the attribute */
    std::string getDocString(AttributeDef& def);
    static void execute(PyOpalObject<C>& pyobject);
    static void registerObject(PyOpalObject<C>& pyobject);
    static boost::python::object getPyOpalElement(PyOpalObject<C>& pyobject);
};

template <class C>
boost::python::object getFieldValue(PyOpalObjectNS::PyOpalObject<C>& pyobject, double x, double y, double z, double t) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->update();
    ElementBase* element = objectPtr->getElement()->removeWrappers();
    Component* component = dynamic_cast<Component*>(element);
    if (component == NULL) {
        throw OpalException("PyElement<C>::getFieldValue",
                            "Failed to deduce Component from ElementBase.");
    }

    Vector_t R(x, y, z);
    Vector_t P(0.0, 0.0, 0.0);
    Vector_t B;
    Vector_t E;
    bool outOfBounds = component->apply(R, P, t, E, B);
    return boost::python::make_tuple(outOfBounds,
                                     B[0], B[1], B[2],
                                     E[0], E[1], E[2]);
}





/** Helper class to handle getting Attributes from python */
template <class C>
struct PyOpalObjectGetProperty : boost::python::default_call_policies {
public:
    /** Constructor
     *  - type: type of the python object
     *  - opalName: name of the Opal Attribute
     */
    PyOpalObjectGetProperty(AttributeType type, std::string opalName): type_m(type), opalName_m(opalName) {}
    /** destructor */
    ~PyOpalObjectGetProperty() {;}

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
    static void setObject(const PyOpalObject<C>* object) {object_m = object;}

private:
    AttributeType type_m;
    std::string opalName_m;
    static const PyOpalObject<C>* object_m;
};

/** Helper class to handle setting Attributes from python */
template <class C>
struct PyOpalObjectSetProperty : boost::python::default_call_policies {
public:
    /** Constructor
     *  - type: type of the python object
     *  - opalName: name of the Opal Attribute
     */
    PyOpalObjectSetProperty(AttributeType type, std::string opalName): type_m(type), opalName_m(opalName) {}
    /** destructor */
    ~PyOpalObjectSetProperty() {;}

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
    static void setObject(PyOpalObject<C>* object) {object_m = object;}

private:
    AttributeType type_m;
    std::string opalName_m;
    static PyOpalObject<C>* object_m;
};

///////// templated implementations for PyElement ////////////////


template <class C>
void PyOpalObject<C>::execute(PyOpalObject<C>& pyobject) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->execute();
}

template <class C>
void PyOpalObject<C>::registerObject(PyOpalObjectNS::PyOpalObject<C>& pyobject) {
    C* wrappedC = pyobject.getOpalShared().get();
    Object* objectPtr = dynamic_cast<Object*>(wrappedC);
    if (objectPtr == NULL) {
        throw OpalException("PyOpalObject<C>::registerObject",
                            "Trying to register something that was not a Opal Object");
    }
    //Object* objectPtr = &(*pyobject.getOpalShared());
    OpalData::getInstance()->define(objectPtr);
}


template <class C>
boost::python::object PyOpalObject<C>::getPyOpalElement(PyOpalObjectNS::PyOpalObject<C>& pyobject) {
    std::shared_ptr<OpalElement> elementPtr =
            std::dynamic_pointer_cast<OpalElement, C>(pyobject.getOpalShared());
    if (elementPtr.get() == NULL) {
        throw OpalException("PyOpalObject<C>::getPyOpalElement",
                            "Wrapped object was not an OpalElement");
    }
    PyOpalObject<OpalElement> element(elementPtr);
    boost::python::object pyelement(element);
    return pyelement;
}


// defined in PyOpalElement
template <>
PyOpalObject<OpalElement>::PyOpalObject();

template <class C>
PyOpalObject<C>::PyOpalObject() : object_m(new C) {}

template <class C>
PyObject* PyOpalObject<C>::getAttribute(AttributeType type, std::string opalName) const {
    if (!object_m) {
        throw OpalException("PyOpalObject<C>::getRealAttribute",
                            "Object was not initialised");       
    }
    Attribute* attribute = object_m->findAttribute(opalName);
    if (attribute == NULL) {
        throw OpalException("PyOpalObject<C>::getRealAttribute",
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
    } else if (type == STRING) {
        std::string value = Attributes::getString(*attribute);
        pyvalue = PyUnicode_FromString(value.c_str());
    } else if (type == BOOL) {
        bool value = Attributes::getBool(*attribute);
        if (value) {
            pyvalue = Py_True;
        } else {
            pyvalue = Py_False;        
        }
    } else if (type == FLOATLIST) {
        std::vector<double> value = Attributes::getRealArray(*attribute);
        pyvalue = PyList_New(value.size());
        for (size_t i = 0; i < value.size(); ++i) {
            PyList_SetItem(pyvalue, i, PyFloat_FromDouble(value[i])); // WARNING check memory...
        }
    } else {
        throw OpalException("PyOpalObject<C>::getAttribute",
                            "Attribute type "+attributeName[type]+" not implemented");
    }
    Py_INCREF(pyvalue);
    return pyvalue;
}


template <class C>
void PyOpalObject<C>::setAttribute(AttributeType type, std::string opalName, PyObject* pyvalue) {
    if (!object_m) {
        throw OpalException("PyOpalObject<C>::setAttribute",
                            "Element was not initialised");       
    }
    Attribute* attribute = object_m->findAttribute(opalName);
    if (attribute == NULL) {
        throw OpalException("PyOpalObject<C>::setAttribute",
                            "Failed to parse attribute "+opalName);
    }
    if (type == DOUBLE) {
        double value = PyFloat_AsDouble(pyvalue);
        Attributes::setReal(*attribute, value);
    } else if (type == INT) {
        double value = PyLong_AsDouble(pyvalue);
        Attributes::setReal(*attribute, value);
    } else if (type == STRING) {
        std::string value = PyUnicode_AsUTF8(pyvalue);
        Attributes::setString(*attribute, value);
    } else if (type == BOOL) {
        bool value = PyObject_IsTrue(pyvalue);
        Attributes::setBool(*attribute, value);
    } else if (type == FLOATLIST) {
        Py_ssize_t listSize = PyList_Size(pyvalue);
        std::vector<double> value(listSize);
        for (Py_ssize_t i = 0; i < listSize; ++i) {
            double value_i = PyFloat_AsDouble(PyList_GetItem(pyvalue, i));
            value[i] = value_i;
        }
        Attributes::setRealArray(*attribute, value);
    } else {
        throw OpalException("PyOpalObject<C>::setAttribute",
                            "Attribute type "+attributeName[type]+" not implemented");
    }
}

template <class C>
PyOpalObject<C>::PyOpalObject(const PyOpalObject<C>& rhs) : object_m(rhs.object_m) { 
}

template <class C>
std::string PyOpalObject<C>::getDocString(AttributeDef& def) {
    Attribute* attribute = object_m->findAttribute(def.opalName_m);
    if (attribute == NULL) {
        throw OpalException("PyOpalObject<C>::getRealAttribute",
                            "Failed to parse attribute "+def.opalName_m);
    }
    std::string docString = def.pyName_m+" ("+attributeName[def.type_m]+"): "+attribute->getHelp();
    if (def.docString_m != "") {
        docString = def.pyName_m+" ("+attributeName[def.type_m]+"): "+def.docString_m;
    }
    return docString;
}

template <class C>
boost::python::class_<PyOpalObject<C> > PyOpalObject<C>::make_class(const char* className) {
    typedef boost::python::class_<PyOpalObject<C> > PyClass;
    boost::python::docstring_options docop(true, true, false); // user_def, py_sig, cpp_sig
    PyClass pyclass = PyClass(className);
    addAttributes(pyclass);
    return pyclass;
}


template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addExecute(PYCLASS& pyclass) {
    pyclass.def("execute", &PyOpalObject<C>::execute);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addRegister(PYCLASS& pyclass) {
    pyclass.def("register", &PyOpalObject<C>::registerObject);
}

template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addGetOpalElement(PYCLASS& pyclass) {
    pyclass.def("get_opal_element", &PyOpalObject<C>::getPyOpalElement);
}



template <class C>
template <class PYCLASS>
void PyOpalObject<C>::addAttributes(PYCLASS& pyclass) {
    for (std::vector<AttributeDef>::iterator iter = attributes.begin(); iter != attributes.end(); ++iter) {
        PyOpalObjectGetProperty<C> getProp(iter->type_m, iter->opalName_m);
        PyOpalObjectSetProperty<C> setProp(iter->type_m, iter->opalName_m);
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
         }  else if (iter->type_m == STRING) {
            pyclass.add_property(pyname.c_str(),
                                 boost::python::make_function(&PyC::dummyGet<std::string>, getProp),
                                 boost::python::make_function(&PyC::dummySet<std::string>, setProp),
                                 docString.c_str()
            );
        }  else if (iter->type_m == BOOL) {
            pyclass.add_property(pyname.c_str(),
                                 boost::python::make_function(&PyC::dummyGet<bool>, getProp),
                                 boost::python::make_function(&PyC::dummySet<bool>, setProp),
                                 docString.c_str()
            );
        }  else if (iter->type_m == FLOATLIST) {
            pyclass.add_property(pyname.c_str(),
                                 boost::python::make_function(&PyC::dummyGet<boost::python::list>, getProp),
                                 boost::python::make_function(&PyC::dummySet<boost::python::list>, setProp),
                                 docString.c_str()
            );
        } else {
            // Looks like exception handling doesn't work properly at module
            // import time so this may not be handled politely - thrown as an 
            // unrecognised SystemError
            throw OpalException("PyOpalObject<C>::addAttributes", "Type not implemented");
        }
    }
}


///////// templated implementations for PyOpalObjectGetProperty ////////////////

template <class C>
template <class ArgumentPackage>
PyObject* PyOpalObjectGetProperty<C>::postcall(ArgumentPackage const&, PyObject* result) {
        Py_DECREF(result);
        result = object_m->getAttribute(type_m, opalName_m);
        return result;
}

template <class C>
const PyOpalObject<C>* PyOpalObjectGetProperty<C>::object_m = NULL;

///////// templated implementations for PyOpalObjectSetProperty ////////////////

template <class C>
template <class ArgumentPackage>
PyObject* PyOpalObjectSetProperty<C>::postcall(ArgumentPackage const& args, PyObject* result) {
    PyObject* value;
    PyObject* pyObject; // this is a direct pointer to the C but I don't know how to unwrap it...
    if (!PyArg_ParseTuple(args, "OO", &pyObject, &value)) {
        return NULL; // ParseTuple sets the error message
    }
    Py_DECREF(result);
    object_m->setAttribute(type_m, opalName_m, value);
    object_m = NULL;
    Py_RETURN_NONE;
}

template <class C>
PyOpalObject<C>* PyOpalObjectSetProperty<C>::object_m = NULL;

} // PyOpalObject
} // PyOpal

#endif // PyOpalObject_H
