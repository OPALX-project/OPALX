#ifndef OPAL_SURFACEPHYSICS_HH
#define OPAL_SURFACEPHYSICS_HH

// ------------------------------------------------------------------------
// $RCSfile: Wake.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Wake
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Solvers/SurfacePhysicsHandler.hh"
class ElementBase;
class Inform;

// Class Wake
// ------------------------------------------------------------------------
/// The WAKE definition.
//  A WAKE definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class SurfacePhysics: public Definition {

public:

    /// Exemplar constructor.
    SurfacePhysics();

    virtual ~SurfacePhysics();

    /// Test if replacement is allowed.
    //  Can replace only by another WAKE.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual SurfacePhysics *clone(const std::string &name);

    /// Check the SURFACEPHYSICS data.
    virtual void execute();

    /// Find named SURFACEPHYSICS.
    static SurfacePhysics *find(const std::string &name);

    /// Update the SURFACEPHYSICS data.
    virtual void update();

    /// Print the TFS descriptors for the surfac physics
    void tfsDescriptors(std::ostream &os) const;

    void print(std::ostream &os) const;

    void initSurfacePhysicsHandler(ElementBase &element, const double &major, const double &minor);

    void updateElement(ElementBase *element);

    SurfacePhysicsHandler *handler_m;

private:

    // Not implemented.
    SurfacePhysics(const SurfacePhysics &);
    void operator=(const SurfacePhysics &);

    // Clone constructor.
    SurfacePhysics(const std::string &name, SurfacePhysics *parent);

    // The particle reference data.
    PartData reference;

    // The conversion from GeV to eV.
    static const double energy_scale;

    // the element the surface physics is attached to
    ElementBase *itsElement_m;
    std::string material_m;

};

inline std::ostream &operator<<(std::ostream &os, const SurfacePhysics &b) {
    b.print(os);
    return os;
}

#endif // OPAL_SURFACEPHYSICS_HH
