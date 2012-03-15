#ifndef SURFACEPHYSICSHANDLER_HH
#define SURFACEPHYSICSHANDLER_HH

#include <string>

class ElementBase;
class PartBunch;

class SurfacePhysicsHandler {
public:
    SurfacePhysicsHandler(std::string name, ElementBase *elref);
    virtual void apply(PartBunch &bunch) = 0;
    virtual const std::string getType() const = 0;

    void updateElement(ElementBase *newref);

protected:
    ElementBase *element_ref_m;

private:
    const std::string name_m;
};

inline SurfacePhysicsHandler::SurfacePhysicsHandler(std::string name, ElementBase *elref):
    name_m(name),
    element_ref_m(elref)
{}


inline void SurfacePhysicsHandler::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

#endif // SURFACEPHYSICS_HH
