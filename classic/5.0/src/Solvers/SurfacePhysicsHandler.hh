#ifndef SURFACEPHYSICSHANDLER_HH
#define SURFACEPHYSICSHANDLER_HH

#include <string>

class ElementBase;
class PartBunch;
class Inform;

class SurfacePhysicsHandler {
public:
    SurfacePhysicsHandler(std::string name, ElementBase *elref);
    virtual ~SurfacePhysicsHandler() { };
    virtual void apply(PartBunch &bunch) = 0;
    virtual const std::string getType() const = 0;
    virtual void print(Inform& os) = 0;
    void updateElement(ElementBase *newref);

protected:
    ElementBase *element_ref_m;

private:
    const std::string name_m;

};

inline SurfacePhysicsHandler::SurfacePhysicsHandler(std::string name, ElementBase *elref):
    element_ref_m(elref),
    name_m(name)
{}


inline void SurfacePhysicsHandler::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

#endif // SURFACEPHYSICS_HH
