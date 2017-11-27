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
    virtual void apply(PartBunch &bunch, size_t numParticlesInSimulation = 0) = 0;
    virtual const std::string getType() const = 0;
    virtual void print(Inform& os) = 0;
    virtual bool stillActive() = 0;
    virtual bool stillAlive(PartBunch &bunch) = 0;
    virtual double getTime() = 0;
    virtual std::string getName() = 0;
    virtual size_t getParticlesInMat() = 0;
    virtual unsigned getRedifused() = 0;
    void AllParticlesIn(bool p);
    void updateElement(ElementBase *newref);
  
protected:
    ElementBase *element_ref_m;
    bool allParticleInMat_m;  
private:
    const std::string name_m;

};

inline SurfacePhysicsHandler::SurfacePhysicsHandler(std::string name, ElementBase *elref):
    element_ref_m(elref),
    allParticleInMat_m(false),
    name_m(name)
{}

inline void SurfacePhysicsHandler::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

inline void SurfacePhysicsHandler::AllParticlesIn(bool p) {
  allParticleInMat_m = p;
} 

#endif // SURFACEPHYSICS_HH
