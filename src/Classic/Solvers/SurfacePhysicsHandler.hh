#ifndef SURFACEPHYSICSHANDLER_HH
#define SURFACEPHYSICSHANDLER_HH

#include <string>
#include "Algorithms/Vektor.h"

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;
class Inform;

class SurfacePhysicsHandler {
public:
    SurfacePhysicsHandler(std::string name, ElementBase *elref);
    virtual ~SurfacePhysicsHandler() { };
    virtual void apply(PartBunchBase<double, 3> &bunch,
                       const std::pair<Vector_t, double> &boundingSphere,
                       size_t numParticlesInSimulation = 0) = 0;
    virtual const std::string getType() const = 0;
    virtual void print(Inform& os) = 0;
    virtual bool stillActive() = 0;
    virtual bool stillAlive(PartBunchBase<double, 3> &bunch) = 0;
    virtual double getTime() = 0;
    virtual std::string getName() = 0;
    virtual size_t getParticlesInMat() = 0;
    virtual unsigned getRedifused() = 0;
    void setFlagAllParticlesIn(bool p);
    bool getFlagAllParticlesIn() const;
    void updateElement(ElementBase *newref);
    ElementBase* getElement();

protected:
    ElementBase *element_ref_m;
    bool allParticleInMat_m;
private:
    const std::string name_m;

};

inline
SurfacePhysicsHandler::SurfacePhysicsHandler(std::string name, ElementBase *elref):
    element_ref_m(elref),
    name_m(name),
    allParticleInMat_m(false)
{}

inline
void SurfacePhysicsHandler::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

inline
ElementBase* SurfacePhysicsHandler::getElement() {
    return element_ref_m;
}

inline
void SurfacePhysicsHandler::setFlagAllParticlesIn(bool p) {
  allParticleInMat_m = p;
}

inline
bool SurfacePhysicsHandler::getFlagAllParticlesIn() const {
    return allParticleInMat_m;
}
#endif // SURFACEPHYSICS_HH