#ifndef PARTICLEMATERINTERACTIONHANDLER_HH
#define PARTICLEMATERINTERACTIONHANDLER_HH

#include <string>
#include "Algorithms/Vektor.h"

class ElementBase;
class PartBunch;
class Inform;

class ParticleMaterInteractionHandler {
public:
    ParticleMaterInteractionHandler(std::string name, ElementBase *elref);
    virtual ~ParticleMaterInteractionHandler() { };
    virtual void apply(PartBunch &bunch,
                       const std::pair<Vector_t, double> &boundingSphere,
                       size_t numParticlesInSimulation = 0) = 0;
    virtual const std::string getType() const = 0;
    virtual void print(Inform& os) = 0;
    virtual bool stillActive() = 0;
    virtual bool stillAlive(PartBunch &bunch) = 0;
    virtual double getTime() = 0;
    virtual std::string getName() = 0;
    virtual size_t getParticlesInMat() = 0;
    virtual unsigned getRediffused() = 0;
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
ParticleMaterInteractionHandler::ParticleMaterInteractionHandler(std::string name, ElementBase *elref):
    element_ref_m(elref),
    allParticleInMat_m(false),
    name_m(name)
{}

inline
void ParticleMaterInteractionHandler::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

inline
ElementBase* ParticleMaterInteractionHandler::getElement() {
    return element_ref_m;
}

inline
void ParticleMaterInteractionHandler::setFlagAllParticlesIn(bool p) {
  allParticleInMat_m = p;
}

inline
bool ParticleMaterInteractionHandler::getFlagAllParticlesIn() const {
    return allParticleInMat_m;
}
#endif // PARTICLEMATERINTERACTION_HH