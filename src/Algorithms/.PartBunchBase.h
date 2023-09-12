//
// Class PartBunch
//   Base class for representing particle bunches.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef PART_BUNCH_BASE_H
#define PART_BUNCH_BASE_H

#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

// /#include "Particle/AbstractParticle.h"
#include "Particle/ParticleAttrib.h"
#include "Particle/ParticleLayout.h"

#include "AbstractObjects/OpalParticle.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/PBunchDefs.h"
#include "Algorithms/Quaternion.h"
 
#include "Distribution/Distribution.h"

#include "Algorithms/Matrix.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Units.h"
#include "Structure/FieldSolver.h"
#include "Utilities/GeneralClassicException.h"
#include "Utility/IpplTimings.h"

#include <memory>
#include <utility>
#include <vector>

class Distribution;
class FieldSolver;
class PartBins;
class PartBinsCyc;
class PartData;

template <class T, unsigned Dim>
class PartBunch : std::enable_shared_from_this<PartBunch<T, Dim>> {
public:
    // typedef typename AbstractParticle<T, Dim>::ParticlePos_t ParticlePos_t;
    // typedef typename AbstractParticle<T, Dim>::ParticleIndex_t ParticleIndex_t;
    // typedef typename AbstractParticle<T, Dim>::UpdateFlags UpdateFlags_t;
    // typedef typename AbstractParticle<T, Dim>::Position_t Position_t;

    typedef std::pair<Vector_t<double, 3>, Vector_t<double, 3>> VectorPair_t;

    static const unsigned Dimension = Dim;

    enum UnitState_t { units = 0, unitless = 1 };

public:
    virtual ~PartBunch() {
    }

    PartBunch(AbstractParticle<T, Dim>* pb, const PartData* ref);

    PartBunch(const PartBunch& rhs) = delete;  // implement if needed

    class ConstIterator {
        friend class PartBunch<T, Dim>;

    public:
        ConstIterator() : bunch_m(nullptr), index_m(0) {
        }
        ConstIterator(PartBunch const* bunch, unsigned int i) : bunch_m(bunch), index_m(i) {
        }

        ~ConstIterator() {
        }

        bool operator==(ConstIterator const& rhs) const {
            return bunch_m == rhs.bunch_m && index_m == rhs.index_m;
        }

        bool operator!=(ConstIterator const& rhs) const {
            return bunch_m != rhs.bunch_m || index_m != rhs.index_m;
        }

        OpalParticle operator*() const {
            if (index_m >= bunch_m->getLocalNum()) {
                throw GeneralClassicException(
                    "PartBunch::ConstIterator::operator*", "out of bounds");
            }
            return OpalParticle(
                bunch_m->ID[index_m], bunch_m->R[index_m], bunch_m->P[index_m], bunch_m->getT(),
                bunch_m->Q[index_m], bunch_m->getM() * Units::eV2MeV);
        }

        ConstIterator operator++() {
            ++index_m;
            return *this;
        }

        ConstIterator operator++(int) {
            ConstIterator it = *this;
            ++index_m;

            return it;
        }

        int operator-(const ConstIterator& other) const {
            return index_m - other.index_m;
        }

    private:
        PartBunch const* bunch_m;
        unsigned int index_m;
    };

    ConstIterator begin() const {
        return ConstIterator(this, 0);
    }

    ConstIterator end() const {
        return ConstIterator(this, getLocalNum());
    }

    /// Return maximum amplitudes.
    //  The matrix [b]D[/b] is used to normalise the first two modes.
    //  The maximum normalised amplitudes for these modes are stored
    //  in [b]axmax[/b] and [b]aymax[/b].

    void get_PBounds(Vector_t<double, 3>& min, Vector_t<double, 3>& max) const;

    Quaternion_t getQKs3D();
    void setQKs3D(Quaternion_t q);
    Vector_t<double, 3> getKs3DRefr();
    void setKs3DRefr(Vector_t<double, 3> r);
    Vector_t<double, 3> getKs3DRefp();
    void setKs3DRefp(Vector_t<double, 3> p);

    void iterateEmittedBin(int binNumber);

    void calcEMean();

    size_t getTotalNum() const;
    void setTotalNum(size_t n);
    void setLocalNum(size_t n);
    size_t getLocalNum() const;

    size_t getDestroyNum() const;
    size_t getGhostNum() const;

    bool getUpdateFlag(UpdateFlags_t f) const;
    void setUpdateFlag(UpdateFlags_t f, bool val);

    ParticleBConds<Position_t, Dimension>& getBConds() {
        return pbase_m->getBConds();
    }

    void setBConds(const ParticleBConds<Position_t, Dimension>& bc) {
        pbase_m->setBConds(bc);
    }

    void resetID();

    void update();
    void update(const ParticleAttrib<char>& canSwap);

    void createWithID(unsigned id);
    void create(size_t M);
    void globalCreate(size_t np);

    void destroy(size_t M, size_t I, bool doNow = false);
    void performDestroy(bool updateLocalNum = false);
    void ghostDestroy(size_t M, size_t I);

protected:
    size_t calcMoments();  // Calculates bunch moments using only emitted particles.

    /* Calculates bunch moments by summing over bins
     * (not accurate when any particles have been emitted).
     */
    void calcMomentsInitial();
    /// angle range [0~2PI) degree
    double calculateAngle(double x, double y);

private:
    virtual void updateDomainLength(Vektor<int, 3>& grid) = 0;

    virtual void updateFields(const Vector_t<double, 3>& hr, const Vector_t<double, 3>& origin);

    void setup(AbstractParticle<T, Dim>* pb);
};

template <class T, unsigned Dim>
typename PartBunch<T, Dim>::ConstIterator begin(PartBunch<T, Dim> const& bunch) {
    return bunch.begin();
}

template <class T, unsigned Dim>
typename PartBunch<T, Dim>::ConstIterator end(PartBunch<T, Dim> const& bunch) {
    return bunch.end();
}

#include "PartBunch.hpp"

#endif
