/*
 *  Copyright (c) 2014, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CLASSIC_ABSBEAMLINE_VariableRFCavity_HH
#define CLASSIC_ABSBEAMLINE_VariableRFCavity_HH

#include "Algorithms/AbstractTimeDependence.h"
#include "Fields/EMField.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "AbsBeamline/Component.h"

class Fieldmap;

/** Class VariableRFCavity
 *
 *  Generates a field like
 *      E = E0(x, y, z)*a(t)*sin{f(t)*t-q(t)}
 *      B = B0(x, y, z)*a(t)*cos{f(t)*t-q(t)}
 *  where E0, B0 are user defined field maps, a(t), f(t), q(t) are time
 *  dependent amplitude, frequency, phase respectively; it is assumed that these
 *  quantities vary sufficiently slowly that Maxwell is satisfied.
 */
class VariableRFCavity: public Component {
  public:
    /// Constructor with given name.
    explicit VariableRFCavity(const std::string &name);
    VariableRFCavity(const VariableRFCavity &);
    VariableRFCavity();
    VariableRFCavity& operator=(const VariableRFCavity &);

    virtual ~VariableRFCavity();

    /// Apply visitor to RFCavity.
    virtual void accept(BeamlineVisitor &) const;

    virtual ElementBase* clone() const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);
    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);
    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);
    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);
    virtual void finalise();
    virtual bool bends() const {return false;}
    virtual void getDimensions(double &zBegin, double &zEnd) const {}

    virtual inline double getAmplitude(double time) const;
    virtual inline double getFrequency(double time) const;
    virtual inline double getPhase(double time) const;

    virtual double getHeight() const {return halfHeight_m*2;}
    virtual double getWidth() const {return halfWidth_m*2;}
    virtual void setHeight(double fullHeight) {halfHeight_m = fullHeight/2;}
    virtual void setWidth(double fullWidth) {halfWidth_m = fullWidth/2;}


    virtual std::shared_ptr<AbstractTimeDependence> getAmplitudeModel() const;
    virtual std::shared_ptr<AbstractTimeDependence> getPhaseModel() const;
    virtual std::shared_ptr<AbstractTimeDependence> getFrequencyModel() const;

    virtual void setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    virtual void setPhaseModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    virtual void setFrequencyModel(std::shared_ptr<AbstractTimeDependence> time_dep);

    virtual void setAmplitudeName(std::string amplitude)
    { amplitudeName_m = amplitude; }
    virtual void setPhaseName(std::string phase)
    { phaseName_m = phase; }
    virtual void setFrequencyName(std::string frequency)
    { frequencyName_m = frequency; }

    virtual void setLength(double length);
    virtual double getLength() const {return _length;}

    virtual StraightGeometry& getGeometry();
    virtual const StraightGeometry& getGeometry() const;

    /// Not implemented
    virtual EMField &getField();
    /// Not implemented
    virtual const EMField &getField() const;
  protected:
    std::shared_ptr<AbstractTimeDependence> _phase_td;
    std::shared_ptr<AbstractTimeDependence> _amplitude_td;
    std::shared_ptr<AbstractTimeDependence> _frequency_td;
    std::string phaseName_m;
    std::string amplitudeName_m;
    std::string frequencyName_m;
    double halfWidth_m;
    double halfHeight_m;
    double _length;
    /// The cavity's geometry.
    StraightGeometry geometry;

  private:
    void initNull();

};

double VariableRFCavity::getAmplitude(double time) const {
    return _amplitude_td->getValue(time);
}

double VariableRFCavity::getPhase(double time) const {
    return _phase_td->getValue(time);
}

double VariableRFCavity::getFrequency(double time) const {
    return _frequency_td->getValue(time);
}


#endif // CLASSIC_VirtualRFCavity_HH
