#ifndef MAGNETICFIELD_H
#define MAGNETICFIELD_H


#include <iostream>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "Utilities/GeneralClassicException.h"
#include "AbsBeamline/Cyclotron.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Fields/NullField.h"

#include "Physics/Physics.h"


class MagneticField : public Cyclotron {
    
public:
    /*!
     * @param fmapfn specifies the fieldmap file. We support only
     * the CARBONCYCL format at the moment.
     */
    MagneticField(const std::string fmapfn, const double& symmetry);
    
    /*!
     * @param type of the magnetic field file (e.g. CARBONCYCL)
     * @param scaleFactor for the magnetic field.
     */
    void read(const std::string& type,const double& scaleFactor);
    
    /*!
     * Get the value of B_{z}, dB_{z}/dr and dB_{z}/dtheta
     * where theta is the azimuth angle and r the radius.
     * @param bint is the magnetic field B_{z} [kG]
     * @param brint is dB_{z}/dr [kG/m]
     * @param btint is dB_{z}/dtheta [kG/deg]
     * @param r is the radial interpolation point [m]
     * @param theta is the azimuthal interpolation point [deg]
     */
    void interpolate(double& bint,
                     double& brint,
                     double& btint,
                     const double& r,
                     const double& theta);
    
    
    /*! Do not call.
     * @returns a dummy value.
     */
    double getSlices() const;
    
    /*! Do not call.
     * @returns a dummy value.
     */
    double getStepsize() const;
    
    /*! Do not call.
     * @returns a null field.
     */
    const EMField &getField() const;
    
    /*! Do not call.
     * @returns a null field.
     */
    EMField &getField();
    
    /*! Do not call.
     * @returns a null pointer.
     */
    ElementBase* clone() const;
    
    /*! Do not call.
     * @returns a null geometry.
     */
    const BGeometryBase &getGeometry() const;
    
    /*! Do not call.
     * @returns a null geometry.
     */
    BGeometryBase &getGeometry();
    
private:
    NullGeometry nullGeom_m;
    NullField nullField_m;
};

#endif