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
    MagneticField();
    
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