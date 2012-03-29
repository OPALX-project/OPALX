#ifndef CLASSIC_FIELDMAP1DPROFILE1_HH
#define CLASSIC_FIELDMAP1DPROFILE1_HH

#include "Fields/Fieldmap.hh"

// Class FM1DProfile1
//---------------------------------------------------------------------------
/// Field definition for 1D representation of bending magnet.
//
// Class FM1DProfile1 defines a 1D field map for use in bending magnets.

class FM1DProfile1: public Fieldmap {

public:

    /// Calculates the normalized (to 1) field strength on axis and the first and
    /// second derivatives (strength[0], stength[1], strength[2]) of the profile
    /// given the position, X.
    virtual bool getFieldstrength(const Vector_t &X, Vector_t &strength, Vector_t &info) const;

    virtual bool getFieldstrength_fdiff(const Vector_t &X, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;

    /// Get field dimensions.
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;

    virtual void swap();

    /// Get field information.
    virtual void getInfo(Inform *);

    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void setExitFaceSlope(const double &);

    /// Set exit edge constants.
    virtual void setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle);

    /// Set magnet gap.
    virtual void setFieldGap(const double &gap);

    /// Set distance between Enge function zeros.
    virtual void setFieldLength(const double &length);

    /// Adjust the extend of the fringe field regions. This is used
    /// to make sure we don't see funny jumps in the field as we
    /// move into and out of the magnet.
    virtual bool adjustFringeFields();

private:

    /// Constructor with field map file name.
    FM1DProfile1(std::string aFilename);

    virtual ~FM1DProfile1();

    /// Read field map.
    virtual void readMap();

    /// Free field map.
    virtual void freeMap();

    /// Enge coefficients for map entry and exit regions.
    double *EngeCoefs_entry_m;
    double *EngeCoefs_exit_m;

    /// Entry position of field map entry region.
    double zbegin_entry_m;

    /// End position of field map entry region.
    double zend_entry_m;

    /// Origin position for Enge function of field map entry region.
    double polynomialOrigin_entry_m;

    /// Enge function order for entry region.
    int polynomialOrder_entry_m;

    /// Entry position of field map exit region.
    double zbegin_exit_m;

    /// End position of field map exit region.
    double zend_exit_m;

    /// Origin position for Enge function of field map entry region.
    double polynomialOrigin_exit_m;

    /// Enge functioin order for entry region.
    int polynomialOrder_exit_m;

    /// Length of map (m).
    double length_m;

    /// Bend full gap height (m).
    double gapHeight_m;

    /// x position in local coordinate system where central trajectory intersects
    /// the exit edge.
    double xExit_m;

    /// z position in local coordinate system where central trajectory intersects
    /// the exit edge.
    double zExit_m;

    /// Cos and sin of the exit edge rotation with respect to the local coordinates.
    double cosExitRotation_m;
    double sinExitRotation_m;

    double exit_slope_m;

    friend class Fieldmap;
};

#endif
