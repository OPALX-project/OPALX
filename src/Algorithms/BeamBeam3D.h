#ifndef CLASSIC_BeamBeam3D_HH
#define CLASSIC_BeamBeam3D_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamBeam3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam3D
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamBeam.h"
#include "Algorithms/Particle.h"
#include "BeamlineGeometry/NullGeometry.h"
#include "Fields/NullField.h"
#include "FixedAlgebra/FMatrix.h"
#include <complex>
#include <vector>

class PartBunch;

template <class T, int N> class FTps;
template <class T, int N> class FVps;

typedef FTps<double, 6> Series;
typedef FVps<double, 6> Map;


// Class BeamBeam3D
// ------------------------------------------------------------------------
/// A concrete representation for a beam-beam interaction.
// The ``strong'' bunch has a Gaussian distribution of the form
// [pre]
// N(x,y,s) = c * exp(- transpose(z - delta) sigma**(-1) (z - delta)),
// [/pre]
// with the definitions
// [ol]
//   [li] [b]z[/b]:     position vector in space,
//   [li] [b]delta[/b]: centroid of strong bunch in global reference system,
//   [li] [b]sigma[/b]: "beam" matrix in the TRANSPORT sense,
//   [li] [b]c[/b]:     a normalising factor, such that the total charge is Q,
//   [li] [b]sigma(i,i)[/b]: standard deviation sigma(i) in direction i,
//   [li] [b]r(i,j)[/b]: = sigma(i,j) / sqrt(sigma(i,i)*sigma(i,j)),
//     correlations between phase space coordinates $i$ and $j$.
// [/ol]
// ***** MISSING ***** This class needs further work, based on the report
// 6D Beam-Beam Kick including Coupled Motion by Leunissen et al.

class BeamBeam3D: public BeamBeam {

public:

    /// Collection of lattice functions and emittances.
    struct Beta {
        double betax, alphax, etax, etapx, emitx;
        double betay, alphay, etay, etapy, emity;
        double sigt, sige;
    };

    /// Constructor with given name.
    explicit BeamBeam3D(const std::string &name);

    BeamBeam3D();
    BeamBeam3D(const BeamBeam3D &right);
    virtual ~BeamBeam3D();

    /// Apply visitor to BeamBeam.
    virtual void accept(BeamlineVisitor &) const;

    /// Return exact copy of the element.
    virtual ElementBase *clone() const;

    /// Return a channel to an attribute.
    virtual Channel *getChannel(const std::string &aKey, bool create = false);

    /// Return the zero electromagnetic field.
    virtual NullField &getField();

    /// Return the zero electromagnetic field. Version for const object.
    virtual const NullField &getField() const;

    /// Return the null geometry.
    virtual NullGeometry &getGeometry();

    /// Return the null geometry. Version for const object.
    virtual const NullGeometry &getGeometry() const;

    /// Return an image of the element.
    virtual ElementImage *getImage() const;

    /// Return type name string.
    virtual ElementBase::ElementType getType() const;

    /// Get the bunch charge.
    //  Return number of particles times the particle charge in the strong
    //  bunch.  Units are proton charges.
    virtual double getBunchCharge() const;

    /// Get the moment matrix for the strong bunch.
    //  Return matrix of second momenta. Units are square metres.
    virtual const Matrix3D &getBunchMoment() const;

    /// Get displacement vector for the strong bunch.
    //  Return the displacement in metres.
    virtual const Vector3D &getBunchDisplacement() const;

    /// Set the pointer to the complex error function.
    //  Changing this pointer allows use of a different, potentially
    //  faster algorithm for determination of the kick.
    void setErrorFunctionPointer
    (std::complex<double> (*fun)(std::complex<double>));

    /// Set the crossing angle.
    void setCrossingAngle(double);

    /// Set the proportionality factor.
    void setBeamBeamParameter(double);

    /// Store the description of the strong beam.
    //  [b]disp[/b] is the displacement of the closed orbit of the strong beam.
    //  [b]latFun[/b] is the collection of lattice functions and emittances of
    //  the strong beam.
    void setBeamDescription(const Vector3D &disp, const Beta &latFun);

    /// Select number of slices for strong beam.
    void setSlices(int);

    /// Track a particle bunch.
    virtual void trackBunch
    (PartBunch &, const PartData &, bool revBeam, bool revTrack) const;

    /// Track a transfer map.
    virtual void trackMap
    (Map &map, const PartData &, bool revBeam, bool revTrack) const;

private:

    // Representation of a single slice in the strong beam.
    struct Slice {
        double xstar, ystar, zstar;      // Slice position.
        double sigx, sigpx, sigy, sigpy; // Standard deviations.
    };

    // Not implemented.
    void operator=(const BeamBeam3D &);

    // Lorenz boost.
    void boost(PartBunch &bunch) const;
    void boost(Map &map) const;

    // Inverse Lorenz boost.
    void boosti(PartBunch &bunch) const;
    void boosti(Map &map) const;

    // Compute internal quantities.
    void computeF();
    void computeSlices();

    // Apply synchro-beam collision.
    void synchroBeamCollision(PartBunch &bunch) const;
    void synchroBeamCollision(Map &map) const;

    // Apply beam-beam kick for one slice.
    void bbf(double sepx, double sepy, double sigxx, double sigyy,
             double &bbfx, double &bbfy, double &bbgx, double &bbgy)
    const;
    void bbf(const Series &sepx, const Series &sepy,
             const Series &sigxx, const Series &sigyy,
             Series &bbfx, Series &bbfy, Series &bbgx, Series &bbgy)
    const;

    // The zero magnetic field.
    NullField field;

    // The geometry (a straight geometry with zero length).
    NullGeometry geometry;

    // The pointer to the complex error function.
    // Can be replaced by a call to setErrorFunctionPointer().
    std::complex<double> (*errorFunction)(std::complex<double>);

    // The crossing angle.
    double phi;
    double cphi, sphi, tphi;

    // Factor describing the effect of one slice.
    double F;

    // The beam-beam parameter.
    double xiyn;

    // The extent of the strong beam.
    double sigx, sigy;

    // The normalising matrix.
    FMatrix<double, 6, 6> D;

    // Displacement of the strong beam.
    Vector3D displacement;

    // Lattice functions for the strong beam.
    Beta lf;

    // The table of slices.
    std::vector<Slice> slices;
    int nsli;

    // Maximum normalised amplitudes of this beam.
    mutable double axmax, aymax;
};

#endif // CLASSIC_BeamBeam3D_HH