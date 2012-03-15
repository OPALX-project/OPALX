#ifndef MAD_IdealMapper_HH
#define MAD_IdealMapper_HH

// ------------------------------------------------------------------------
// $RCSfile: IdealMapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: IdealMapper
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/LinearMapper.h"
#include "FixedAlgebra/LinearMap.h"


class BMultipoleField;
class Euclid3D;

template <class T, int M, int N> class FMatrix;
template <class T, int N> class FTps;
template <class T, int N> class FVps;


// Class IdealMapper
// ------------------------------------------------------------------------
/// Build a map using the linear map around the design orbit for each element.
// Ignores the following:
// [ul]
// [li]Any misalignments,
// [li]Any field errors,
// [li]Any momentum errors.
// [/ul]
// Phase space coordinates numbering:
// [tab 3 b]
// [row]number [&]name        [&]unit  [/row]
// [row]0      [&]x           [&]metres [/row]
// [row]1      [&]p_x/p_r     [&]1      [/row]
// [row]2      [&]y           [&]metres [/row]
// [row]3      [&]p_y/p_r     [&]1      [/row]
// [row]4      [&]v*delta_t   [&]metres [/row]
// [row]5      [&]delta_p/p_r [&]1      [/row]
// [/tab][p]
// Where $p_r$ is the constant reference momentum defining the reference
// frame velocity, $m$ is the rest mass of the particles, and $v$ is the
// instantaneous velocity of the particle.
// [p]
// Other units used:
// [tab 2 b]
// [row]quantity             [&]unit           [/row]
// [row]reference momentum   [&]electron-volts [/row]
// [row]velocity             [&]metres/second  [/row]
// [row]accelerating voltage [&]volts          [/row]
// [row]separator voltage    [&]volts          [/row]
// [row]frequencies          [&]hertz          [/row]
// [row]phase lags           [&]$2*pi$         [/row]
// [/tab][p]
// Approximations used:
// [ul]
// [li] All elements are represented by maps for finite-length elements.
// [li] Geometric transformations ignore rotations about transverse axes and
//   translations along the design orbit and truncate after first order.
// [li] Beam-beam elements are ignored.
// [/ul]

class IdealMapper: public LinearMapper {

public:

  /// Constructor.
  //  The beam line to be tracked is [b]bl[/b].
  //  The particle reference data are taken from [b]data[/b].
  //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
  //  If [b]revTrack[/b] is true, we track against the beam.
  IdealMapper(const Beamline &beamline, const PartData &,
	      bool revBeam, bool revTrack);

  virtual ~IdealMapper();

  /// Return the linear part of the accumulated map.
  virtual void getMatrix(FMatrix<double,6,6> &) const;

  /// Reset the linear part of the accumulated map for restart.
  virtual void setMatrix(const FMatrix<double,6,6> &);

  /// Apply the algorithm to a Corrector.
  //  Override to ignore corrector kick.
  virtual void visitCorrector(const Corrector &);

  /// Apply the algorithm to a patch.
  //  Override to ignore patch.
  virtual void visitPatch(const Patch &pat);

  /// Apply the algorithm to a Separator.
  //  Override to ignore separator kick.
  virtual void visitSeparator(const Separator &);
  
  /// Apply the algorithm to an offset beamline object wrapper.
  //  Override to ignore misalignment.
  virtual void visitAlignWrapper(const AlignWrapper &);
  
  /// Apply the algorithm to an corrector wrapper..
  virtual void visitCorrectorWrapper(const CorrectorWrapper &);
  
  /// Apply the algorithm to an multipole wrapper..
  virtual void visitMultipoleWrapper(const MultipoleWrapper &);
  
  /// Apply the algorithm to an RBend wrapper..
  virtual void visitRBendWrapper(const RBendWrapper &);
  
  /// Apply the algorithm to an SBend wrapper..
  virtual void visitSBendWrapper(const SBendWrapper &);

protected:

  // Apply thin multipole kick (integrated over length) to ideal map.
  virtual void applyMultipoleBody
  (double length, double refLength, const BMultipoleField &field, double scale);

  // Apply thin multipole kick (integrated over length) to all particles.
  virtual void applySBendBody
  (double length, double refLength, double h,
   const BMultipoleField &field, double scale);

  /// Thin multipole kick.
  //  Apply a thin multipole kick (integrated over length) to current map.
  virtual void applyThinMultipole
  (const BMultipoleField &field, double factor);

  /// Thin SBend kick.
  //  Special kick routine for thin SBend.
  virtual void applyThinSBend
  (const BMultipoleField &field, double scale, double h);

  /// Apply transform.
  //  Propagate current map through a geometric transformation.
  virtual void applyTransform(const Euclid3D &, double refLength);

private:

  // Not implemented.
  IdealMapper();
  IdealMapper(const IdealMapper &);
  void operator=(const IdealMapper &);

  /// Apply linear map, defined by the linear coefficients.
  void applyLinearMap(double length, double kx, double ks, double ky);

  /// Helper function for finding first-order coefficients.
  static void makeFocus(double k, double L, double &c, double &s);
};

#endif // MAD_IdealMapper_HH
