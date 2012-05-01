#ifndef CLASSIC_ParallelPlate_HH
#define CLASSIC_ParallelPlate_HH

// ------------------------------------------------------------------------
// $RCSfile: ParallelPlate.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelPlate
//   Defines the abstract interface for parallel plate element.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2010/10/12 $
// $Author: ChuanWang $
//
// ------------------------------------------------------------------------


#include "AbsBeamline/Component.h"

class PartBunch;
class Fieldmap;

// Class ParallelPlate
// ------------------------------------------------------------------------
/// Interface for RF cavity.
//  Class ParallelPlate defines the abstract interface for RF cavities.


class ParallelPlate: public Component {

public:

    //enum CavityType { SW, SGSW };
    /// Constructor with given name.
    explicit ParallelPlate(const string &name);

    ParallelPlate();
    ParallelPlate(const ParallelPlate &);
    virtual ~ParallelPlate();

    /// Apply visitor to ParallelPlate.
    virtual void accept(BeamlineVisitor &) const;



    void getDimensions(double &zBegin, double &zEnd) const;

    const string &getType() const;

    string getFieldMapFN() const;

    void setAmplitude(double vPeak);
    double getAmplitude() const ;

    void setFrequency(double freq);
    double getFrequency() const ;

    void setPhase(double phase);
    double getPhase() const ;

    void setElementLength(double length);
    double getElementLength() const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();





private:
    string filename_m;             /**< The name of the inputfile*/

    double scale_m;              /**< scale multiplier*/
    double phase_m;              /**< phase shift of time varying field(degrees)*/
    double frequency_m;          /**< Read in frequency of time varying field(MHz)*/
    double length_m;             /**< Read in distance/length of Parallel Plate*/
    double startField_m;
    double endField_m;
    double ElementEdge_m;
    double ptime_m;

    // Not implemented.
    void operator=(const ParallelPlate &);
};

#endif // CLASSIC_ParallelPlate_HH
