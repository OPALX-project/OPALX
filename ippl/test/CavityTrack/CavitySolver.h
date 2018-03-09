/***************************************************************************
                          cavitysolver.hh  -  description
                             -------------------
    begin                : Fri Juli 20 2004
    author               : Zenon Mathews
    email                : mathews.zenon@psi.ch
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                  This is the interface to the Cavity Solver,            *
 *                  which calculates the electric and magnetic fields      *
 *                  of particles in the cavity (if particle is outside     *
 *                  the cavity then the returned elec. and mag. fields     *
 *                  are zero)                                              *
 ***************************************************************************/

#ifndef CAVITY_SOLVER_H_
#define CAVITY_SOLVER_H_

// Trilinos includes
#include <Epetra_ConfigDefs.h>
#include <Teuchos_ParameterList.hpp>
// femaXX includes
#include "vector.h"
#include "myrlog.h"
#include "femaxxdriver.h"

class Epetra_Comm;

class CavitySolver 
{
public:
    /** Constructor, all the parameters need for femaxx (meshfile) and
        for the coordinate transformation(7 parameters) are stored
        into member variables. */
    CavitySolver(string meshfilename, 
                 double cavitylength,
                 double xtrans,
                 double ytrans,
                 double ztrans,
                 double xrot, 
                 double yrot, 
                 double zrot, 
                 double scaling, 
                 double peakvoltage, 
                 double phase);

    /** Destructor */
    ~CavitySolver();

    /** Runnning the femaxx solver, which computes the eigenvectors
        and the eigenvalues. */
    void runSolver();
 
    /** Compute the electric field of particle k inside the cavity k
        is the eigenvector with which the field value will be
        calculated!  (and not the particle nr) the rescaled, time
        dependant Efield is returned! */
    void getECavity(double* pos, double* ECavity, unsigned long k, double t);
    
    /** Compute the magnetic field of particle k inside the cavity the
        rescaled, time dependant Bfield is returned! */
    void getBCavity(double* pos, double* BCavity, unsigned long k, double t);

    /** This function transforms the FSTest coordinates into femaxx
        coordinates using the parameters read into member variables
        during the construction of this object. */
    void transform_coords(double* coords);
    
    /** 
     * Return the length of the cavity. 
     */
    double getlength() { return cavitylength_; }
    /**
     * Return k-th eigenvalue.
     */
    double getLambda(int k) { return eigenvalues_[k]; }
protected:
    /** Store k-th eigenpair in lambda, q on all processors.
      */
    void get_eigenpair_all(int k,     
                           Epetra_SerialDenseVector& L,
                           Epetra_MultiVector& Q,
                           double& lambda, 
                           colarray::Vector<double>& q);
private:
    const Epetra_Comm* comm_;
    FemaxxDriver* driver_;
    std::vector<double> eigenvalues_;
    std::vector< colarray::Vector<double> > eigenvectors_;

    string meshfile_;
    /** Gap length. Only used for scaling the fields. */
    double cavitylength_;
    /** Peak voltage. Only used for scaling the fields. The E-field
        is scaled by peakvoltage_/cavitylength_. This ensures that
        the gap voltage is equal to peakvoltage_, given that the
        (amplitude of) E-field is constant along the gap. */
    double peakvoltage_;
    /** Scaling factor (either 1.0 or -1.0) forcing the orientation
        of the E- and B-fields. */
    double factor_;
    /** Phase of time-harmonic the E- and B-field. */
    double phase_;
    double xtrans_;
    double ytrans_;
    double ztrans_;
    double xrot_;
    double yrot_;
    double zrot_;
    /** Scaling factor for coordinate transformation. */
    double scaling_;
    /** Rlog node: for output to terminal. */
    rlog::MyStdioNode* stdLog_;
    /** Rlog node: for output to log file. */
    rlog::StdioNode* dbgLog_;
};

#endif
