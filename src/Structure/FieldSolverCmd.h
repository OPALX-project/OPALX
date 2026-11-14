//
// Class FieldSolverCmd
//   The class for the OPAL FIELDSOLVER command.
//   Stores parsed mesh, boundary, and algorithm attributes. SpaceCharge owns
//   their conversion and solver compatibility validation.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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
#ifndef OPAL_FieldSolver_HH
#define OPAL_FieldSolver_HH

#include <string>
#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "PartBunch/BCHandler.hpp"
#include "Structure/BinningCmd.h"
#include "Track/ObjectValues.h"

#include "Ippl.h"

enum class FieldSolverCmdType : short { NONE = -1, FFT = 0, OPEN = 1, CG = 2, P3M = 3, FFT2D5 = 4 };

// The attributes of class FieldSolverCmd.
namespace FIELDSOLVER {
    enum {
        TYPE,          // The field solver name
        BINS,          // Name of BINNING definition or NONE
        NX,            // mesh size in x
        NY,            // mesh size in y
        NZ,            // mesh size in z
        PARFFTX,       // parallelized grid in x
        PARFFTY,       // parallelized grid in y
        PARFFTZ,       // parallelized grid in z
        BCFFTX,        // Poisson-domain boundary in x
        BCFFTY,        // Poisson-domain boundary in y
        BCFFTZ,        // Poisson-domain boundary in z
        GREENSF,       // Green function for OPEN; P3M selects its kernel internally
        P3MRCUT,       // P3M particle-particle cutoff radius [m]
        P3MRCUTCELLS,  // P3M cutoff in units of the largest mesh spacing
        PERIODICBOX,   // edge length of a fixed all-periodic cube [m]; zero follows the particles
        BBOXINCR,      // how much the boundingbox is increased
        FFT2D5MODEL,   // Explicit FFT2D5 model identity; ORIGINAL is currently supported
        PIPEMODE,      // One of OPEN, CIRCULAR, PLATES, NONE [FFT2D5 only]
        BEAMR,         // Beam radius in metres [FFT2D5 only]
        CLOSEDRING,    // TRUE if the ring is closed [FFT2D5 only]
        PIPESIZEX,     // Size of the pipe in meters in the transverse direction [FFT2D5 only]
        PIPESIZEY,     // Size of the pipe in meters in the vertical direction [FFT2D5 only]
        REFPATHFNAME,  // Reference path file name [FFT2D5 only]
        SCATTERLONGITUDINALLY,  // Scatter charge between longitudinal slices [FFT2D5 only]
        SIZE
    };
}

/** @brief Parsed FIELDSOLVER attributes; solver compatibility is validated in SpaceCharge. */
class FieldSolverCmd : public Definition {
public:
    /// Exemplar constructor.
    FieldSolverCmd();

    virtual ~FieldSolverCmd();

    /// Find named FieldSolverCmd.
    static FieldSolverCmd* find(const std::string& name);

    static std::string getType(const ObjectValues& values);
    static std::string getBinsName(const ObjectValues& values);
    BinningCmd* getBinningCmd() const;
    static std::string getGreensFunction(const ObjectValues& values);
    static double getP3MCutoff(const ObjectValues& values);
    /// P3M cutoff in units of the largest mesh spacing (default 2 when RCUT is not given).
    static double getP3MCutoffCells(const ObjectValues& values);
    static double getPeriodicBox(const ObjectValues& values);
    /// Whether RCUT or RCUT_CELLS was set in the input, rather than left at its default.
    static bool hasP3MCutoff(const ObjectValues& values);
    static bool hasP3MCutoffCells(const ObjectValues& values);

    /// Returns solver boundary conditions handler object.
    static BCHandler<3> constructBCHandler(const ObjectValues& values);

    /// Return meshsize
    static double getNX(const ObjectValues& values);

    /// Return meshsize
    static double getNY(const ObjectValues& values);

    /// Return meshsize
    static double getNZ(const ObjectValues& values);

    void setNX(double);

    void setNY(double);

    void setNZ(double);

    static double getBoxIncr(const ObjectValues& values);

    /// Execute (init) the field solver data.
    virtual void execute();

    static FieldSolverCmdType getFieldSolverCmdType(const ObjectValues& values);

    static ippl::Vector<bool, 3> getDomainDecomposition(const ObjectValues& values);

    static Inform& printInfo(Inform& os, const ObjectValues& values);

    /// Return the configured FFT2D5 identity; ORIGINAL is the default and only supported model.
    static std::string getFFT2D5Model(const ObjectValues& values);
    /// Whether FFT2D5MODEL was explicitly supplied rather than inherited from its default.
    static bool hasFFT2D5Model(const ObjectValues& values);
    /** @brief Set a predefined FFT2D5 model identity.
     * @param model ORIGINAL; invalid predefined values are rejected by attribute validation. */
    void setFFT2D5Model(const std::string& model);
    static std::string getPipeMode(const ObjectValues& values);
    static double getBeamRadius(const ObjectValues& values);
    static bool getClosedRing(const ObjectValues& values);
    static bool getScatterLongitudinally(const ObjectValues& values);
    static double getPipeSizeX(const ObjectValues& values);
    static double getPipeSizeY(const ObjectValues& values);
    static std::string getRefPathFileName(const ObjectValues& values);
    void setPipeMode(const std::string& pipeMode);
    void setBeamRadius(double beamRadius);
    void setClosedRing(bool closedRing);
    void setScatterLongitudinally(bool val);
    void setPipeSizeX(double pipeSizeX);
    void setPipeSizeY(double pipeSizeY);
    void setRefPathFileName(const std::string& refPathFileName);

private:
    // Not implemented.
    FieldSolverCmd(const FieldSolverCmd&);
    void operator=(const FieldSolverCmd&);
};

#endif  // OPAL_FieldSolverCmd_HH
