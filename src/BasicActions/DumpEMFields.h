/*
 *  Copyright (c) 2017, Chris Rogers
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

#ifndef OPAL_BASICACTIONS_DUMPEMFIELDS_HH
#define OPAL_BASICACTIONS_DUMPEMFIELDS_HH

#include <unordered_set>
#include <string>
#include "AbstractObjects/Action.h"

namespace interpolation {
  class NDGrid;
}

class Component;

/** DumpEMFields dumps the dynamically changing fields of a Ring in a user-
 *  defined grid.
 *
 *  The idea is to print out the field map across a 4D grid in space-time for
 *  debugging purposes. The problem is to manage the DumpEMFields object through
 *  three phases of program execution; initial construction, parsing and then
 *  actual field map writing (where we need to somehow let DumpFields know what
 *  the field maps are). So for each  DumpFields object created, we store in a
 *  set. When the execute() method is called, DumpFields builds a grid using
 *  the parsed information.
 *
 *  When the ParallelCyclotronTracker is about to start tracking, it calls
 *  writeFields method which loops over the static set of DumpFields and writes
 *  each one. It is not the cleanest implementation, but I can't see a better
 *  way.
 *
 *  The DumpEMFields themselves operate by iterating over a NDGrid object
 *  and looking up the field/writing it out on each grid point.
 *
 */
class DumpEMFields : public Action {
  public:
    enum CoordinateSystem {CYLINDRICAL, CARTESIAN};

    /** Constructor */
    DumpEMFields();

    /** Destructor deletes grid_m and if in the dumps set, take it out */
    virtual ~DumpEMFields();

    /** Make a clone (overloadable copy-constructor).
     *    @param name not used
     *  If this is in the dumpsSet_m, so will the clone. Not sure how the
     *  itsAttr stuff works, so this may not get properly copied?
     */
    virtual DumpEMFields *clone(const std::string &name);

    /** Builds the grid but does not write the field map
     *
     *  Builds a grid of points in x-y-z space using the NDGrid algorithm.
     *  Checks that X_STEPS, Y_STEPS, Z_STEPS are integers or throws
     *  OpalException.
     */
    virtual void execute();

    /** Write the fields for all defined DumpEMFields objects
     *    @param field borrowed reference to the Component object that holds the
     *    field map; caller owns the memory.
     *  Iterates over the DumpEMFields in the dumpsSet_m and calls writeFieldThis
     *  on each DumpEMFields. This writes each field map in turn. Format is:
     *  <number of rows>
     *  <column 1> <units>
     *  <column 2> <units>
     *  <column 3> <units>
     *  <column 4> <units>
     *  <column 5> <units>
     *  <column 6> <units>
     *  0
     *  <field map data>
     */
    static void writeFields(Component* field);

  private:
    virtual void writeFieldThis(Component* field);
    virtual void buildGrid();
    void parseCoordinateSystem();
    static void checkInt(double value, std::string name, double tolerance = 1e-9);
    void writeHeader(std::ofstream& fout) const;
    void writeFieldLine(Component* field,
                        const Vector_t& point,
                        const double& time,
                        std::ofstream& fout) const;

    interpolation::NDGrid* grid_m;
    std::string filename_m;
    CoordinateSystem coordinates_m = CARTESIAN;

    static std::unordered_set<DumpEMFields*> dumpsSet_m;
    static std::string dumpemfields_docstring;

    constexpr static double DEGREE = 180./M_PI;
    DumpEMFields(const DumpEMFields& dump);  // disabled
    DumpEMFields& operator=(const DumpEMFields& dump);  // disabled
};

#endif // ifdef OPAL_DUMPFIELDS_HH

