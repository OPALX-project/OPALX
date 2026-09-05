//
// Class Ring
//   The RING definition. A RING uses LINE list syntax, but declares a closed
//   topology and does not permit repetition or reflection modifiers.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
#ifndef OPAL_Ring_HH
#define OPAL_Ring_HH

#include "Lines/Line.h"

class Ring : public Line {
public:
    /// Exemplar constructor.
    Ring();

    ~Ring() override = default;

    /// Make an empty clone which is filled by the parser.
    Ring* clone(const std::string& name) override;

    /// Make a complete copy, including the sequence members.
    Ring* copy(const std::string& name) override;

    /// Parse a closed-topology sequence definition.
    void parse(Statement& stat) override;

    /// Isolate and tag all occurrences owned by this RING.
    void prepareForTracking() override;

    /**
     * @brief Design circumference in metres, summed over nominal occurrence arc lengths.
     *
     * Drifts fill the nominal gaps and are counted once; field-support extensions never
     * contribute. Rectangular-bend chord lengths are converted by Geometry::getArcLength().
     * This is not the length of a tracked or closed orbit and does not validate 3D closure.
     */
    double getLength() const override;

protected:
    const char* getSequenceKeyword() const override;

private:
    Ring(const std::string& name, Ring* parent);

    Ring(const Ring&)           = delete;
    void operator=(const Ring&) = delete;
};

#endif  // OPAL_Ring_HH
