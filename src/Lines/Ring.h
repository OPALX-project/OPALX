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

/**
 * @brief RING sequence declaration: ordered lattice members with periodic topology.
 *
 * A RING is not a physical element. opalx::run::describeLine() gives it the RING topology,
 * and opalx::run::buildElements() gives its elements the RING membership; TrackRun then
 * tracks periodically. This declaration alone does not enforce geometrical closure or find
 * a closed orbit. The circumference is opalx::run::ringLength() of the built elements.
 */
class Ring : public Line {
public:
    /// Exemplar constructor.
    Ring();

    ~Ring() override = default;

    /// Make an empty clone which is filled by the parser.
    Ring* clone(const std::string& name) override;

    /// Make complete copy, the member list included.
    Ring* copy(const std::string& name) override;

    /// Parse a closed-topology sequence definition.
    void parse(Statement& stat) override;

protected:
    const char* getSequenceKeyword() const override;

private:
    Ring(const std::string& name, Ring* parent);

    Ring(const Ring&)           = delete;
    void operator=(const Ring&) = delete;
};

#endif  // OPAL_Ring_HH
