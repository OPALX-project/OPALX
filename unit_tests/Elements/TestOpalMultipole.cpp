//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include <mpi.h>

#include "Ippl.h"

#include "Attributes/Attributes.h"
#include "Elements/OpalMultipole.h"
#include "Utilities/OpalException.h"
#include "ValueDefinitions/RealVariable.h"
#include "gtest/gtest.h"

#include <memory>
#include <string>
#include <vector>

// MULTIPOLE is a straight element from quadrupole up to decapole: the first (dipole) entry of
// KN, KS, DKN and DKS must be 0, and each array has at most 5 entries.

namespace {

    class TestOpalMultipole : public testing::Test {
    protected:
        // The element allocates Kokkos views for its coefficients.
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        void SetUp() override {
            // update() reads OpalData::getP0(), which needs the P0 variable (see TestOpalBend).
            [[maybe_unused]] static RealVariable* p0Owner = new RealVariable();
        }

        /// Update a clone of the exemplar (the checks skip the builtin prototype) with the given
        /// arrays; an empty array is left unset. Returns the error text, or "" if none.
        std::string updateError(
                const std::vector<double>& kn, const std::vector<double>& ks = {},
                const std::vector<double>& dkn = {}) {
            std::unique_ptr<OpalMultipole> mult(exemplar_m.clone("M1"));
            Attributes::setReal(*mult->findAttribute("L"), 0.3);
            Attributes::setReal(*mult->findAttribute("ELEMEDGE"), 0.5);
            Attributes::setRealArray(*mult->findAttribute("KN"), kn);
            if (!ks.empty()) {
                Attributes::setRealArray(*mult->findAttribute("KS"), ks);
            }
            if (!dkn.empty()) {
                Attributes::setRealArray(*mult->findAttribute("DKN"), dkn);
            }
            try {
                mult->update();
            } catch (const OpalException& e) {
                return e.what();
            }
            return "";
        }

    private:
        /// Outlives every clone made from it.
        OpalMultipole exemplar_m;
    };

}  // namespace

TEST_F(TestOpalMultipole, AcceptsQuadrupoleUpToDecapole) {
    EXPECT_EQ(updateError({0, 2, 1, 0.5, 0.2}, {0, 0.3}), "");
}

TEST_F(TestOpalMultipole, RejectsDipoleEntry) {
    EXPECT_NE(updateError({0.05}).find("first KN entry is the dipole"), std::string::npos);
    EXPECT_NE(updateError({0, 2}, {0.1, 0}).find("first KS entry"), std::string::npos);
    EXPECT_NE(updateError({0, 2}, {}, {0.01}).find("first DKN entry"), std::string::npos);
}

TEST_F(TestOpalMultipole, RejectsOrdersAboveDecapole) {
    EXPECT_NE(updateError({0, 2, 0, 0, 0, 0}).find("KN has 6 entries"), std::string::npos);
}
