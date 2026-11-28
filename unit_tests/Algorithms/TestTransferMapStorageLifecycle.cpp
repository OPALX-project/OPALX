#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>
#include "Algorithms/LinearTransferMapBuilder.h"
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Physics/Physics.h"
#include "Structure/LinearTransferMap.h"
#include "Utilities/OpalException.h"

int main() {
    try {
        LinearTransferMapHandle retained;
        std::shared_ptr<const TransferMapPass> pass;
        {
            Kokkos::ScopeGuard runtime;
            LinearTransferMapDraft draft;
            draft.entrance.momentum = draft.exit.momentum = {0, 0, 1};
            for (unsigned i = 0; i < 6; ++i)
                draft.matrix(i, i) = 1;
            retained = std::make_shared<const LinearTransferMap>(draft, "explicit-release");
            if (!retained->available() || retained->matrix()(0, 0) != 1) return 1;
            releaseTransferMapStorage();
            releaseTransferMapStorage();
            if (retained->available()) return 2;
            bool rejected = false;
            try {
                (void)retained->matrix();
            } catch (const OpalException&) {
                rejected = true;
            }
            if (!rejected) return 3;
            // Fresh payloads after explicit release must still be covered by the finalize hook.
            retained = std::make_shared<const LinearTransferMap>(draft, "finalize-release");
            OpalBeamline beamline;
            PartData particle(1, 9.382720813e8, 1e6);
            auto entrance    = draft.entrance;
            auto exit        = entrance;
            exit.position(2) = exit.pathLength = .1;
            exit.time                          = .1 * std::sqrt(2.) / Physics::c;
            LinearTransferMapBuilder builder(beamline, particle, 1e-11);
            pass = builder.build({{entrance}, {exit}}, 0).pass;
            if (!pass || !pass->available()) return 4;
        }
        if (retained->available() || pass->available() || pass->segments().front()->available())
            return 5;
        if (retained->id() != "finalize-release" || retained->metadata().entrance.momentum(2) != 1)
            return 6;
        bool rejected = false;
        try {
            (void)pass->matrix();
        } catch (const OpalException&) {
            rejected = true;
        }
        if (!rejected) return 7;
        pass.reset();
        retained.reset();  // Must not destroy any Kokkos allocation after finalize.
        std::cout << "Persistent map teardown and late handle destruction passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 8;
    }
}
