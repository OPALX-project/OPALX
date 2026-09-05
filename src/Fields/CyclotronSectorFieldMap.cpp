#include "Fields/CyclotronSectorFieldMap.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include "Utilities/OpalException.h"

std::shared_ptr<const CyclotronSectorFieldMap> CyclotronSectorFieldMap::read(
        const std::string& filename) {
    static std::mutex mutex;
    static std::map<std::string, std::weak_ptr<const CyclotronSectorFieldMap>> cache;
    std::lock_guard<std::mutex> lock(mutex);
    const auto key = std::filesystem::canonical(filename).string();
    if (auto found = cache[key].lock()) return found;
    auto map   = std::shared_ptr<const CyclotronSectorFieldMap>(new CyclotronSectorFieldMap(key));
    cache[key] = map;
    return map;
}

CyclotronSectorFieldMap::CyclotronSectorFieldMap(const std::string& filename) {
    std::ifstream in(filename);
    auto fail = [&]() {
        throw OpalException("CyclotronSectorFieldMap", "Invalid PSI map: " + filename);
    };
    auto skip = [&](int n) {
        std::string token;
        for (int i = 0; i < n; ++i)
            if (!(in >> token)) fail();
    };
    if (!(in >> rmin >> dr >> thetaMin >> dtheta)) fail();
    if (dr < 0) dr = -1 / dr;
    if (dtheta < 0) dtheta = -1 / dtheta;
    rmin *= 0.001;
    dr *= 0.001;
    thetaMin *= std::acos(-1.0) / 180;
    dtheta *= std::acos(-1.0) / 180;
    skip(13);
    if (!(in >> nr >> nt) || nr < 5 || nt < 5 || nr > 100000 || nt > 100000
        || static_cast<long long>(nr) * nt > 100000000 || !(rmin > 0) || !(dr > 0) || !(dtheta > 0)
        || !std::isfinite(rmin + dr + thetaMin + dtheta))
        fail();
    skip(5);
    int count;
    if (!(in >> count) || count < 0 || count > 10000000) fail();
    skip(4);
    for (int i = 0; i < count; ++i) {
        double unused;
        if (!(in >> unused)) fail();
    }
    skip(6);
    std::string token;
    bool found = false;
    for (int i = 0; i < 10000 && in >> token; ++i)
        if (token == "LREC=") {
            found = true;
            break;
        }
    if (!found) fail();
    skip(5);
    data = View("cyclotron_map", nr, nt + 1, 3);
    host = Kokkos::create_mirror_view(data);
    for (int i = 0; i < nr; ++i) {
        if (i) skip(6);
        for (int channel = 0; channel < 4; ++channel)
            for (int j = 0; j < nt; ++j) {
                double value;
                if (!(in >> value) || !std::isfinite(value)) fail();
                if (channel == 0) host(i, j, 0) = value * 0.1;  // kG -> T
            }
    }
    // Legacy OPAL uses one-sided stencils near both boundaries, including theta.
    const double weights[5][5] = {
            {-50, 96, -72, 32, -6},
            {-6, -20, 36, -12, 2},
            {2, -16, 0, 16, -2},
            {-2, 12, -36, 20, 6},
            {6, -32, 72, -96, 50}};
    for (int i = 0; i < nr; ++i)
        for (int j = 0; j < nt; ++j) {
            int ri = std::clamp(i - 2, 0, nr - 5), tj = std::clamp(j - 2, 0, nt - 5);
            double radial = 0, angular = 0;
            for (int k = 0; k < 5; ++k) {
                radial += weights[i - ri][k] * host(ri + k, j, 0);
                angular += weights[j - tj][k] * host(i, tj + k, 0);
            }
            host(i, j, 1) = radial / (24 * dr);
            host(i, j, 2) = angular / (24 * dtheta);
        }
    for (int i = 0; i < nr; ++i)
        for (int k = 0; k < 3; ++k)
            host(i, nt, k) = host(i, 0, k);
    Kokkos::deep_copy(data, host);
}
