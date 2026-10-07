#include "Fields/FM3DDynamic.h"

#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

namespace {
    class MapReader {
    public:
        explicit MapReader(const std::string& filename) : filename_m(filename), file_m(filename) {
            if (!file_m) fail("Cannot open field map");
        }

        bool next(std::string& record) {
            while (std::getline(file_m, record)) {
                ++line_m;
                const auto comment = record.find('#');
                if (comment != std::string::npos) record.erase(comment);
                if (record.find_first_not_of(" \t\r\n") != std::string::npos) return true;
            }
            if (file_m.bad()) fail("Error reading field map");
            return false;
        }

        std::string record() {
            std::string result;
            if (!next(result)) fail("Unexpected end of field map");
            return result;
        }

        template <typename... Values>
        void read(Values&... values) {
            std::istringstream stream(record());
            if (!(stream >> ... >> values)) fail("Invalid field map record");
            stream >> std::ws;
            if (!stream.eof()) fail("Unexpected extra values in field map record");
        }

        [[noreturn]] void fail(const std::string& message) const {
            throw GeneralOpalException(
                    "FM3DDynamic",
                    message + " in '" + filename_m + "' at line " + std::to_string(line_m));
        }

    private:
        std::string filename_m;
        std::ifstream file_m;
        std::size_t line_m = 0;
    };
}  // namespace

FM3DDynamic::FM3DDynamic(const std::string& filename) : Fieldmap(filename) {
    Type = T3DDynamic;
    MapReader reader(filename);
    std::istringstream descriptor(reader.record());
    std::string identifier, normalization, extra;
    descriptor >> identifier;
    if (identifier != "3DDynamic") reader.fail("Expected 3DDynamic header");
    if (descriptor >> normalization) {
        normalization = Util::toUpper(normalization);
        if (normalization != "TRUE" && normalization != "FALSE") {
            reader.fail("Normalization must be TRUE or FALSE");
        }
        normalize_m = normalization == "TRUE";
        if (descriptor >> extra) reader.fail("Unexpected extra values in header");
    }

    reader.read(frequency_m);
    frequency_m *= Physics::two_pi * Units::MHz2Hz;
    if (!std::isfinite(frequency_m) || frequency_m < 0.0) {
        reader.fail("Frequency must be finite and nonnegative");
    }

    const auto readAxis = [&](double& begin, double& end, double& spacing, int& points) {
        int intervals;
        reader.read(begin, end, intervals);
        begin *= Units::cm2m;
        end *= Units::cm2m;
        if (!std::isfinite(begin) || !std::isfinite(end) || !(end > begin) || intervals < 1
            || intervals == std::numeric_limits<int>::max()) {
            reader.fail("Each axis needs finite increasing bounds and a positive interval count");
        }
        spacing = (end - begin) / intervals;
        if (!std::isfinite(spacing) || !(spacing > 0.0)) reader.fail("Invalid grid spacing");
        points = intervals + 1;
    };
    readAxis(xbegin_m, xend_m, hx_m, num_gridpx_m);
    readAxis(ybegin_m, yend_m, hy_m, num_gridpy_m);
    readAxis(zbegin_m, zend_m, hz_m, num_gridpz_m);

    std::size_t size = 1;
    for (const int points : {num_gridpx_m, num_gridpy_m, num_gridpz_m}) {
        if (size > std::numeric_limits<std::size_t>::max() / (6 * sizeof(double)) / points) {
            reader.fail("Grid size exceeds addressable field storage");
        }
        size *= points;
    }
}

FM3DDynamic::~FM3DDynamic() = default;

void FM3DDynamic::readMap() {
    if (FieldstrengthEz_m.extent(0) != 0) return;

    MapReader reader(Filename_m);
    for (int headerLine = 0; headerLine < 5; ++headerLine)
        reader.record();

    const std::size_t size = static_cast<std::size_t>(num_gridpx_m) * num_gridpy_m * num_gridpz_m;
    const std::array<std::string, 6> labels{"FM3DDynamic::Ex", "FM3DDynamic::Ey",
                                            "FM3DDynamic::Ez", "FM3DDynamic::Bx",
                                            "FM3DDynamic::By", "FM3DDynamic::Bz"};
    std::array<Kokkos::DualView<double*>, 6> fields;
    for (std::size_t component = 0; component < fields.size(); ++component) {
        fields[component] = Kokkos::DualView<double*>(labels[component], size);
    }

    for (std::size_t index = 0; index < size; ++index) {
        std::array<double, 6> values;
        reader.read(values[0], values[1], values[2], values[3], values[4], values[5]);
        for (std::size_t component = 0; component < fields.size(); ++component) {
            if (!std::isfinite(values[component])) reader.fail("Field values must be finite");
            fields[component].view_host()(index) = values[component];
        }
    }
    std::string extra;
    if (reader.next(extra)) reader.fail("Unexpected data after the last grid point");

    double peakEz = 1.0;
    if (normalize_m) {
        peakEz = 0.0;
        for (const double field : sampleOnaxisEz(fields[2].view_host())) {
            peakEz = std::max(peakEz, std::abs(field));
        }
        if (!(peakEz > 0.0) || !std::isfinite(peakEz)) {
            reader.fail(
                    "Normalization needs a finite nonzero on-axis Ez; use 3DDynamic FALSE "
                    "otherwise");
        }
    }

    for (std::size_t component = 0; component < fields.size(); ++component) {
        auto host         = fields[component].view_host();
        const double unit = component < 3 ? Units::MVpm2Vpm : Physics::mu_0;
        for (std::size_t index = 0; index < size; ++index) {
            host(index) = (host(index) / peakEz) * unit;
            if (!std::isfinite(host(index))) reader.fail("Field conversion exceeds finite range");
        }
        fields[component].modify_host();
        fields[component].sync_device();
    }

    FieldstrengthEx_m = fields[0];
    FieldstrengthEy_m = fields[1];
    FieldstrengthEz_m = fields[2];
    FieldstrengthBx_m = fields[3];
    FieldstrengthBy_m = fields[4];
    FieldstrengthBz_m = fields[5];

    Inform message("FM3DDynamic::readMap");
    message << level3 << "Read in fieldmap '" << Filename_m << "'" << endl;
}

void FM3DDynamic::freeMap() {
    FieldstrengthEx_m = Kokkos::DualView<double*>();
    FieldstrengthEy_m = Kokkos::DualView<double*>();
    FieldstrengthEz_m = Kokkos::DualView<double*>();
    FieldstrengthBx_m = Kokkos::DualView<double*>();
    FieldstrengthBy_m = Kokkos::DualView<double*>();
    FieldstrengthBz_m = Kokkos::DualView<double*>();
}

void FM3DDynamic::requireLoaded() const {
    if (FieldstrengthEz_m.extent(0) == 0) {
        throw GeneralOpalException(
                "FM3DDynamic", "Field map '" + Filename_m + "' has not been loaded");
    }
}

void FM3DDynamic::applyField(std::shared_ptr<ParticleContainer_t> particles, double, double) {
    applyRFField(particles, 1.0, 1.0, zbegin_m, zend_m);
}

void FM3DDynamic::applyRFField(
        std::shared_ptr<ParticleContainer_t> particles, double electricScale, double magneticScale,
        double startField, double endField) {
    requireLoaded();
    const double xBegin = xbegin_m, xEnd = xend_m;
    const double yBegin = ybegin_m, yEnd = yend_m;
    const double zBegin = zbegin_m, zEnd = zend_m;
    const double spacingX = hx_m, spacingY = hy_m, spacingZ = hz_m;
    const int pointsX = num_gridpx_m, pointsY = num_gridpy_m, pointsZ = num_gridpz_m;
    const auto fieldEx        = FieldstrengthEx_m.view_device();
    const auto fieldEy        = FieldstrengthEy_m.view_device();
    const auto fieldEz        = FieldstrengthEz_m.view_device();
    const auto fieldBx        = FieldstrengthBx_m.view_device();
    const auto fieldBy        = FieldstrengthBy_m.view_device();
    const auto fieldBz        = FieldstrengthBz_m.view_device();
    const auto positions      = particles->R.getView();
    const auto electricFields = particles->E.getView();
    const auto magneticFields = particles->B.getView();

    Kokkos::parallel_for(
            "FM3DDynamic::applyRFField", particles->getLocalNum(),
            KOKKOS_LAMBDA(const std::size_t particle) {
                const auto& position = positions(particle);
                if (position(0) >= xBegin && position(0) < xEnd && position(1) >= yBegin
                    && position(1) < yEnd && position(2) >= zBegin && position(2) < zEnd
                    && position(2) >= startField && position(2) < endField) {
                    Vector_t<double, 3> electricField(0.0), magneticField(0.0);
                    computeField(
                            position, electricField, magneticField, fieldEx, fieldEy, fieldEz,
                            fieldBx, fieldBy, fieldBz, spacingX, spacingY, spacingZ, xBegin, yBegin,
                            zBegin, xEnd, yEnd, zEnd, pointsX, pointsY, pointsZ);
                    electricFields(particle) += electricScale * electricField;
                    magneticFields(particle) += magneticScale * magneticField;
                }
            });
}

bool FM3DDynamic::getFieldstrength(
        const Vector_t<double, 3>& position, Vector_t<double, 3>& electricField,
        Vector_t<double, 3>& magneticField) const {
    if (!isInside(position)) return true;
    requireLoaded();
    computeField(
            position, electricField, magneticField, FieldstrengthEx_m.view_host(),
            FieldstrengthEy_m.view_host(), FieldstrengthEz_m.view_host(),
            FieldstrengthBx_m.view_host(), FieldstrengthBy_m.view_host(),
            FieldstrengthBz_m.view_host(), hx_m, hy_m, hz_m, xbegin_m, ybegin_m, zbegin_m, xend_m,
            yend_m, zend_m, num_gridpx_m, num_gridpy_m, num_gridpz_m);
    return false;
}

bool FM3DDynamic::getFieldDerivative(
        const Vector_t<double, 3>&, Vector_t<double, 3>&, Vector_t<double, 3>&,
        const DiffDirection&) const {
    throw GeneralOpalException("FM3DDynamic::getFieldDerivative", "not implemented");
}

void FM3DDynamic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void FM3DDynamic::getFieldDimensions(
        double& xBegin, double& xEnd, double& yBegin, double& yEnd, double& zBegin,
        double& zEnd) const {
    xBegin = xbegin_m;
    xEnd   = xend_m;
    yBegin = ybegin_m;
    yEnd   = yend_m;
    getFieldDimensions(zBegin, zEnd);
}

void FM3DDynamic::swap() {}

void FM3DDynamic::getInfo(Inform* message) {
    (*message) << Filename_m << " (3D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m
               << " m;" << endl;
}

double FM3DDynamic::getFrequency() const { return frequency_m; }

void FM3DDynamic::setFrequency(double frequency) {
    if (!std::isfinite(frequency) || frequency < 0.0) {
        throw GeneralOpalException(
                "FM3DDynamic::setFrequency", "Frequency must be finite and nonnegative");
    }
    frequency_m = frequency;
}

std::vector<double> FM3DDynamic::sampleOnaxisEz(
        const Kokkos::DualView<double*>::t_host& fieldEz) const {
    if (!(xbegin_m <= 0.0 && xend_m >= 0.0 && ybegin_m <= 0.0 && yend_m >= 0.0)) {
        throw GeneralOpalException(
                "FM3DDynamic::sampleOnaxisEz",
                "Field map '" + Filename_m + "' does not contain x = y = 0; "
                "on-axis sampling is unavailable and normalization requires 3DDynamic FALSE");
    }
    const double gridX          = std::clamp(-xbegin_m / hx_m, 0.0, double(num_gridpx_m - 1));
    const double gridY          = std::clamp(-ybegin_m / hy_m, 0.0, double(num_gridpy_m - 1));
    const int indexX            = std::min(static_cast<int>(gridX), num_gridpx_m - 2);
    const int indexY            = std::min(static_cast<int>(gridY), num_gridpy_m - 2);
    const double fractionX      = gridX - indexX;
    const double fractionY      = gridY - indexY;
    const std::size_t strideX   = static_cast<std::size_t>(num_gridpy_m) * num_gridpz_m;
    const std::size_t strideY   = static_cast<std::size_t>(num_gridpz_m);
    const std::size_t baseIndex = static_cast<std::size_t>(indexX) * strideX + indexY * strideY;
    std::vector<double> result(num_gridpz_m, 0.0);
    for (int corner = 0; corner < 4; ++corner) {
        const int offsetX = (corner >> 1) & 1;
        const int offsetY = corner & 1;
        const double weight =
                (offsetX ? fractionX : 1.0 - fractionX) * (offsetY ? fractionY : 1.0 - fractionY);
        const std::size_t index = baseIndex + offsetX * strideX + offsetY * strideY;
        for (int indexZ = 0; indexZ < num_gridpz_m; ++indexZ) {
            result[indexZ] += weight * fieldEz(index + indexZ);
        }
    }
    return result;
}

void FM3DDynamic::getOnaxisEz(std::vector<std::pair<double, double>>& field) {
    requireLoaded();
    const auto samples = sampleOnaxisEz(FieldstrengthEz_m.view_host());
    field.resize(samples.size());
    for (std::size_t index = 0; index < samples.size(); ++index) {
        field[index] = {hz_m * index, samples[index] / Units::MVpm2Vpm};
    }
}
