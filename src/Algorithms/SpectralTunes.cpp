// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/SpectralTunes.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "AbsBeamline/CyclotronSector.h"
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"
extern Inform* gmsg;
namespace {
    void require(bool ok, const std::string& message) {
        if (!ok) throw OpalException("SpectralTunes", message);
    }
}  // namespace
SpectralTunes::Spectrum SpectralTunes::analyze(const std::vector<double>& values, unsigned turns) {
    require(values.size() >= 12 && turns > 0,
            "At least 12 samples and positive turns are required.");
    double mean = 0;
    for (double v : values) {
        require(std::isfinite(v), "Non-finite signal.");
        mean += v;
    }
    mean /= values.size();
    double ss = 0, residual = 0;
    for (double v : values) {
        const double d = v - mean;
        ss += d * d;
        residual += d;
    }
    const double variance = (ss - residual * residual / values.size()) / (values.size() - 1);
    require(std::isfinite(variance) && variance > 0, "Zero or invalid signal variance.");
    Spectrum result{{0, -1, 1.0 / (4.0 * turns)}, {}, {}};
    const size_t n = values.size(), bins = static_cast<size_t>(1.6 * n);
    for (size_t k = 1; k <= bins; ++k) {
        const double omega = 2 * std::acos(-1.0) * k / (4 * (n - 1)), dc = std::cos(omega),
                     ds = std::sin(omega);
        double c = 0, s = 0, cc = 0, cs = 0, sn = 0, yc = 0, ys = 0;
        for (size_t i = 0; i < n; ++i) {
            // Reanchor recurrence to bound roundoff in sinusoidal least squares.
            if (i % 256 == 0) {
                c = std::cos(omega * (double(i) - 0.5 * (n - 1)));
                s = std::sin(omega * (double(i) - 0.5 * (n - 1)));
            }
            cc += c * c;
            sn += s * s;
            cs += c * s;
            yc += (values[i] - mean) * c;
            ys += (values[i] - mean) * s;
            const double next = c * dc - s * ds;
            s                 = s * dc + c * ds;
            c                 = next;
        }
        const double determinant = cc * sn - cs * cs;
        require(determinant > 0, "Singular sinusoidal least-squares basis.");
        const double power =
                (sn * yc * yc - 2 * cs * yc * ys + cc * ys * ys) / (2 * variance * determinant);
        const double tune = k * result.peak.spacing;
        result.tune.push_back(tune);
        result.power.push_back(power);
        if (power >= result.peak.power) result.peak = {tune, power, result.peak.spacing};
    }
    return result;
}
void SpectralTunes::run(
        OpalBeamline& beamline, const PartData& reference, const std::vector<double>& triples,
        const Settings& settings, const std::string& base) {
    require(!triples.empty() && triples.size() % 3 == 0,
            "TUNEINITIAL requires (kinetic GeV, radius m, pr/mc) triples.");
    require(settings.dt > 0 && std::isfinite(settings.dt) && settings.turns > 0
                    && settings.stepsPerTurn > 0 && settings.sampleEvery > 0,
            "Invalid sampling controls.");
    require(reference.getM() > 0 && std::isfinite(reference.getM()) && reference.getQ() > 0,
            "Spectral rays require positive mass and charge.");
    require(std::isfinite(settings.radialOffset) && std::isfinite(settings.verticalOffset)
                    && settings.radialOffset != 0 && settings.verticalOffset != 0,
            "Both tune displacements must be finite and nonzero.");
    std::shared_ptr<ElementBase> anchor;
    for (const auto& element : beamline.getElements()) {
        require(element->getType() == ElementType::CYCLOTRONSECTOR,
                "Only CYCLOTRONSECTOR elements are supported (no RF).");
        if (element->getName() == settings.sector) anchor = element;
    }
    require(bool(anchor), "TUNESECTOR was not found: " + settings.sector);
    const auto& pose          = anchor->getCSTrafoGlobal2Local();
    const double designRadius = 1 / anchor->getGeometry().getCurvature();
    const auto centre         = pose.transformFrom(Vector_t<double, 3>(-designRadius, 0, 0));
    const auto radial         = pose.rotateFrom(Vector_t<double, 3>(1, 0, 0));
    const auto tangent        = pose.rotateFrom(Vector_t<double, 3>(0, 0, 1));
    const auto vertical       = pose.rotateFrom(Vector_t<double, 3>(0, 1, 0));
    ExternalFieldRayTracker tracker(beamline, reference, settings.integrator);
    std::ofstream summary(base + "-tunes.csv");
    require(bool(summary), "Cannot create tune summary.");
    summary << "energy_MeV,radius_m,pr,nu_r,nu_y,grid_spacing,legacy_turns,samples,max_relative_"
               "energy_drift\n"
            << std::setprecision(16);
    const size_t steps = size_t(settings.turns) * settings.stepsPerTurn;
    for (size_t row = 0; row < triples.size() / 3; ++row) {
        const double energy = triples[3 * row], radius = triples[3 * row + 1],
                     pr    = triples[3 * row + 2];
        const double gamma = 1 + energy * 1e9 / reference.getM(), p2 = gamma * gamma - 1;
        require(std::isfinite(energy + radius + pr) && std::isfinite(p2) && energy > 0 && radius > 0
                        && pr * pr < p2,
                "Invalid TUNEINITIAL row.");
        ExternalFieldRayTracker::State a, b;
        a.position = centre + radius * radial;
        a.momentum = pr * radial + std::sqrt(p2 - pr * pr) * tangent;
        b          = a;
        b.position += settings.radialOffset * radial + settings.verticalOffset * vertical;
        std::vector<double> r, y;
        std::ofstream samples(base + "-tune-" + std::to_string(row) + "-samples.csv");
        require(bool(samples), "Cannot create tune samples.");
        samples << "time_s,delta_r_m,y_m,reference_angle_rad\n" << std::setprecision(16);
        double angle = 0, oldAngle = 0, maxDrift = 0;
        unsigned sampledTurns = 1;
        for (size_t step = 0; step < steps; ++step) {
            require(!beamline.getElements(a.position).empty()
                            && !beamline.getElements(b.position).empty(),
                    "A tune ray left sector support.");
            const Vector_t<double, 3> ra = a.position - centre, rb = b.position - centre;
            const double dr = std::hypot(dot(rb, radial), dot(rb, tangent))
                              - std::hypot(dot(ra, radial), dot(ra, tangent));
            const double z                 = dot(rb, vertical);
            a                              = tracker.advance(a, settings.dt);
            b                              = tracker.advance(b, settings.dt);
            const Vector_t<double, 3> next = a.position - centre;
            const double theta             = std::atan2(dot(next, tangent), dot(next, radial));
            angle += std::remainder(theta - oldAngle, 2 * std::acos(-1.0));
            oldAngle = theta;
            require(std::isfinite(angle) && angle >= -1e-8,
                    "Unexpected reference rotation direction.");
            for (const auto& ray : {a, b})
                maxDrift = std::max(
                        maxDrift, std::abs(
                                          (std::sqrt(1 + dot(ray.momentum, ray.momentum)) - gamma)
                                          / (gamma - 1)));
            // Old SEO saves pre-step motion and the post-step turn count, initially 1.
            if (step % settings.sampleEvery == 0) {
                r.push_back(dr);
                y.push_back(z);
                // Legacy isTurnDone() counts nominal steps, not geometric crossings.
                sampledTurns = 1 + static_cast<unsigned>((step + 1) / settings.stepsPerTurn);
                samples << step * settings.dt << ',' << dr << ',' << z << ',' << angle << '\n';
            }
        }
        require(!beamline.getElements(a.position).empty()
                        && !beamline.getElements(b.position).empty(),
                "A tune ray left sector support on the final step.");
        auto sr = analyze(r, sampledTurns), sy = analyze(y, sampledTurns);
        std::ofstream spectrum(base + "-tune-" + std::to_string(row) + "-spectrum.csv");
        require(bool(spectrum), "Cannot create tune spectrum.");
        spectrum << "tune,radial_power,vertical_power\n" << std::setprecision(16);
        for (size_t i = 0; i < sr.tune.size(); ++i)
            spectrum << sr.tune[i] << ',' << sr.power[i] << ',' << sy.power[i] << '\n';
        summary << energy * 1000 << ',' << radius << ',' << pr << ',' << sr.peak.tune << ','
                << sy.peak.tune << ',' << sr.peak.spacing << ',' << sampledTurns << ',' << r.size()
                << ',' << maxDrift << std::endl;
        *gmsg << "SpectralTunes K=" << energy * 1000 << " MeV nu_r=" << sr.peak.tune
              << " nu_y=" << sy.peak.tune << " turns=" << sampledTurns << endl;
    }
}
