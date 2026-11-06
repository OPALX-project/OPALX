/**
 * Standalone FFT2D5 backend experiment; no production solver or IPPL changes.
 * Compare the actual solveBatched() with existing heFFTe public plans and
 * backend executor factories. Only x and y are transformed; z is a batch index.
 * Setup/warmup is excluded. Validate all potential halos, all E components,
 * empty/vacated slices, neutral sources, remainders and rectangular grids.
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include "Physics/Physics.h"
#include "SpaceCharge/FFT2D5/FFT2D5Poisson.h"

namespace {
    using P       = FFT2D5Poisson;
    using C       = Kokkos::complex<double>;
    using Mem     = typename P::ScalarField2_t::view_type::memory_space;
    using HB      = ippl::fft::HeffteR2C<double, 2, Mem>;
    using Backend = typename HB::backend_t;
    using Exec    = typename heffte::one_dim_backend<Backend>::executor;
    using RExec   = typename heffte::one_dim_backend<Backend>::executor_r2c;
    using Box     = heffte::box3d<long long>;
    using CV      = Kokkos::View<C*, Mem>;
    using RV      = Kokkos::View<double*, Mem>;
    using IV      = Kokkos::View<size_t*, Mem>;
    using Policy  = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
    Policy range(size_t a, size_t b, size_t c) {
        return Policy(
                {0, 0, 0},
                {static_cast<long long>(a), static_cast<long long>(b), static_cast<long long>(c)});
    }
    using Clock  = std::chrono::steady_clock;
    using Stream = typename heffte::backend::device_instance<
            typename heffte::backend::buffer_traits<Backend>::location>::stream_type;
    Stream stream() {
#ifdef KOKKOS_ENABLE_CUDA
        return Kokkos::DefaultExecutionSpace().cuda_stream();
#else
        return nullptr;
#endif
    }
    Box box(size_t x, size_t y, size_t z) {
        return Box(
                {0, 0, 0}, {static_cast<long long>(x - 1), static_cast<long long>(y - 1),
                            static_cast<long long>(z - 1)});
    }
    std::complex<double>* ptr(CV v) {
        static_assert(sizeof(C) == sizeof(std::complex<double>));
        return reinterpret_cast<std::complex<double>*>(v.data());
    }
    enum Mode { PublicReorder, PublicNoReorder, RowsTranspose, MergedC2C };
    const char* label(Mode m) {
        switch (m) {
            case PublicReorder:
                return "public_reorder_scale_fused";
            case PublicNoReorder:
                return "public_no_reorder_scale_fused";
            case RowsTranspose:
                return "executor_r2c_rows_transpose";
            default:
                return "executor_c2c_merged_xy";
        }
    }

    struct Plans {
        std::unique_ptr<typename HB::heffte_t> publicPlan;
        std::unique_ptr<RExec> x;
        std::unique_ptr<Exec> y, xy;
        typename HB::workspace_t workspace;
        CV scratch;
        Plans(Mode mode, size_t nx2, size_t ny2, int count) {
            size_t ncx = nx2 / 2 + 1;
            if (mode == PublicReorder || mode == PublicNoReorder) {
                auto opts          = heffte::default_options<Backend>();
                opts.use_pencils   = true;
                opts.use_reorder   = mode == PublicReorder;
                opts.use_gpu_aware = true;
                publicPlan         = std::make_unique<typename HB::heffte_t>(
                        box(nx2, ny2, 1), box(ncx, ny2, 1), 0, MPI_COMM_SELF, opts);
                workspace = typename HB::workspace_t(publicPlan->size_workspace() * count);
            } else if (mode == RowsTranspose) {
                // Flatten independent plane rows into one executor batch. Columns are
                // transposed into contiguous y rows, also batched over every plane.
                x = heffte::make_executor_r2c<Backend>(stream(), box(nx2, ny2 * count, 1), 0);
                y = heffte::make_executor<Backend>(stream(), box(ny2, ncx * count, 1), 0);
                scratch =
                        CV("executorScratch",
                           std::max({size_t(1), x->workspace_size(), y->workspace_size()}));
            } else {
                // This is a heFFTe 2-D executor with count planes, not a 3-D FFT.
                xy      = heffte::make_executor<Backend>(stream(), box(nx2, ny2, count), 0, 1);
                scratch = CV("executorScratch", std::max(size_t(1), xy->workspace_size()));
            }
        }
    };

    struct Candidate {
        Mode mode;
        size_t nx, ny, nx2, ny2, ncx, nr, nc;
        int capacity;
        double area, inv2dx, inv2dy;
        RV real;
        CV spectrum, columns, full, greenRow, greenColumn, greenFull;
        IV slices;
        std::map<int, std::unique_ptr<Plans>> plans;
        Candidate(
                Mode m, size_t x, size_t y, int b, const P::ComplexField2_t& normalizedGreen,
                double dx, double dy)
            : mode(m),
              nx(x),
              ny(y),
              nx2(2 * x),
              ny2(2 * y),
              ncx(x + 1),
              nr(nx2 * ny2),
              nc(ncx * ny2),
              capacity(b),
              area(dx * dy),
              inv2dx(.5 / dx),
              inv2dy(.5 / dy),
              real("batchReal", nr * b),
              spectrum("batchSpectrum", nc * b),
              greenRow("greenRow", nc),
              greenColumn("greenColumn", nc),
              greenFull("greenFull", nr) {
            if (mode == RowsTranspose) columns = CV("batchColumns", nc * b);
            if (mode == MergedC2C) full = CV("batchFull", nr * b);
            initializeGreen(normalizedGreen);
        }
        void initializeGreen(const P::ComplexField2_t& normalizedGreen) {
            auto g       = normalizedGreen.getView();
            const int ng = normalizedGreen.getNghost();
            auto gr = greenRow, gc = greenColumn, gf = greenFull;
            const size_t cx = ncx, xx = nx2, yy = ny2;
            Kokkos::parallel_for(
                    "benchmarkGreenLayouts", nr, KOKKOS_LAMBDA(size_t p) {
                        size_t i = p % xx, j = p / xx;
                        if (i < cx) {
                            gr(i + cx * j) = g(i + ng, j + ng);
                            gc(i * yy + j) = g(i + ng, j + ng);
                        }
                        if (i < cx)
                            gf(p) = g(i + ng, j + ng);
                        else {
                            auto value = g(xx - i + ng, (yy - j) % yy + ng);
                            gf(p)      = C(value.real(), -value.imag());
                        }
                    });
            Kokkos::fence();
        }
        template <class Work>
        static void phase(const char* name, Work&& work) {
            opalx::spacecharge::fft2d5_detail::timed(name, [&] {
                work();
                // Match the explicit phase fence inside the current solveBatched().
                Kokkos::fence();
            });
        }
        Plans& plan(int count) {
            auto& p = plans[count];
            if (!p) p = std::make_unique<Plans>(mode, nx2, ny2, count);
            return *p;
        }
        void solve(
                const P::ScalarField3_t& density, P::VectorField3_t& electric,
                std::span<const size_t> occupied, double coupling, P::ScalarField3_t* potential) {
            Kokkos::deep_copy(electric.getView(), P::Vector3_t(0.));
            if (potential) Kokkos::deep_copy(potential->getView(), 0.);
            if (occupied.empty()) return;
            if (slices.extent(0) < occupied.size()) slices = IV("benchmarkSlices", occupied.size());
            auto host = Kokkos::create_mirror_view(slices);
            for (size_t i = 0; i < occupied.size(); ++i)
                host(i) = occupied[i];
            Kokkos::deep_copy(slices, host);
            const auto source = density.getView();
            auto e            = electric.getView();
            const int sg = density.getNghost(), eg = electric.getNghost();
            auto r = real;
            auto s = spectrum, c = columns, f = full;
            auto ids = slices;
            auto gr = greenRow, gc = greenColumn, gf = greenFull;
            const size_t x = nx, y = ny, xx = nx2, yy = ny2, cx = ncx, rr = nr, cc = nc;
            const bool complex = mode == MergedC2C;
            const double scale = area, idx = inv2dx, idy = inv2dy;
            for (size_t first = 0; first < occupied.size(); first += capacity) {
                int count = static_cast<int>(std::min(size_t(capacity), occupied.size() - first));
                auto& p   = plan(count);
                phase("2D5 Pack", [&] {
                    Kokkos::parallel_for(
                            "benchmarkPack", range(count, xx, yy),
                            KOKKOS_LAMBDA(size_t b, size_t i, size_t j) {
                                double value = (i < x && j < y)
                                                       ? source(i + sg, j + sg, ids(first + b) + sg)
                                                                 * coupling
                                                       : 0.;
                                if (complex)
                                    f(b * rr + i + xx * j) = C(value, 0.);
                                else
                                    r(b * rr + i + xx * j) = value;
                            });
                });
                phase("2D5 FFT Rho", [&] {
                    if (mode == MergedC2C)
                        p.xy->forward(ptr(full), ptr(p.scratch));
                    else if (mode == RowsTranspose) {
                        p.x->forward(real.data(), ptr(spectrum), ptr(p.scratch));
                        Kokkos::parallel_for(
                                "benchmarkTransposeToColumns", size_t(count) * cc,
                                KOKKOS_LAMBDA(size_t v) {
                                    size_t b = v / cc, k = v % cc, kx = k % cx, ky = k / cx;
                                    c(b * cc + kx * yy + ky) = s(v);
                                });
                        p.y->forward(ptr(columns), ptr(p.scratch));
                    } else
                        p.publicPlan->forward(
                                count, real.data(), spectrum.data(), p.workspace.data(),
                                heffte::scale::none);
                });
                phase("2D5 Convolution", [&] {
                    if (mode == MergedC2C) {
                        Kokkos::parallel_for(
                                "benchmarkMultiplyFull", size_t(count) * rr,
                                KOKKOS_LAMBDA(size_t v) { f(v) *= gf(v % rr) * scale; });
                    } else if (mode == RowsTranspose) {
                        Kokkos::parallel_for(
                                "benchmarkMultiplyColumns", size_t(count) * cc,
                                KOKKOS_LAMBDA(size_t v) { c(v) *= gc(v % cc) * scale; });
                    } else {
                        Kokkos::parallel_for(
                                "benchmarkMultiplyRows", size_t(count) * cc,
                                KOKKOS_LAMBDA(size_t v) { s(v) *= gr(v % cc) * scale; });
                    }
                });
                phase("2D5 FFT Back", [&] {
                    if (mode == MergedC2C)
                        p.xy->backward(ptr(full), ptr(p.scratch));
                    else if (mode == RowsTranspose) {
                        p.y->backward(ptr(columns), ptr(p.scratch));
                        Kokkos::parallel_for(
                                "benchmarkTransposeToRows", size_t(count) * cc,
                                KOKKOS_LAMBDA(size_t v) {
                                    size_t b = v / cc, k = v % cc, kx = k % cx, ky = k / cx;
                                    s(v) = c(b * cc + kx * yy + ky);
                                });
                        p.x->backward(ptr(spectrum), real.data(), ptr(p.scratch));
                    } else
                        p.publicPlan->backward(
                                count, spectrum.data(), real.data(), p.workspace.data(),
                                heffte::scale::none);
                });
                phase("2D5 Field", [&] {
                    Kokkos::parallel_for(
                            "benchmarkField", range(count, x, y),
                            KOKKOS_LAMBDA(size_t b, size_t i, size_t j) {
                                const size_t a  = b * rr + (i + 1) % xx + xx * j,
                                             bx = b * rr + (i + xx - 1) % xx + xx * j;
                                const size_t d  = b * rr + i + xx * ((j + 1) % yy),
                                             dy = b * rr + i + xx * ((j + yy - 1) % yy);
                                e(i + eg, j + eg, ids(first + b) + eg) = P::Vector3_t(
                                        -(complex ? f(a).real() - f(bx).real() : r(a) - r(bx))
                                                * idx,
                                        -(complex ? f(d).real() - f(dy).real() : r(d) - r(dy))
                                                * idy,
                                        0.);
                            });
                    if (potential) {
                        auto phi = potential->getView();
                        int pg   = potential->getNghost();
                        Kokkos::parallel_for(
                                "benchmarkPotential", range(count, phi.extent(0), phi.extent(1)),
                                KOKKOS_LAMBDA(size_t b, size_t i, size_t j) {
                                    size_t ix = (int(i) - pg + int(xx)) % int(xx),
                                           iy = (int(j) - pg + int(yy)) % int(yy),
                                           v  = b * rr + ix + xx * iy;
                                    phi(i, j, ids(first + b) + pg) = complex ? f(v).real() : r(v);
                                });
                    }
                });
            }
        }
    };

    struct Problem {
        size_t nx, ny, nz;
        double dx = .001, dy = .0013;
        ippl::NDIndex<3> domain;
        P::Mesh3_t mesh;
        ippl::FieldLayout<3> layout;
        P::ScalarField3_t density, potential, referencePotential;
        P::VectorField3_t field, referenceField;
        P::NDIndex2_t domain2;
        P::Mesh2_t mesh2;
        P::Layout2_t layout2;
        P::ScalarField2_t slice;
        P gold;
        std::vector<size_t> occupied;
        typename P::ScalarField3_t::view_type original;
        Problem(size_t x, size_t y, size_t z)
            : nx(x),
              ny(y),
              nz(z),
              domain{int(x), int(y), int(z)},
              mesh(domain, {dx, dy, .05}, {0., 0., 0.}),
              layout(MPI_COMM_WORLD, domain, {false, false, false}),
              density(mesh, layout),
              potential(mesh, layout),
              referencePotential(mesh, layout),
              field(mesh, layout),
              referenceField(mesh, layout),
              domain2{int(x), int(y)},
              mesh2(domain2, {dx, dy}, {0., 0.}),
              layout2(MPI_COMM_SELF, domain2, {false, false}),
              slice(mesh2, layout2) {
            gold.initialiseViews(slice);
        }
        void fill(int sourceCase) {
            occupied.clear();
            if (sourceCase == 0) {
                for (size_t k = 0; k < nz; ++k)
                    occupied.push_back(k);
            } else if (sourceCase == 1) {
                for (size_t k = 0; k < nz; ++k)
                    if (k % 11) occupied.push_back(k);
            } else if (sourceCase == 2) {
                occupied = {2, 7, nz - 1};
            }
            auto h = Kokkos::create_mirror_view(density.getView());
            Kokkos::deep_copy(h, 0.);
            for (size_t k : occupied)
                for (size_t i = 0; i < nx; ++i)
                    for (size_t j = 0; j < ny; ++j) {
                        double v = 1.e-9
                                   * ((1. + .01 * k)
                                              * std::exp(
                                                      -.03 * (double(i) - nx / 2)
                                                              * (double(i) - nx / 2)
                                                      - .05 * (double(j) - ny / 2)
                                                                * (double(j) - ny / 2))
                                      + .15 * std::sin(.7 * i + .9 * j));
                        if (sourceCase == 2 && k == 2)
                            v = (i == 0 && j == 0)             ? 1.e-9
                                : (i == nx - 1 && j == ny - 1) ? -1.e-9
                                                               : 0.;
                        h(i + 1, j + 1, k + 1) = v;
                    }
            Kokkos::deep_copy(density.getView(), h);
            original = decltype(original)("sourceCopy", nx + 2, ny + 2, nz + 2);
            Kokkos::deep_copy(original, density.getView());
        }
        void reference(int batch) {
            gold.solveBatched(
                    density, referenceField, occupied, batch, 1. / Physics::epsilon_0,
                    &referencePotential);
        }
        std::pair<double, double> check() {
            auto a = field.getView(), b = referenceField.getView();
            auto p = potential.getView(), q = referencePotential.getView();
            auto src = density.getView(), before = original;
            auto policy  = range(a.extent(0), a.extent(1), a.extent(2));
            double ediff = 0., emax = 0., pdiff = 0., pmax = 0., sourceDiff = 0.;
            Kokkos::parallel_reduce(
                    "fieldDifference", policy,
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k, double& m) {
                        for (int d = 0; d < 3; ++d)
                            m = Kokkos::max(m, Kokkos::abs(a(i, j, k)[d] - b(i, j, k)[d]));
                    },
                    Kokkos::Max<double>(ediff));
            Kokkos::parallel_reduce(
                    "fieldScale", policy,
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k, double& m) {
                        for (int d = 0; d < 3; ++d)
                            m = Kokkos::max(m, Kokkos::abs(b(i, j, k)[d]));
                    },
                    Kokkos::Max<double>(emax));
            Kokkos::parallel_reduce(
                    "potentialDifference", policy,
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k, double& m) {
                        m = Kokkos::max(m, Kokkos::abs(p(i, j, k) - q(i, j, k)));
                    },
                    Kokkos::Max<double>(pdiff));
            Kokkos::parallel_reduce(
                    "potentialScale", policy,
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k, double& m) {
                        m = Kokkos::max(m, Kokkos::abs(q(i, j, k)));
                    },
                    Kokkos::Max<double>(pmax));
            Kokkos::parallel_reduce(
                    "unchangedSource", policy,
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k, double& m) {
                        m = Kokkos::max(m, Kokkos::abs(src(i, j, k) - before(i, j, k)));
                    },
                    Kokkos::Max<double>(sourceDiff));
            double e = ediff / std::max(emax, 1.e-10), v = pdiff / std::max(pmax, 1.e-10);
            if (e > 1.e-12 || v > 1.e-12 || sourceDiff != 0.)
                throw std::runtime_error(
                        "Numerical mismatch: E=" + std::to_string(e) + " phi=" + std::to_string(v));
            return {e, v};
        }
    };

    template <class Work>
    double measure(Work work, int repeat) {
        Kokkos::fence();
        auto t = Clock::now();
        for (int i = 0; i < repeat; ++i)
            work();
        Kokkos::fence();
        return std::chrono::duration<double>(Clock::now() - t).count() / repeat;
    }
    void run(bool validateOnly, bool quick) {
        std::ofstream checks("validation.csv"), times("samples.csv");
        checks << std::setprecision(17);
        times << std::setprecision(17);
        checks << "method,nx,ny,nz,batch,source_case,occupied,relative_E_error,relative_phi_"
                  "error\n";
        times << "method,nx,ny,nz,batch,round,seconds_per_solve\n";
        // The pinned heFFTe no-reorder R2C path fails multi-batch equivalence.
        // Retain its single-slice measurement; never execute that path with batch > 1.
        const std::vector<Mode> modes{PublicReorder, PublicNoReorder, RowsTranspose, MergedC2C};
        for (auto shape : std::vector<std::pair<int, int>>{{6, 8}, {8, 6}}) {
            Problem p(shape.first, shape.second, 73);
            for (int batch : {1, 8, 32, 64}) {
                std::vector<std::unique_ptr<Candidate>> candidates;
                for (auto m : modes)
                    if (m != PublicNoReorder || batch == 1)
                        candidates.push_back(
                                std::make_unique<Candidate>(
                                        m, p.nx, p.ny, batch, p.gold.getGreensFnTr(), p.dx, p.dy));
                // Reuse instances across dense, neutral/sparse and empty sources: stale fields
                // fail.
                for (int source : {1, 2, 3}) {
                    p.fill(source);
                    p.reference(batch);
                    for (auto& c : candidates) {
                        c->solve(
                                p.density, p.field, p.occupied, 1. / Physics::epsilon_0,
                                &p.potential);
                        auto e = p.check();
                        checks << label(c->mode) << ',' << p.nx << ',' << p.ny << ',' << p.nz << ','
                               << batch << ',' << source << ',' << p.occupied.size() << ','
                               << e.first << ',' << e.second << '\n';
                    }
                }
            }
        }
        if (validateOnly) return;
        for (int nz : (quick ? std::vector<int>{256} : std::vector<int>{256, 512})) {
            Problem p(64, 64, nz);
            p.fill(0);
            for (int batch : (quick ? std::vector<int>{64} : std::vector<int>{1, 8, 32, 64})) {
                p.reference(batch);
                std::vector<std::unique_ptr<Candidate>> candidates;
                std::vector<Mode> selected;
                for (auto m : modes) {
                    if (m == PublicNoReorder && batch > 1) continue;
                    selected.push_back(m);
                    auto c = std::make_unique<Candidate>(
                            m, 64, 64, batch, p.gold.getGreensFnTr(), p.dx, p.dy);
                    c->solve(p.density, p.field, p.occupied, 1. / Physics::epsilon_0, &p.potential);
                    auto e = p.check();
                    checks << label(m) << ",64,64," << nz << ',' << batch << ",0,"
                           << p.occupied.size() << ',' << e.first << ',' << e.second << '\n';
                    candidates.push_back(std::move(c));
                }
                // Retained implementations have the same output and synchronized phase boundaries.
                auto operation = [&](int method) {
                    if (method == 0)
                        p.gold.solveBatched(
                                p.density, p.field, p.occupied, batch, 1. / Physics::epsilon_0,
                                &p.potential);
                    else
                        candidates[method - 1]->solve(
                                p.density, p.field, p.occupied, 1. / Physics::epsilon_0,
                                &p.potential);
                };
                for (int i = 0; i < 3; ++i)
                    for (int method = 0; method < int(candidates.size()) + 1; ++method)
                        operation(method);
                for (int round = 0; round < (quick ? 3 : 5); ++round)
                    for (int offset = 0; offset < int(candidates.size()) + 1; ++offset) {
                        int method     = (round + offset) % (int(candidates.size()) + 1);
                        double elapsed = measure(
                                [&] {
                                    operation(method);
                                },
                                quick ? 3 : 5);
                        times
                                << (method == 0 ? "current_solveBatched"
                                                : label(selected[method - 1]))
                                << ",64,64," << nz << ',' << batch << ',' << round << ',' << elapsed
                                << '\n';
                        times.flush();
                    }
                for (auto& c : candidates) {
                    c->solve(p.density, p.field, p.occupied, 1. / Physics::epsilon_0, &p.potential);
                    p.check();
                }
                std::cout << "COMPLETE nz=" << nz << " batch=" << batch << std::endl;
            }
        }
    }
}  // namespace
int main(int argc, char** argv) {
    bool validate = false, quick = false;
    for (int i = 1; i < argc; ++i) {
        validate |= std::string(argv[i]) == "--validate-only";
        quick |= std::string(argv[i]) == "--quick";
    }
    ippl::initialize(argc, argv);
    int result = 0;
    try {
        if (ippl::Comm->size() != 1)
            throw std::runtime_error("Local FFT microbenchmark requires one rank");
        run(validate, quick);
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << std::endl;
        result = 1;
    }
    ippl::finalize();
    return result;
}
