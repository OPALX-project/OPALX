#ifndef OPAL_FIELD_SOLVER_H
#define OPAL_FIELD_SOLVER_H

#ifdef __CUDACC__
#pragma push_macro("__cpp_consteval")
#pragma push_macro("_NODISCARD")
#pragma push_macro("__builtin_LINE")

#define __cpp_consteval 201811L

#ifdef _NODISCARD
    #undef _NODISCARD
    #define _NODISCARD
#endif

#define consteval constexpr

#include <source_location>

#undef consteval
#pragma pop_macro("__cpp_consteval")
#pragma pop_macro("_NODISCARD")
#else
#include <source_location>
#endif

#include <memory>
#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"

template <typename T = double, unsigned Dim = 3>
using CGSolver_t = ippl::PoissonCG<Field<T, Dim>, Field_t<Dim>>;

using ippl::detail::ConditionalType, ippl::detail::VariantFromConditionalTypes;

template <typename T = double, unsigned Dim = 3>
using FFTSolver_t = ConditionalType<
    Dim == 2 || Dim == 3, ippl::FFTPeriodicPoissonSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using P3MSolver_t = ConditionalType<Dim == 3, ippl::P3MSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using OpenSolver_t =
    ConditionalType<Dim == 3, ippl::FFTOpenPoissonSolver<VField_t<T, Dim>, Field_t<Dim>>>;

// Define the FieldSolver class
template <typename T, unsigned Dim>
class FieldSolver : public ippl::FieldSolverBase<T, Dim> {
private:
    Field_t<Dim>* rho_m;
    VField_t<T, Dim>* E_m;
    Field_t<Dim>* phi_m;
    unsigned int call_counter_m;
public:
    FieldSolver(std::string solver, Field_t<Dim>* rho, VField_t<T, Dim>* E, Field<T, Dim>* phi)
        : ippl::FieldSolverBase<T, Dim>(solver), rho_m(rho), E_m(E), phi_m(phi), call_counter_m(0) {
        setPotentialBCs();
    }

    ~FieldSolver() {
    }

    void dumpScalField(std::string what);
    void dumpVectField(std::string what);

    Field_t<Dim>* getRho() {
        return rho_m;
    }
    void setRho(Field_t<Dim>* rho) {
        rho_m = rho;
    }

    VField_t<T, Dim>* getE() const {
        return E_m;
    }
    void setE(VField_t<T, Dim>* E) {
        E_m = E;
    }

    Field<T, Dim>* getPhi() const {
        return phi_m;
    }
    void setPhi(Field<T, Dim>* phi) {
        phi_m = phi;
    }

    void initOpenSolver();

    void initSolver() override ;

    void setPotentialBCs();

    void runSolver() override;
    
    template <typename Solver>
    void initSolverWithParams(const ippl::ParameterList& sp) {
std::cout << "1" << std::endl;
        this->getSolver().template emplace<Solver>();
        Solver& solver = std::get<Solver>(this->getSolver());
std::cout << "2" << std::endl;
        solver.mergeParameters(sp);
std::cout << "3" << std::endl;
        solver.setRhs(*rho_m);
std::cout << "4" << std::endl;
        if constexpr (std::is_same_v<Solver, CGSolver_t<T, Dim>>) {
            // The CG solver computes the potential directly and
            // uses this to get the electric field
std::cout << "5" << std::endl;
            solver.setLhs(*phi_m);
std::cout << "6" << std::endl;
            solver.setGradient(*E_m);
        } else if constexpr (std::is_same_v<Solver, OpenSolver_t<T, Dim>>) {
            // The periodic Poisson solver, Open boundaries solver,
            // and the P3M solver compute the electric field directly
std::cout << "52" << std::endl;
            solver.setLhs(*E_m);
std::cout << "62" << std::endl;
            solver.setGradFD();
        }
        call_counter_m = 0;
    }

    void initNullSolver();
    
    void initFFTSolver() {
    ippl::ParameterList sp;
    sp.add("output_type", FFTSolver_t<double, 3>::GRAD);
    sp.add("use_heffte_defaults", false);
    sp.add("use_pencils", true);
    sp.add("use_reorder", false);
    sp.add("use_gpu_aware", true);
    sp.add("comm", ippl::p2p_pl);
    sp.add("r2c_direction", 0);
    initSolverWithParams<FFTSolver_t<double, 3>>(sp);
    }
    
    void initCGSolver() { }

    void initP3MSolver() { }

};

// Explicit specialization declaration
template<>
void FieldSolver<double, 3>::initNullSolver();

#endif
