#ifndef AMR_MULTI_GRID_H
#define AMR_MULTI_GRID_H

#include <vector>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>

#include "AmrMultiGridCore.h"

#include "AmrMultiGridLevel.h"

#define AMR_MG_TIMER 1

class AmrMultiGrid {
    
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::dmap_t   dmap_t;
    typedef amr::comm_t   comm_t;
    
    typedef AmrMultiGridLevel<matrix_t, vector_t> AmrMultiGridLevel_t;
    
    typedef AmrMultiGridLevel_t::AmrField_t AmrField_t;
    typedef AmrMultiGridLevel_t::AmrGeometry_t AmrGeometry_t;
    typedef AmrMultiGridLevel_t::AmrField_u AmrField_u;
    typedef AmrMultiGridLevel_t::AmrField_s AmrField_s;
    typedef AmrMultiGridLevel_t::AmrIntVect_t AmrIntVect_t;
    typedef AmrMultiGridLevel_t::indices_t indices_t;
    typedef AmrMultiGridLevel_t::coefficients_t coefficients_t;
    
    typedef amrex::BoxArray boxarray_t;
    typedef amrex::Box box_t;
    typedef amrex::BaseFab<int> basefab_t;
    typedef amrex::FArrayBox farraybox_t;
    
    typedef std::map<AmrIntVect_t,
                     std::list<std::pair<int, int> >
                     > map_t;
    
    typedef AmrSmoother::Smoother Smoother;
    
    enum Interpolater {
        TRILINEAR = 0,
        LAGRANGE,
        PIECEWISE_CONST
    };
    
    enum BaseSolver {
        // all Belos
        BICGSTAB,
        MINRES,
        PCPG,
        CG,
        GMRES,
        STOCHASTIC_CG,
        RECYCLING_CG,
        RECYCLING_GMRES
        // .. add others
    };
    
    enum Boundary {
        DIRICHLET = 0,
        OPEN
    };
    
public:
    
    AmrMultiGrid(Boundary bc = Boundary::DIRICHLET,
                 Interpolater interp = Interpolater::TRILINEAR,
                 Interpolater interface = Interpolater::LAGRANGE,
                 BaseSolver solver = BaseSolver::CG,
                 Preconditioner precond = Preconditioner::NONE,
                 Smoother smoother = Smoother::JACOBI);
    
    void solve(const amrex::Array<AmrField_u>& rho,
               amrex::Array<AmrField_u>& phi,
               amrex::Array<AmrField_u>& efield,
               const amrex::Array<AmrGeometry_t>& geom,
               int lbase, int lfine, bool previous = false);
    
    void setNumberOfSweeps(std::size_t nSweeps) {
        nSweeps_m = nSweeps;
    }
    
    
    std::size_t getNumIters() {
        return nIter_m;
    }
    
    
private:
    
    void residual_m(Teuchos::RCP<vector_t>& r,
                    const Teuchos::RCP<vector_t>& b,
//                     const Teuchos::RCP<matrix_t>& A,
                    const Teuchos::RCP<vector_t>& x,
                    int level);
    
    void relax_m(int level);
    
    void residual_no_fine_m(Teuchos::RCP<vector_t>& result,
                            const Teuchos::RCP<vector_t>& rhs,
                            const Teuchos::RCP<vector_t>& crs_rhs,
                            const Teuchos::RCP<vector_t>& b,
                            int level);
                           
    
    
    double l2error_m();
    
    double lInfError_m();
    
    /*!
     * Build all matrices and vectors, i.e. AMReX to Trilinos
     * @param matrices if we need to build matrices as well or only vectors.
     */
    void setup_m(const amrex::Array<AmrField_u>& rho,
                 const amrex::Array<AmrField_u>& phi,
                 bool matrices = true);
    
    /*!
     * Set matrix and vector pointer
     * @param level for which we fill matrix + vector
     * @param matrices if we need to set matrices as well or only vectors.
     */
    void open_m(int level, bool matrices);
    
    /*!
     * Call fill complete
     * @param level for which we filled matrix
     * @param matrices if we set matrices as well.
     */
    void close_m(int level, bool matrices);
    
    /*!
     * Build the Poisson matrix for a level assuming no finer level (i.e. the whole fine mesh
     * is taken into account).
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * @param level for which we build the Poisson matrix
     */
    void buildNoFinePoissonMatrix_m(const AmrIntVect_t& iv,
                                    const basefab_t& mfab, int level);
    
    /*!
     * Build the Poisson matrix for a level that got refined (it does not take the covered
     * cells (covered by fine cells) into account). The finest level does not build such a matrix.
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * @param level for which we build the special Poisson matrix
     */
    void buildCompositePoissonMatrix_m(const AmrIntVect_t& iv,
                                       const basefab_t& mfab,
                                       const boxarray_t& ba,
                                       int level);
    
    /*!
     * Build a matrix that averages down the data of the fine cells down to the
     * corresponding coarse cells. The base level does not build such a matrix.
     * \f[
     *      x^{(l-1)} = R\cdot x^{(l)}
     * \f]
     * @param level for which to build restriction matrix
     */
    void buildRestrictionMatrix_m(const AmrIntVect_t& iv,
                                  const boxarray_t& crse_fine_ba,
                                  D_DECL(int& ii, int& jj, int& kk),
                                  int level);
    
    /*!
     * Interpolate data from coarse cells to appropriate refined cells. The interpolation
     * scheme is allowed only to have local cells in the stencil, i.e.
     * 2D --> 4, 3D --> 8.
     * \f[
     *      x^{(l+1)} = I\cdot x^{(l)}
     * \f]
     * @param level for which to build the interpolation matrix. The finest level
     * does not build such a matrix.
     */
    void buildInterpolationMatrix_m(const AmrIntVect_t& iv,
                                    int level);
    
    /*!
     * The boundary values at the crse-fine-interface need to be taken into account properly.
     * This matrix is used to update the fine boundary values from the coarse values, i.e.
     * \f[
     *      x^{(l+1)} = B_{crse}\cdot x^{(l)}
     * \f]
     * Dirichlet boundary condition
     * @param level the base level is omitted
     */
    void buildCrseBoundaryMatrix_m(const AmrIntVect_t& iv,
                                   map_t& cells, int level);
    
    void fillCrseBoundaryMatrix_m(map_t& cells,
                                  const boxarray_t& crse_fine_ba,
                                  int level);
    
    /*!
     * 
     * @param level the finest level is omitted
     */
    void buildFineBoundaryMatrix_m(const AmrIntVect_t& iv,
                                   map_t& cells,
                                   const boxarray_t& crse_fine_ba,
                                   int level);
    
    void fillFineBoundaryMatrix_m(map_t& cells,
                                  const boxarray_t& crse_fine_ba,
                                  int level);
    
    void buildDensityVector_m(const AmrField_t& rho,
                              int level);
    
    void buildPotentialVector_m(const AmrField_t& phi,
                                int level);
    
//     /*!
//      * The smoother matrix is used in the relaxation step. The base level
//      * does not require a smoother matrix.
//      * @param level for which to build matrix
//      */
//     void buildSmootherMatrix_m(int level);
    
    
    /*!
     * Gradient matrix is used to compute the electric field
     */
    void buildGradientMatrix_m(const AmrIntVect_t& iv,
                               const basefab_t& mfab, int level);
    
    /*!
     * Data transfer from AMReX to Trilinos.
     * @param mf is the multifab of a level
     * @param comp component to copy
     * @param mv is the vector to be filled
     * @param level where to perform
     */
    void amrex2trilinos_m(const AmrField_t& mf, int comp, Teuchos::RCP<vector_t>& mv, int level);
    
    /*!
     * Data transfer from Trilinos to AMReX.
     * @param mf is the multifab to be filled
     * @param comp component to copy
     * @param mv is the corresponding Trilinos vector
     */
    void trilinos2amrex_m(AmrField_t& mf, int comp, const Teuchos::RCP<vector_t>& mv);
    
    /*!
     * Some indices might occur several times due to boundary conditions, etc.
     * This function removes all duplicates and sums the coefficients up.
     * @param indices in matrix
     * @param values are the coefficients
     */
    void unique_m(indices_t& indices,
                  coefficients_t& values);
    
    
    void smooth_m(Teuchos::RCP<vector_t>& e,
                  Teuchos::RCP<vector_t>& r,
                  int level);
    
    void restrict_m(int level);
    
    void prolongate_m(int level);
    
    void averageDown_m(int level);
    
    
    void initInterpolater_m(const Interpolater& interp);
    
    void initCrseFineInterp_m(const Interpolater& interfaces);
    
    void initBaseSolver_m(const BaseSolver& solver,
                          const Preconditioner& precond);
    
    void writeYt_m(const amrex::Array<AmrField_u>& rho,
                   amrex::Array<AmrField_u>& phi,
                   amrex::Array<AmrField_u>& efield,
                   const amrex::Array<AmrGeometry_t>& geom);
    
private:
    Teuchos::RCP<comm_t> comm_mp;
    Teuchos::RCP<amr::node_t> node_mp;
    
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interp_mp;
    
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interface_mp;
    
    std::size_t nIter_m;
    std::size_t nSweeps_m;
    Smoother smootherType_m;
    
    std::vector<std::unique_ptr<AmrMultiGridLevel_t > > mglevel_m;
    
    std::shared_ptr<BottomSolver<Teuchos::RCP<matrix_t>, Teuchos::RCP<vector_t> > > solver_mp;
    
    std::vector<std::unique_ptr<AmrSmoother> > smoother_m;
    
    int lbase_m;
    int lfine_m;
    
    Boundary bc_m;
    
#if AMR_MG_TIMER
    IpplTimings::TimerRef buildTimer_m;
    IpplTimings::TimerRef restrictTimer_m;
    IpplTimings::TimerRef smoothTimer_m;
    IpplTimings::TimerRef interpTimer_m;
    IpplTimings::TimerRef residnofineTimer_m;
    IpplTimings::TimerRef bottomTimer_m;
#endif
    
};

#endif
