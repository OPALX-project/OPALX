#ifndef AMR_MULTI_GRID_H
#define AMR_MULTI_GRID_H

#include <vector>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include "AmrMultiGridCore.h"

#include "AmrMultiGridLevel.h"

#define AMR_MG_TIMER 1

class AmrMultiGrid {
    
public:
    typedef amr::matrix_t         matrix_t;
    typedef amr::vector_t         vector_t;
    typedef amr::multivector_t    mv_t;
    typedef amr::dmap_t           dmap_t;
    typedef amr::comm_t           comm_t;
    typedef amr::local_ordinal_t  lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::scalar_t         scalar_t;
    
    typedef AmrMultiGridLevel<matrix_t, vector_t> AmrMultiGridLevel_t;
    
    typedef AmrMultiGridLevel_t::AmrField_t     AmrField_t;
    typedef AmrMultiGridLevel_t::AmrGeometry_t  AmrGeometry_t;
    typedef AmrMultiGridLevel_t::AmrField_u     AmrField_u;
    typedef AmrMultiGridLevel_t::AmrField_s     AmrField_s;
    typedef AmrMultiGridLevel_t::AmrIntVect_t   AmrIntVect_t;
    typedef AmrMultiGridLevel_t::indices_t      indices_t;
    typedef AmrMultiGridLevel_t::coefficients_t coefficients_t;
    typedef AmrMultiGridLevel_t::umap_t         umap_t;

    typedef amrex::BoxArray boxarray_t;
    typedef amrex::Box box_t;
    typedef amrex::BaseFab<int> basefab_t;
    typedef amrex::FArrayBox farraybox_t;
            
    typedef AmrSmoother::Smoother Smoother;
    
    /// Supported interpolaters for prolongation operation
    enum Interpolater {
        TRILINEAR = 0,
        LAGRANGE,
        PIECEWISE_CONST
    };
    
    /// Supported bottom solvers
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
        // all Amesos2
#ifdef HAVE_AMESOS2_KLU2
        , KLU2
#endif
#ifdef HAVE_AMESOS2_SUPERLU
        , SUPERLU
#endif
#ifdef HAVE_AMESOS2_UMFPACK  
        , UMFPACK
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
        , PARDISO_MKL
#endif
#ifdef HAVE_AMESOS2_MUMPS
        , MUMPS
#endif
#ifdef HAVE_AMESOS2_LAPACK
        , LAPACK
#endif
        // add others ...
    };
    
    /// Supported physical boundaries
    enum Boundary {
        DIRICHLET = 0,
        OPEN,
        PERIODIC
    };
    
    /// Supported convergence criteria
    enum Norm {
        L1,
        L2,
        LINF
    };
    
public:
    
    /*!
     * Instantiation.
     * @param bc physical boundary condition
     * @param interp interpolater from coarse to fine grids without taking care of
     * coarse-fine interface
     * @param interface interpolater taking care of coarse-fine interface
     * @param solver for bottom level
     * @param preconditioner for bottom solver if available
     * @param smoother for error
     * @param norm convergence criteria
     */
    AmrMultiGrid(Boundary bc = Boundary::DIRICHLET,
                 Interpolater interp = Interpolater::TRILINEAR,
                 Interpolater interface = Interpolater::LAGRANGE,
                 BaseSolver solver = BaseSolver::CG,
                 Preconditioner precond = Preconditioner::NONE,
                 Smoother smoother = Smoother::JACOBI,
                 Norm norm = Norm::LINF
                );
    
    /*!
     * Compute the potential and the electric field
     * @param rho is the right-hand side
     * @param phi is the potential
     * @param efield is the electric field
     * @param geom specifies the geometry
     * @param lbase is the start level (currently only lbase = 0 supported)
     * @param lfine is the finest level
     * @param previous solution as initial guess
     */
    void solve(const amrex::Array<AmrField_u>& rho,
               amrex::Array<AmrField_u>& phi,
               amrex::Array<AmrField_u>& efield,
               const amrex::Array<AmrGeometry_t>& geom,
               int lbase, int lfine, bool previous = false);
    
    /*!
     * Specify the number of smoothing steps
     * @param nSweeps for each smoothing step
     */
    void setNumberOfSweeps(std::size_t nSweeps) {
        nSweeps_m = nSweeps;
    }
    
    
    /*!
     * Obtain some convergence info
     * @returns the number of iterations till convergence
     */
    std::size_t getNumIters() {
        return nIter_m;
    }
    
    
private:
    
    /*!
     * Instantiate boundary object
     */
    void initPhysicalBoundary_m();
    
    /*!
     * Instantiate all levels and set boundary conditions
     * @param rho is the charge density
     * @param geom is the geometry
     */
    void initLevels_m(const amrex::Array<AmrField_u>& rho,
                      const amrex::Array<AmrGeometry_t>& geom);
    
    /*!
     * Clear masks (required to build matrices) no longer needed.
     */
    void clearMasks_m();
    
    /*!
     * Reset potential to zero (currently)
     * @param phi is the potential
     * @param previous solution as initial guess
     */
    void initGuess_m(amrex::Array<AmrField_u>& phi, bool previous);
    
    /*!
     * Actual solve.
     */
    void iterate_m();
    
    /*!
     * Compute composite residual of a level
     * @param r is the residual to compute
     * @param b is the right-hand side
     * @param x is the left-hand side
     * @param level to solve for
     */
    void residual_m(const lo_t& level,
		    Teuchos::RCP<vector_t>& r,
                    const Teuchos::RCP<vector_t>& b,
                    const Teuchos::RCP<vector_t>& x);
    
    /*!
     * Recursive call.
     * @param level to relax
     */
    void relax_m(const lo_t& level);
    
    /*!
     * Compute the residual of a level without considering refined level.
     * @param result is computed
     * @param rhs is the right-hand side
     * @param crs_rhs is the coarse right-hand side for internal boundary
     * @param b is the left-hand side
     * @param level to solver for
     */
    void residual_no_fine_m(const lo_t& level,
			    Teuchos::RCP<vector_t>& result,
                            const Teuchos::RCP<vector_t>& rhs,
                            const Teuchos::RCP<vector_t>& crs_rhs,
                            const Teuchos::RCP<vector_t>& b);
                           
    /*!
     * @returns the maximum norm over all levels using the norm specified
     * by the user
     */
    scalar_t residualNorm_m();
    
    /*!
     * Vector norm computation.
     * @param x is the vector for which we compute the norm
     * @returns the evaluated norm of a level
     */
    scalar_t evalNorm_m(const Teuchos::RCP<const vector_t>& x);
    
    /*!
     * Initial convergence criteria values.
     * @param maxResidual maximum norm of residual over all levels
     * @param maxRho maximum norm of right-hand side over all levels
     */
    void initResidual_m(scalar_t& maxResidual, scalar_t& maxRho);
    
    /*!
     * @param efield to compute
     */
    void computeEfield_m(amrex::Array<AmrField_u>& efield);
    
    /*!
     * Build all matrices and vectors, i.e. AMReX to Trilinos
     * @param matrices if we need to build matrices as well or only vectors.
     */
    void setup_m(const amrex::Array<AmrField_u>& rho,
                 const amrex::Array<AmrField_u>& phi,
                 const bool& matrices = true);
    
    /*!
     * Build all matrices and vectors needed for single-level computation
     */
    void buildSingleLevel_m(const amrex::Array<AmrField_u>& rho,
			    const amrex::Array<AmrField_u>& phi,
			    const bool& matrices = true);
    
    /*!
     * Build all matrices and vectors needed for multi-level computation
     */
    void buildMultiLevel_m(const amrex::Array<AmrField_u>& rho,
			   const amrex::Array<AmrField_u>& phi,
			   const bool& matrices = true);

    /*!
     * Set matrix and vector pointer
     * @param level for which we fill matrix + vector
     * @param matrices if we need to set matrices as well or only vectors.
     */
    void open_m(const lo_t& level, const bool& matrices);
    
    /*!
     * Call fill complete
     * @param level for which we filled matrix
     * @param matrices if we set matrices as well.
     */
    void close_m(const lo_t& level, const bool& matrices);
    
    /*!
     * Build the Poisson matrix for a level assuming no finer level (i.e. the whole fine mesh
     * is taken into account).
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param level for which we build the Poisson matrix
     */
    void buildNoFinePoissonMatrix_m(const lo_t& level,
				    const go_t& gidx,
				    const AmrIntVect_t& iv,
                                    const basefab_t& mfab,
				    const scalar_t* invdx2);
    
    /*!
     * Build the Poisson matrix for a level that got refined (it does not take the covered
     * cells (covered by fine cells) into account). The finest level does not build such a matrix.
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param rfab is the mask between levels
     * @param level for which we build the special Poisson matrix
     */
    void buildCompositePoissonMatrix_m(const lo_t& level,
				       const go_t& gidx,
				       const AmrIntVect_t& iv,
                                       const basefab_t& mfab,
				       const basefab_t& rfab,
				       const basefab_t& cfab,
				       const scalar_t* invdx2);
    
    /*!
     * Build a matrix that averages down the data of the fine cells down to the
     * corresponding coarse cells. The base level does not build such a matrix.
     * \f[
     *      x^{(l)} = R\cdot x^{(l+1)}
     * \f]
     * @param iv is the current cell
     * @param rfab is the mask between levels
     * @param level for which to build restriction matrix
     */
    void buildRestrictionMatrix_m(const lo_t& level,
				  const go_t& gidx,
				  const AmrIntVect_t& iv,
				  D_DECL(const go_t& ii,
					 const go_t& jj,
					 const go_t& kk),
				  const basefab_t& rfab);
    
    /*!
     * Interpolate data from coarse cells to appropriate refined cells. The interpolation
     * scheme is allowed only to have local cells in the stencil, i.e.
     * 2D --> 4, 3D --> 8.
     * \f[
     *      x^{(l)} = I\cdot x^{(l-1)}
     * \f]
     * @param iv is the current cell
     * @param level for which to build the interpolation matrix. The finest level
     * does not build such a matrix.
     */
    void buildInterpolationMatrix_m(const lo_t& level,
				    const go_t& gidx,
				    const AmrIntVect_t& iv);
    
    /*!
     * The boundary values at the crse-fine-interface need to be taken into account properly.
     * This matrix is used to update the fine boundary values from the coarse values, i.e.
     * \f[
     *      x^{(l)} = B_{crse}\cdot x^{(l-1)}
     * \f]
     * Dirichlet boundary condition
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param cells all fine cells that are at the crse-fine interface
     * @param level the base level is omitted
     */
    void buildCrseBoundaryMatrix_m(const lo_t& level,
				   const go_t& gidx,
				   const AmrIntVect_t& iv,
				   const basefab_t& mfab,
				   const basefab_t& cfab,
				   const scalar_t* invdx2);
    
    /*!
     * The boundary values at the crse-fine-interface need to be taken into account properly.
     * This matrix is used to update the coarse boundary values from fine values, i.e.
     * \f[
     *      x^{(l)} = B_{fine}\cdot x^{(l+1)}
     * \f]
     * Dirichlet boundary condition. Flux matching.
     * @param iv is the current cell
     * @param cells all coarse cells that are at the crse-fine interface but are
     * not refined
     * @param crse_fine_ba coarse cells that got refined
     * @param level the finest level is omitted
     */
    void buildFineBoundaryMatrix_m(const lo_t& level,
				   const go_t& gidx,
				   const AmrIntVect_t& iv,
                                   const basefab_t& mfab,
				   const basefab_t& rfab,
				   const basefab_t& cfab);
    
    /*!
     * Copy data from AMReX to Trilinos
     * @param rho is the charge density
     * @param level for which to copy
     */
    void buildDensityVector_m(const lo_t& level,
			      const AmrField_t& rho);
    
    /*!
     * Copy data from AMReX to Trilinos
     * @param phi is the potential
     * @param level for which to copy
     */
    void buildPotentialVector_m(const lo_t& level,
				const AmrField_t& phi);
    
    /*!
     * Gradient matrix is used to compute the electric field
     * @param iv is the current cell
     * @param mfab is the mask (internal cell, boundary cell, ...) of that level
     * @param level for which to compute
     */
    void buildGradientMatrix_m(const lo_t& level,
			       const go_t& gidx,
			       const AmrIntVect_t& iv,
                               const basefab_t& mfab,
			       const scalar_t* invdx);
    
    /*!
     * Data transfer from AMReX to Trilinos.
     * @param mf is the multifab of a level
     * @param comp component to copy
     * @param mv is the vector to be filled
     * @param level where to perform
     */
    void amrex2trilinos_m(const lo_t& level,
			  const lo_t& comp,
			  const AmrField_t& mf,
			  Teuchos::RCP<vector_t>& mv);
    
    /*!
     * Data transfer from Trilinos to AMReX.
     * @param mf is the multifab to be filled
     * @param comp component to copy
     * @param mv is the corresponding Trilinos vector
     */
    void trilinos2amrex_m(const lo_t& comp,
			  AmrField_t& mf,
			  const Teuchos::RCP<vector_t>& mv);
    
    /*!
     * Some indices might occur several times due to boundary conditions, etc. We
     * avoid this by filling a map and then copying the data to a vector for filling
     * the matrices. The map gets cleared inside the function.
     * @param indices in matrix
     * @param values are the coefficients
     */
    void map2vector_m(umap_t& map, indices_t& indices,
		      coefficients_t& values);
    
    /*!
     * Perform one smoothing step
     * @param e error to update (left-hand side)
     * @param r residual (right-hand side)
     * @param level on which to relax
     */
    void smooth_m(const lo_t& level,
		  Teuchos::RCP<vector_t>& e,
                  Teuchos::RCP<vector_t>& r);
    
    /*!
     * Restrict coarse level residual based on fine level residual
     * @param level to restrict
     */
    void restrict_m(const lo_t& level);
    
    /*!
     * Update error of fine level based on error of coarse level
     * @param level to update
     */
    void prolongate_m(const lo_t& level);
    
    /*!
     * Average data from fine level to coarse
     * @param level finest level is omitted
     */
    void averageDown_m(const lo_t& level);
    
    /*!
     * Instantiate interpolater
     * @param interp interpolater type
     */
    void initInterpolater_m(const Interpolater& interp);
    
    /*!
     * Instantiate interface interpolater
     * @param interface handler
     */
    void initCrseFineInterp_m(const Interpolater& interface);
    
    /*!
     * Instantiate a bottom solver
     * @param solver type
     * @param precond preconditioner
     */
    void initBaseSolver_m(const BaseSolver& solver,
                          const Preconditioner& precond);
    
private:
    Teuchos::RCP<comm_t> comm_mp;       ///< communicator
    Teuchos::RCP<amr::node_t> node_mp;  ///< kokkos node
    
    /// interpolater without coarse-fine interface
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interp_mp;
    
    /// interpolater for coarse-fine interface
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interface_mp;
    
    std::size_t nIter_m;            ///< number of iterations till convergence
    std::size_t nSweeps_m;          ///< number of smoothing iterations
    Smoother smootherType_m;        ///< type of smoother
    
    /// container for levels
    std::vector<std::unique_ptr<AmrMultiGridLevel_t > > mglevel_m;
    
    /// bottom solver
    std::shared_ptr<BottomSolver<Teuchos::RCP<matrix_t>, Teuchos::RCP<mv_t> > > solver_mp;
    
    /// error smoother
    std::vector<std::shared_ptr<AmrSmoother> > smoother_m;
    
    int lbase_m;            ///< base level (currently only 0 supported)
    int lfine_m;            ///< fineste level
    int nlevel_m;           ///< number of levelss
    
    Boundary bcType_m;          ///< physical boundary type
    std::shared_ptr<AmrBoundary<AmrMultiGridLevel_t> > bc_m;
    
    Norm norm_m;            ///< norm for convergence criteria (l1, l2, linf)
    
#if AMR_MG_TIMER
    IpplTimings::TimerRef buildTimer_m;         ///< timer for matrix and vector construction
    IpplTimings::TimerRef restrictTimer_m;      ///< timer for restriction operation
    IpplTimings::TimerRef smoothTimer_m;        ///< timer for all smoothing steps
    IpplTimings::TimerRef interpTimer_m;        ///< prolongation timer
    IpplTimings::TimerRef residnofineTimer_m;   ///< timer for no-fine residual computation
    IpplTimings::TimerRef bottomTimer_m;        ///< bottom solver timer
#endif

    IpplTimings::TimerRef bopen_m;
    IpplTimings::TimerRef bclose_m;
    IpplTimings::TimerRef bclear_m;
    IpplTimings::TimerRef bRestict_m;
    IpplTimings::TimerRef bInterp_m;
    IpplTimings::TimerRef bCompo_m;
    IpplTimings::TimerRef bPoiss_m;
    IpplTimings::TimerRef bBf_m;
    IpplTimings::TimerRef bBc_m;
    IpplTimings::TimerRef bG_m;
    IpplTimings::TimerRef bSmoother_m;
};

#endif
