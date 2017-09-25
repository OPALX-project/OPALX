#ifndef AMR_MULTI_GRID_H
#define AMR_MULTI_GRID_H

#include <vector>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>


#include "AmrMultiGridLevel.h"

#include "AmrTrilinearInterpolater.h"


#include "TrilinosSolver.h"


class AmrMultiGrid {
    
public:
    typedef TrilinosSolver::matrix_t matrix_t;
    typedef TrilinosSolver::vector_t vector_t;
    
    typedef AmrMultiGridLevel<matrix_t, vector_t> AmrMultiGridLevel_t;
    
    typedef AmrMultiGridLevel_t::AmrField_t AmrField_t;
    typedef AmrMultiGridLevel_t::AmrGeometry_t AmrGeometry_t;
    typedef AmrMultiGridLevel_t::AmrField_u AmrField_u;
    typedef AmrMultiGridLevel_t::AmrField_s AmrField_s;
    typedef AmrMultiGridLevel_t::AmrIntVect_t AmrIntVect_t;
    typedef AmrMultiGridLevel_t::indices_t indices_t;
    typedef AmrMultiGridLevel_t::coefficients_t coefficients_t;
    
    typedef std::map<AmrIntVect_t,
                     std::list<std::pair<int, int> >
                     > map_t;
    
    enum Interpolater {
        TRILINEAR = 0 //,
//         TRICUBIC  = 1,
    };
    
public:
    
    AmrMultiGrid(Interpolater interp = Interpolater::TRILINEAR);
    
    void solve(const amrex::Array<AmrField_u>& rho,
               amrex::Array<AmrField_u>& phi,
               amrex::Array<AmrField_u>& efield,
               const amrex::Array<AmrGeometry_t>& geom,
               int lbase, int lfine, bool previous = false);
    
    
private:
    
    void residual_m(Teuchos::RCP<vector_t>& r,
                    const Teuchos::RCP<vector_t>& b,
                    const Teuchos::RCP<matrix_t>& A,
                    const Teuchos::RCP<vector_t>& x,
                    bool temporary=false);
    
    void relax_m(int level);
    
    void residual_no_fine_m(Teuchos::RCP<vector_t>& result,
                            const Teuchos::RCP<vector_t>& rhs,
                            const Teuchos::RCP<vector_t>& crs_rhs,
                            const Teuchos::RCP<vector_t>& b,
                            int level, bool plus = false);
                           
    
    
    double l2error_m();    
    
    /*!
     * Build the Poisson matrix for a level assuming no finer level (i.e. the whole fine mesh
     * is taken into account).
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * FIXME The number of entries due to the boundary conditions has to be fixed.
     * @param level for which we build the Poisson matrix
     */
    void buildPoissonMatrix_m(int level);
    
    /*!
     * Build the Poisson matrix for a level that got refined (it does not take the covered
     * cells (covered by fine cells) into account). The finest level does not build such a matrix.
     * It takes care of physical boundaries (i.e. mesh boundary).
     * Internal boundaries (i.e. boundaries due to crse-fine interfaces) are treated by the
     * boundary matrix.
     * FIXME The number of entries due to the boundary conditions has to be fixed.
     * @param level for which we build the special Poisson matrix
     */
    void buildSpecialPoissonMatrix_m(int level);
    
    /*!
     * Build a matrix that averages down the data of the fine cells down to the
     * corresponding coarse cells. The base level does not build such a matrix.
     * \f[
     *      x^{(l-1)} = R\cdot x^{(l)}
     * \f]
     * @param level for which to build restriction matrix
     */
    void buildRestrictionMatrix_m(int level);
    
    /*!
     * Interpolate data from coarse cells to appropriate refined cells. The interpolation
     * scheme is allowed only to have local cells in the stencil, i.e.
     * 2D --> 4, 3D --> 8.
     * \f[
     *      x^{(l+1)} = I\cdot x^{(l)}
     * \f]
     * @param level for which to build the interpolation matrix. The finest level
     * does not build such a matrix.
     * FIXME The number of entries due to the interpolation stencil has to be fixed.
     */
    void buildInterpolationMatrix_m(int level);
    
    /*!
     * The boundary values at the crse-fine-interface need to be taken into account properly.
     * This matrix is used to update the fine boundary values from the coarse values, i.e.
     * \f[
     *      x^{(l+1)} = B_{crse}\cdot x^{(l)}
     * \f]
     * Dirichlet boundary condition
     * FIXME The number of entries due to the interpolation stencil has to be fixed.
     * @param level the base level is omitted
     */
    void buildCrseBoundaryMatrix_m(int level);
    
    /*!
     * 
     * @param level the finest level is omitted
     */
    void buildFineBoundaryMatrix_m(int level);
    
    void buildDensityVector_m(const AmrField_t& rho, int level);
    
    void buildPotentialVector_m(const AmrField_t& phi, int level);
    
    /*!
     * The smoother matrix is used in the relaxation step. The base level
     * does not require a smoother matrix.
     * @param level for which to build matrix
     */
    void buildSmootherMatrix_m(int level);
    
    /*!
     * Data transfer from AMReX to Trilinos.
     * @param mf is the multifab of a level
     * @param mv is the vector to be filled
     * @param level where to perform
     */
    void amrex2trilinos_m(const AmrField_t& mf, Teuchos::RCP<vector_t>& mv, int level);
    
    /*!
     * Data transfer from Trilinos to AMReX.
     * @param mf is the multifab to be filled
     * @param mv is the corresponding Trilinos vector
     */
    void trilinos2amrex_m(AmrField_t& mf, const Teuchos::RCP<vector_t>& mv);
    
    /*!
     * Some indices might occur several times due to boundary conditions, etc.
     * This function removes all duplicates and sums the coefficients up.
     * @param indices in matrix
     * @param values are the coefficients
     */
    void unique_m(indices_t& indices,
                  coefficients_t& values);
    
    
    void checkCrseBoundary_m(Teuchos::RCP<matrix_t>& B,
                             int level,
                             const amrex::BaseFab<int>& mfab,
                             const AmrIntVect_t& lo,
                             const AmrIntVect_t& hi);
    
private:
    Epetra_MpiComm epetra_comm_m;
    
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interp_mp;
    
    int nIter_m;
    
    
    std::vector<std::unique_ptr<AmrMultiGridLevel_t > > mglevel_m;
    
//     std::unique_ptr<LinearSolver<Teuchos::RCP<matrix_t>, Teuchos::RCP<vector_t> > > solver_mp;
    TrilinosSolver solver_m;
    
    AmrIntVect_t rr_m;      //TODO move to level
    
    int lbase_m;
    int lfine_m;
    
};

#endif
