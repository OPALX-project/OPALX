template <class Level>
AmesosBottomSolver<Level>::AmesosBottomSolver(std::string solvertype)
    : solvertype_m(solvertype)
{ }


template <class Level>
AmesosBottomSolver<Level>::~AmesosBottomSolver() {
    solver_mp = Teuchos::null;
}


template <class Level>
void AmesosBottomSolver<Level>::solve(const Teuchos::RCP<mv_t>& x,
                                      const Teuchos::RCP<mv_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    solver_mp->solve(x.get(), b.get());
}


template <class Level>
void AmesosBottomSolver<Level>::setOperator(const Teuchos::RCP<matrix_t>& A,
                                            Level* level_p)
{
    try {
        solver_mp = Amesos2::create<matrix_t, mv_t>(solvertype_m, A);
    } catch(const std::invalid_argument& ex) {
        //TODO change to gmsg IPPL when built-in in OPAL
        std::cerr << ex.what() << std::endl;
    }
    
    solver_mp->symbolicFactorization();
    solver_mp->numericFactorization();
}


template <class Level>
std::size_t AmesosBottomSolver<Level>::getNumIters() {
    return 1;   // direct solvers do only one step
}
