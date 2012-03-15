#ifdef HAVE_ML_SOLVER

#include "MGPoissonSolver.h"
#include "Physics/Physics.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include <algorithm>

#include "EpetraExt_RowMatrixOut.h"
#include <Epetra_Import.h>

#include "TaperDomain.h"

MGPoissonSolver::MGPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::vector<BoundaryGeometry*> geometries, std::string itsolver, std::string interpl, double tol, int maxiters, std::string precmode):
    layout_m(fl),
    mesh_m(mesh),
    LHS(0),
    geometries_m(geometries),
    tol_m(tol),
    maxiters_m(maxiters),
    Comm(MPI_COMM_WORLD)
{ 
    domain_m = layout_m->getDomain();
    e_dim_tag decomp[3]; 

    for(int i=0; i<3; i++) {
        decomp[i] = layout_m->getRequestedDistribution(i);
        hr_m[i] = mesh_m->get_meshSpacing(i);
        orig_nr_m[i] = domain_m[i].length();
    }

    //TODO: incorporate changes from stand-alone solver: x,y parallel
    if(decomp[0] == 1 || decomp[1] == 1) {
        throw OpalException("MGPoissonSolver","MGPoissonSolver only works parallel in z-direction");
    }

    if(itsolver == "CG") itsolver_m = CG;
    else if(itsolver == "BiCGStab") itsolver_m = BiCGStab;
    else if(itsolver == "GMRES") itsolver_m = GMRES;
    else itsolver_m = NONE;
    
    precmode_m = STD_PREC;
    if(precmode == "STD") precmode_m = STD_PREC;
    else if(precmode == "HIERARCHY") precmode_m = REUSE_HIERARCHY;
    else if(precmode == "REUSE") precmode_m = REUSE_PREC;


    //TODO: find CURRENT geometry
    //IFF: at the moment we only have 1 geometry 
    currentGeometry = geometries_m[0];

    if(currentGeometry->getFilename() == "") {
        //bp = new RectangularDomain(currentGeometry->getA(), currentGeometry->getB(), nr_m, hr_m);
        bp = new EllipticalDomain(currentGeometry->getA(), currentGeometry->getB(), orig_nr_m, hr_m, interpl);
    } /*else {
        bp = new ArbitraryDomain(currentGeometry->getFilename(), nr_m, hr_m, interpl);
    }*/

    //IFF: TAPER TEST
    //bp = new TaperDomain(0.02, 0.02, orig_nr_m, hr_m, interpl, 1.4);

    Map = 0;
    A = 0;
    LHS = 0;
    MLPrec = 0;
    
    // set defaults
    ML_Epetra::SetDefaults("SA", MLList);
    MLList.set("max levels", 5);
    MLList.set("increasing or decreasing", "decreasing");
    MLList.set("prec type", "MGV");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","Chebyshev");
    MLList.set("smoother: sweeps",3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type", "Amesos-KLU");
    //MLList.set("ML output", 10);
    //MLList.set("coarse: max size", 1024);
    // try to optimize mem for xt3
    //MLList.set("low memory usage", true);

    //all timers used here
    FunctionTimer1_m = IpplTimings::getTimer("CreateStencil");
    FunctionTimer2_m = IpplTimings::getTimer("PrepareRHS");
    FunctionTimer3_m = IpplTimings::getTimer("CG");
    FunctionTimer4_m = IpplTimings::getTimer("LHStoRho");
    FunctionTimer5_m = IpplTimings::getTimer("recalcMap");
    FunctionTimer6_m = IpplTimings::getTimer("ML");
}


MGPoissonSolver::MGPoissonSolver(PartBunch &beam):
    layout_m(&beam.getFieldLayout()),
    mesh_m(&beam.getMesh()),
    itsBunch_m(&beam),
    LHS(0),
    Comm(MPI_COMM_WORLD)
{
    domain_m = layout_m->getDomain();
    e_dim_tag decomp[3]; 

    for(int i=0; i<3; i++) {
        decomp[i] = layout_m->getRequestedDistribution(i);
        hr_m[i] = mesh_m->get_meshSpacing(i);
        orig_nr_m[i] = domain_m[i].length();
    }

    if(decomp[0] == 1 || decomp[1] == 1) {
        throw OpalException("MGPoissonSolver","MGPoissonSolver only works parallel in z-direction");
    }

    /*
    if(itsolver == "CG") itsolver_m = CG;
    else if(itsolver == "BiCGStab") itsolver_m = BiCGStab;
    else if(itsolver == "GMRES") itsolver_m = GMRES;
    else itsolver_m = NONE;

    //TODO: find CURRENT geometry
    //IFF: at the moment we only have 1 geometry 
    currentGeometry = geometries_m[0];

    if(currentGeometry->getFilename() == "") {
        bp = new EllipticalDomain(currentGeometry->getA(), currentGeometry->getB(), nr_m, hr_m, interpl);
    }*/ /*else {
        bp = new ArbitraryDomain(currentGeometry->getFilename(), nr_m, hr_m, interpl);
    }*/

    Map = 0;
    A = 0;
    LHS = 0;
    MLPrec = 0;
    
    // set defaults
    ML_Epetra::SetDefaults("SA", MLList);
    MLList.set("max levels", 5);
    MLList.set("increasing or decreasing", "decreasing");
    MLList.set("prec type", "MGV");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","Chebyshev");
    MLList.set("smoother: sweeps",3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type", "Amesos-KLU");
    //MLList.set("ML output", 10);
    //MLList.set("coarse: max size", 1024);
    // try to optimize mem for xt3
    //MLList.set("low memory usage", true);

    //all timers used here
    FunctionTimer1_m = IpplTimings::getTimer("Create Stencil");
    FunctionTimer2_m = IpplTimings::getTimer("Prepare RHS");
    FunctionTimer3_m = IpplTimings::getTimer("CG");
    FunctionTimer4_m = IpplTimings::getTimer("LHS to rho");
    FunctionTimer5_m = IpplTimings::getTimer("recalc Map");
    FunctionTimer6_m = IpplTimings::getTimer("ML");
}

////////////////////////////////////////////////////////////////////////////
// destructor
MGPoissonSolver::~MGPoissonSolver()
{
    //cleanup
    if(precmode_m != STD_PREC)
        delete MLPrec;
    delete LHS;
    delete A;
    delete Map;
}

////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation
// TODO: use matrix stencil in computation directly (no Epetra, define operators
// on IPPL GRID)
void MGPoissonSolver::computePotential(Field_t &rho, Vector_t hr)
{
	//TODO: needed?
    nr_m[0] = orig_nr_m[0];
    nr_m[1] = orig_nr_m[1];
    nr_m[2] = orig_nr_m[2];
    bp->setNr(nr_m);

    //FIXME: we need to do this only when geometry changes!: DO THIS IN APPR.
    //BC CLASS
    bp->Compute(hr);

    IpplTimings::startTimer(FunctionTimer5_m);
    NDIndex<3> localidx =  layout_m->getLocalNDIndex();
    //FIXME: when do we need to recreate the map! when geometry or parallel decomposition changes
    if(Map != 0 /* || hasParallelDecompositionChanged()*/)
        delete Map;
    Map = IPPLToMap3D(localidx);

    //we also have to redistribute LHS
    if(LHS == 0) {
        LHS = new Epetra_Vector(*Map);
        LHS->PutScalar(0.0);
    } else { //if(hasParallelDecompositionChanged()) { //FIXME: only when parallel decompisition has changed
        Epetra_Vector *tmp = new Epetra_Vector(*Map);
        Epetra_Import importer(*Map, LHS->Map());
        tmp->Import(*LHS, importer, Add);
        delete LHS;
        LHS = tmp;
    }
    IpplTimings::stopTimer(FunctionTimer5_m);

    IpplTimings::startTimer(FunctionTimer2_m);
    // get charge densities from ippl field -> Epetra vector (RHS)
    Epetra_Vector RHS(*Map);
    register int idx = 0;
    RHS.PutScalar(0.0);
    for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {

                if(bp->isInside(x,y,z)) {
                    RHS.Values()[idx++] = rho[x][y][z].get();
                }
            }
        }
    }
    IpplTimings::stopTimer(FunctionTimer2_m);

    IpplTimings::startTimer(FunctionTimer1_m);
    //REUSE HIERARCHY
    //FIXME: if GEOMETRY CHANGES A CHANGES!
    //if(A == 0 /* || geometryHasChanged()*/ ) 
    //   A = new Epetra_CrsMatrix(Copy, *Map,  7);
    
    //STD VERSION: BUILD PREC IN EVERY STEP
    //if(A != 0)
    //    delete A;
    //A = new Epetra_CrsMatrix(Copy, *Map,  7);
    //A = Stencil3DGeometry(hr, RHS);

    
    switch(precmode_m) {
        case STD_PREC:
            //if(A != 0)
            delete A;
            A = new Epetra_CrsMatrix(Copy, *Map,  7, true);
            break;
   
        case REUSE_PREC:
        case REUSE_HIERARCHY: // BUT RECOMPUTE PRECONDITIONER
            //REUSE HIERARCHY
            //FIXME: if GEOMETRY CHANGES A CHANGES!
            if(A == 0 /* || geometryHasChanged()*/ ) 
                A = new Epetra_CrsMatrix(Copy, *Map,  7, true);
            break;
    }

    // build discretization matrix
    Stencil3DGeometry(hr, A, RHS);
    IpplTimings::stopTimer(FunctionTimer1_m);

    // debug output: discretization stencil
#ifdef DBG_STENCIL 
    EpetraExt::RowMatrixToMatlabFile("DiscrStencil.dat", *A);
#endif

    IpplTimings::startTimer(FunctionTimer6_m);

    //STD VERSION
    //ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
   
    //REUSE PREC
    //if(MLPrec == 0)
    //    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
    
    //REUSE HIERARCHY BUT RECOMPUTE PRECONDITIONER
    //if(MLPrec == 0)
       //MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
    //else
       //MLPrec->ReComputePreconditioner();
    
    switch(precmode_m) {
        case REUSE_PREC:
            if(MLPrec == 0)
                MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
            break;
    
        case REUSE_HIERARCHY: // BUT RECOMPUTE PRECONDITIONER
            if(MLPrec == 0)
                MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
            else
                MLPrec->ReComputePreconditioner();
            break;
        
        case STD_PREC:
            MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);
            break;
    }

    IpplTimings::stopTimer(FunctionTimer6_m);

    // setup preconditioned iterativ solver
    // use old LHS solution as initial guess
    AztecOO solver(A, LHS, &RHS);
    solver.SetPrecOperator(MLPrec);

    //solver.SetAllAztecOptions(options);
    //solver.SetAllAztecParams(params);
    //solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
    //solver.SetAztecOption(AZ_precond, AZ_none);
    //solver.SetAztecParam(AZ_kspace, maxiters_m);

    switch(itsolver_m) {
        case CG: solver.SetAztecOption(AZ_solver, AZ_cg); break;
        case BiCGStab: solver.SetAztecOption(AZ_solver, AZ_bicgstab); break;
        case GMRES: solver.SetAztecOption(AZ_solver, AZ_gmres); break;
        default: throw OpalException("MGPoissonSolver","No valid iterative solver attached"); break;
    }

    solver.SetAztecOption(AZ_conv, AZ_rhs);
    solver.SetAztecOption(AZ_output, AZ_all);
    //solver.SetAztecOption(AZ_output, AZ_none);
    IpplTimings::startTimer(FunctionTimer3_m);
    solver.Iterate(maxiters_m, tol_m);
    IpplTimings::stopTimer(FunctionTimer3_m);

    //if( Comm.MyPID()==0 ) 
    //    cout << "\t\t||b-Ax||_2 = " << solver.TrueResidual() << endl;

    //now transfer solution back to IPPL grid
    IpplTimings::startTimer(FunctionTimer4_m);
    idx = 0;
    rho = 0.0;
    for(int x = localidx[0].first(); x <= localidx[0].last(); x++) {
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++) {
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {

                NDIndex<3> l(Index(x,x),Index(y,y),Index(z,z));
                if(bp->isInside(x,y,z)) 
                    rho.localElement(l) = LHS->Values()[idx++];
            }
        }
    }
    IpplTimings::stopTimer(FunctionTimer4_m);
    
    switch(precmode_m) {
        case STD_PREC:
            delete MLPrec;
            break;
    }

}

///////////////////////////////////////// MAP CREATION //////////////////////////////////////////

//IFF FIXME: think i should move this to the appropriate BoundaryPoints class
inline Epetra_Map* MGPoissonSolver::IPPLToMap3D(NDIndex<3> localidx)
{

    // Version for Elliptic Boundary Points (?)
    int NumMyElements = 0; 
    vector<int> MyGlobalElements;

    for(int x = localidx[0].first(); x <= localidx[0].last(); x++)
        for(int y = localidx[1].first(); y <= localidx[1].last(); y++)
            for(int z = localidx[2].first(); z <= localidx[2].last(); z++) {
                if(bp->isInside(x,y,z)) {
                    MyGlobalElements.push_back(bp->getIdx(x,y,z));
                    NumMyElements++;
                }
            }
    
    return(new Epetra_Map (-1, NumMyElements, &MyGlobalElements[0], 0, Comm));

}

///////////////////////////////////// STENCIL CREATION //////////////////////////////////////////

inline void MGPoissonSolver::Stencil3DGeometry(Vector_t hr, Epetra_CrsMatrix*& Matrix, Epetra_Vector& RHS)
{
    Matrix->PutScalar(0.0);

    int NumMyElements = Map->NumMyElements();
    int* MyGlobalElements = Map->MyGlobalElements();
    int i=0;

    vector<double> Values(6);
    vector<int> Indices(6);

    for (int i = 0 ; i < NumMyElements ; ++i) 
    {

        int NumEntries = 0;
        double WV, EV, SV, NV, FV, BV, CV, scaleFactor;
        double W, E, S, N, F, B;
        bp->getBoundaryStencil(MyGlobalElements[i], WV, EV, SV, NV, FV, BV, CV, scaleFactor);
        bp->getNeighbours(MyGlobalElements[i], W, E, S, N, F, B);

        // scale RHS
        RHS[i] *= scaleFactor;

        if(E != -1) {
            Indices[NumEntries] = E;
            Values[NumEntries++] = EV;
        }
        if(W != -1) {
            Indices[NumEntries] = W;
            Values[NumEntries++] = WV;
        }
        if(S != -1) {
            Indices[NumEntries] = S;
            Values[NumEntries++] = SV;
        }
        if(N != -1) {
            Indices[NumEntries] = N;
            Values[NumEntries++] = NV;
        }
        if(F != -1) {
            Indices[NumEntries] = F;
            Values[NumEntries++] = FV;
        }
        if(B != -1) {
            Indices[NumEntries] = B;
            Values[NumEntries++] = BV;
        }

		// if matrix has allready been filled (FillComplete()) we can only
		// replace entries
        if(Matrix->Filled()) {
            // off-diagonal entries
            Matrix->ReplaceGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
            // diagonal entry
            Matrix->ReplaceGlobalValues(MyGlobalElements[i], 1, &CV, MyGlobalElements + i);
        } else {
            // off-diagonal entries
            Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
            // diagonal entry
            Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &CV, MyGlobalElements + i);
        }


    }

    // convert GIDs to LIDs
    Matrix->FillComplete();
    Matrix->OptimizeStorage();

}

Inform &MGPoissonSolver::print(Inform &os) const
{
    os << "* *************** M G P o i s s o n S o l v e r ************************************ " << endl;
    os << "* h " << hr_m << '\n';
    os << "* ********************************************************************************** " << endl;
}

#endif /* HAVE_ML_SOLVER */
