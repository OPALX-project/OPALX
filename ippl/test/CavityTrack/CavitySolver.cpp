#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <vector>

#include "rlog/rlog.h"
#include "rlog/StdioNode.h"
#include "rlog/RLogChannel.h"

#include "pbe.h"

#include "Epetra_ConfigDefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Vector.h>
#include "Epetra_Export.h"
#include <ml_epetra_preconditioner.h>

#include "myrlog.h"
#include "nedelecmesh.h"
#include "femaxmesh.h"
#include "femaxxdriver.h"

#include "Const.hh"

#include "CavitySolver.h"

CavitySolver::CavitySolver(string meshfilename, double cavitylength,double 
			   xtrans,double ytrans,double ztrans,double 
			   xrot, double yrot, double zrot, 
			   double scaling,double peakvoltage,double phase) : 
#ifdef HAVE_MPI
    comm_(new Epetra_MpiComm(MPI_COMM_WORLD))
#else
    comm_(new Epetra_SerialComm())
#endif 
{
    meshfile_ = meshfilename;
    cavitylength_ = cavitylength;
    xtrans_ = xtrans;
    ytrans_ = ytrans;
    ztrans_ = ztrans;
    xrot_ = xrot;
    yrot_ = yrot;
    zrot_ = zrot;
    scaling_ = scaling;
    peakvoltage_ = peakvoltage;
    phase_ = phase;
    factor_ = 1.0; // this variable is used to set the eigenfields look in the correct direction
	
    //
    // Rlog initialisation
    //
    const bool enable_log = true;
    stdLog_ = new rlog::MyStdioNode();
    dbgLog_ = 0;

    DEF_CHANNEL("warning/all", rlog::Log_Warning);
    DEF_CHANNEL("error/all", rlog::Log_Error);
    DEF_CHANNEL("info/all", rlog::Log_Info);
    DEF_CHANNEL("debug/all", rlog::Log_Debug);

    if (comm_->MyPID() == 0) {
        stdLog_->subscribeTo(RLOG_CHANNEL("warning"));
        stdLog_->subscribeTo(RLOG_CHANNEL("error"));
        stdLog_->subscribeTo(RLOG_CHANNEL("info"));
    } else {
        stdLog_->subscribeTo(RLOG_CHANNEL("error"));
        stdLog_->subscribeTo(RLOG_CHANNEL("warning/all"));
        stdLog_->subscribeTo(RLOG_CHANNEL("info/all"));
    }

    string logname(get_log_file_name("cavity_solver_"));
    if (enable_log && comm_->MyPID() == 0) {
        int fid = ::open(logname.c_str(), O_WRONLY|O_CREAT|O_TRUNC, 0666);
        rAssert(fid != -1);
        dbgLog_ = new rlog::StdioNode(fid);
        dbgLog_->subscribeTo(RLOG_CHANNEL("debug"));
        dbgLog_->subscribeTo(RLOG_CHANNEL("info"));
        dbgLog_->subscribeTo(RLOG_CHANNEL("warning"));
        dbgLog_->subscribeTo(RLOG_CHANNEL("error"));
    }
    if (enable_log)
        rInfo("Logging debug output to %s", logname.c_str());
}

CavitySolver::~CavitySolver() {
    rInfo("Femaxx/CavitySolver shutting down...");
    delete driver_;
    delete stdLog_;
    delete dbgLog_;
    delete comm_;
}

void CavitySolver::runSolver() {
    rInfo("Femaxx/CavitySolver running on %d processor(s)...", comm_->NumProc());
    rLog(RLOG_CHANNEL("info/all"), "Hi from processor %d.", comm_->MyPID());

    //
    // Build femaXX parameter list.
    //
    Teuchos::ParameterList* params = new Teuchos::ParameterList();
           
    // femaXX default parameter values
    FemaxxDriver::set_defaults(*params);

    // Custom parameters
    params->set("mesh_file_name", meshfile_ );
    params->set("sigma", 1.5);
    params->set("kmax", 2);
    params->sublist("jd_params").set("tol", 1e-5);
#define FINITE_ELEMENT_ORDER 2
#if FINITE_ELEMENT_ORDER == 1
    params->sublist("asigma_precon").set("type", "ml");
#else
    params->set("element_order", 2);
    params->sublist("asigma_precon").set("type", "2level");
    params->sublist("asigma_precon").set("solver11", "ml");
#endif

    // Log parameter list.
    if (comm_->MyPID() == 0) {
        ostringstream buf;
        buf << "Input parameters:" << endl;
        params->print(buf, 8);
        rDebug("%s", buf.str().c_str());
    }

    //
    // Run eigenmode solver.
    //
    driver_ = new FemaxxDriver(*comm_, *params); 
    driver_->run();
  
    // Get computed eigenvectors and eigenpairs.
    //
    // Note that the Epetra_MultiVector returned by this function shares data
    // with class-internal data structures. Do not access the Epetra_MultiVector after
    // FemaxxDriver has been destroyed.
    Epetra_MultiVector Q(driver_->get_eigenvectors());
    Epetra_SerialDenseVector L(driver_->get_eigenvalues());

    //
    // Replicate eigenpairs on all processors.
    //
    int nof_converged = Q.NumVectors();
    int nof_dofs = Q.GlobalLength();
    rAssert(nof_converged > 0);
    eigenvectors_.reserve(nof_converged);
    eigenvalues_.reserve(nof_converged);
    for (int k = 0; k < nof_converged; ++ k) {
        eigenvectors_[k] = colarray::Vector<double>(nof_dofs);
        get_eigenpair_all(k, L, Q, eigenvalues_[k], eigenvectors_[k]);
    }
}

void CavitySolver::get_eigenpair_all(int k,     
                                     Epetra_SerialDenseVector& L,
                                     Epetra_MultiVector& Q,
                                     double& lambda, 
                                     colarray::Vector<double>& q) 
{
    assert(0 <= k && k < Q.NumVectors());
    assert(comm_->MyPID() > 0 || q._n == Q.GlobalLength());

    //
    // Bring k-th eigenvector to processor 0
    //

    // create a target map, for which all the elements are on proc 0
    Epetra_MultiVector q_dist(View, Q, k, 1);
    int NumMyElements_target = 0;
    if (comm_->MyPID() == 0)
        NumMyElements_target = Q.GlobalLength();
    Epetra_Map TargetMap(-1, NumMyElements_target, 0, *comm_);

    // redistribute vector to processor 0
    Epetra_Vector q0(TargetMap);
    q0.Export(q_dist, Epetra_Export(q_dist.Map(), TargetMap), Add);

    // store eigenpair on processor 0
    if (comm_->MyPID() == 0) {
        lambda = L[k];
        for (int i = 0; i < Q.GlobalLength(); i ++) {
            q(i) =  q0[i];
        }
    }

    //
    // Broadcast eigenpair to all processors.
    //
    comm_->Broadcast(&lambda, 1, 0);
    comm_->Broadcast(q._v, Q.GlobalLength(), 0);    
}

void CavitySolver::getECavity(double* pos, double* ECavity, unsigned long k, double t) {
    double lambda = eigenvalues_[k];
    colarray::Vector<double>& q(eigenvectors_[k]);

    // now transform pos into the femaxx coordinates from FStest
    // ones using the parameters given this is not needed at the
    // moment	
    // transform_coords(pos);
	
    // make a vector3 of the position
    mesh::Vector3 coord(pos[0], pos[1], pos[2]);
    
    FemaxMesh& femaxmesh = driver_->get_femax_mesh();
    NedelecMesh& nedelec_mesh = femaxmesh.get_nedelec_mesh();
    mesh::Vector3 e_field = nedelec_mesh.eval(coord, q);
    
    // this is the hack to get the eigenvector in the right direction (only for box cavity)
    if (e_field.x <= 0) {
        factor_ = -1.0;
    }
    e_field *= factor_;

    // Rescale the E field of the cavity: This ensures that the gap
    // voltage is equal to peakvoltage_, provided that the (amplitude
    // of) E-field is constant along the gap.
    e_field *= peakvoltage_ / cavitylength_;
    
    // Consider the time-dependency of the E field
    double omega = sqrt(lambda * pow(CLIGHT,2));    
    e_field *= cos((omega * t) + phase_);

    ECavity[0] = e_field.x;
    ECavity[1] = e_field.y;
    ECavity[2] = e_field.z;
}

void CavitySolver::getBCavity(double* pos, double* BCavity, unsigned long k, double t){
    double lambda = eigenvalues_[k];
    colarray::Vector<double>& q(eigenvectors_[k]);

    // now transform pos into the femaxx coordinates from FStest
    // ones using the parameters given this is not needed at the
    // moment    
    // transform_coords(pos);
    
    // make a vector3 from of the position
    mesh::Vector3 coord(pos[0],pos[1],pos[2]);
    
    FemaxMesh* femaxmesh = &driver_->get_femax_mesh();
    NedelecMesh* nedelec_mesh = &femaxmesh->get_nedelec_mesh();
    mesh::Vector3 b_field = nedelec_mesh->eval_curl(coord, q);
    
    // this is the hack for getting the right direction (provided the getEfield was called before)
    b_field *= factor_;
    
    // Rescale the B field of the cavity
    double omega = sqrt(lambda * pow(CLIGHT,2));	
    b_field *= (CLIGHT/omega) * (peakvoltage_ / cavitylength_) * MUE0;
    
    // Consider the time-dependency of the B field
    b_field *= sin((omega * t) + phase_);

    BCavity[0] = b_field.x;
    BCavity[1] = b_field.y;
    BCavity[2] = b_field.z;
}

void CavitySolver::transform_coords(double* coords) {
    // translation
    coords[0]=coords[0] + xtrans_;
    coords[1]=coords[1] + ytrans_;
    coords[2]=coords[2] + ztrans_;


    // uniform homogenous scaling
    coords[0]=coords[0] * scaling_;
    coords[1]=coords[1] * scaling_;
    coords[2]=coords[2] * scaling_;

    // rotation about the axes x, y and z with given rotation angles xrot_, yrot_, zrot_

    // around the x axis
    coords[1]=coords[1] * cos(xrot_) - coords[2] * sin(xrot_);
    coords[2]=coords[1] * sin(xrot_) + coords[2] * cos(xrot_);

    // around the y axis
    coords[2]=coords[2] * cos(yrot_) - coords[0] * sin(yrot_);
    coords[0]=coords[2] * sin(yrot_) + coords[0] * cos(yrot_);

    // around the z-axis
    coords[0]=coords[0] * cos(zrot_) - coords[1] * sin(zrot_);
    coords[1]=coords[0] * sin(zrot_) + coords[1] * cos(zrot_);

}
