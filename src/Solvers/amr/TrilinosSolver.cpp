#include <EpetraExt_MatrixMatrix.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <ParallelDescriptor.H>

#include "TrilinosSolver.h"

void
Solver::SetupProblem(PArray<MultiFab>& rhs , PArray<MultiFab>& soln) 
{

    // Now we create the map
    // Note: 1) it is important to set  numMyGridPoints = 0 outside the MFIter loop
    //       2) "bx" here is the cell-centered non-overlapping box
    int nlevs = rhs.size();
    int numMyGridPoints = 0;
    std::vector<int> MyGlobalElements;

    // This offset is the same for all grids at all levels on this processor
    int my_offset = proc_offset[ParallelDescriptor::MyProc()];

    for (int lev = 0; lev < nlevs; lev++) 
    {
       for (MFIter mfi(rhs[lev]); mfi.isValid(); ++mfi)
       {
          IArrayBox& fab_mask =  mask[lev][mfi];
          IArrayBox& fab_imap =  imap[lev][mfi];

          const Box& bx = mfi.validbox();

          for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; z++) 
          for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; y++) 
          for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; x++) 
          {
             // Only add this point to the list of unknowns if mask == 0
             IntVect iv(x,y,z);
             if (fab_mask(iv,0) == INSIDE)
             {
                 MyGlobalElements.push_back(fab_imap(iv,0));
                 numMyGridPoints++;
             }
         }
       }
    }

    // Define the Map
    Map = new Epetra_Map(-1, numMyGridPoints, &MyGlobalElements[0], 0, Comm_m);

    // Allocate the RHS and LHS with the new Epetra Map
    RHS = rcp(new Epetra_Vector(*Map));
    LHS = rcp(new Epetra_Vector(*Map));
    
    int local_idx;

    // Copy the values from rhs into RHS->Values()
    // Copy the values from soln into LHS->Values()
    // Note: "bx" here is the cell-centered non-overlapping box
    for (int lev = 0; lev < nlevs; lev++) 
    {
       for (MFIter mfi(rhs[lev]); mfi.isValid(); ++mfi)
       {
          const Box& bx = mfi.validbox();
          FArrayBox& fab_rhs  =  rhs[lev][mfi];
          FArrayBox& fab_lhs  = soln[lev][mfi];
          IArrayBox& fab_mask = mask[lev][mfi];
          IArrayBox& fab_imap = imap[lev][mfi];

          for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; z++) 
          for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; y++) 
          for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; x++) 
          {
             IntVect iv(x,y,z);

             if (fab_mask(iv,0) == INSIDE)
             {
                 local_idx = fab_imap(iv,0) - my_offset;
    
                 // Set rhs and lhs
                 RHS->Values()[local_idx] = fab_rhs(iv,0);
                 LHS->Values()[local_idx] = fab_lhs(iv,0);
             }
          }
       }
    }

#if 0
    if(verbose_m)
        this->printLoadBalanceStats();
#endif

    // We change the 7 to a 10 because at a coarse-fine boundary
    // the coarse cell sees 4 fine cells instead of 1 coarse cell.
    A = rcp(new Epetra_CrsMatrix(Copy, *Map, 10));
    ComputeStencil();
    std::cout << "Made it out of ComputeStencil " << std::endl;
}

void Solver::extrapolateLHS() {
// Aitken-Neville
// Pi0 (x) := yi , i = 0 : n
// Pik (x) := (x − xi ) Pi+1,k−1(x) − (x − xi+k ) Pi,k−1(x) /(xi+k − xi )
// k = 1, . . . , n, i = 0, . . . , n − k.

    std::deque< Epetra_Vector >::iterator it = OldLHS.begin();

    if(nLHS_m == 0)
        LHS->PutScalar(1.0);
    else if(OldLHS.size() == 1)
        *LHS = *it;
    else if(OldLHS.size() == 2){
        LHS->Update (2.0, *it++, -1.0, *it, 0.0);
    }
    else if(OldLHS.size() > 0)
    {
        int n = OldLHS.size();
        for(int i=0; i<n; ++i){
            *(*P)(i) = *it++;
        }
        for(int k = 1; k < n; ++k){// x==0 & xi==i+1
            for(int i = 0; i < n-k; ++i){
                (*P)(i)->Update(-(i+1)/(float)k, *(*P)(i+1), (i+k+1)/(float)k);//TODO test
            }
        }
        *LHS = *(*P)(0);
    }
    else
    {
        std::cout << "Invalid number of old LHS: " + OldLHS.size() << std::endl;
    }
}

void Solver::CopySolution(PArray<MultiFab>& soln)
{
    int nlevs     = soln.size();
    int my_offset = proc_offset[ParallelDescriptor::MyProc()];

    for (int lev = 0; lev < nlevs; lev++) 
    {
       // Initialize to zero since we only copy when cells are INSIDE
       soln[lev].setVal(0.0);

       for (MFIter mfi(soln[lev]); mfi.isValid(); ++mfi)
       {

          const Box& bx = mfi.validbox();
          FArrayBox& fab_lhs  = soln[lev][mfi];
          IArrayBox& fab_mask = mask[lev][mfi];
          IArrayBox& fab_imap = imap[lev][mfi];

          for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; z++) 
          for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; y++) 
          for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; x++) 
          {
             IntVect iv(x,y,z);
             if (fab_mask(iv,0) == INSIDE)
             {
                 int local_idx = fab_imap(iv,0) - my_offset;
                 fab_lhs(iv,0) = LHS->Values()[local_idx];
             }
         }
       }
    }
}

int 
Solver::getNumIters()
{
    return solver->getNumIters();
}

void 
Solver::Compute()
{
    //LHS->Random();
    //LHS->PutScalar(1.0);
    extrapolateLHS();
    std::cout << "Done with extrapolateLHS " << std::endl; 
    std::cout << " .... and starting to solve " << std::endl; 

    // create the preconditioner object and compute hierarchy
    // true -> create the multilevel hirarchy
    // ML allows the user to cheaply recompute the preconditioner. You can
    // simply uncomment the following line:
    //
    // MLPrec->ReComputePreconditioner();
    //
    // It is supposed that the linear system matrix has different values, but
    // **exactly** the same structure and layout. The code re-built the
    // hierarchy and re-setup the smoothers and the coarse solver using
    // already available information on the hierarchy. A particular
    // care is required to use ReComputePreconditioner() with nonzero
    // threshold.

    if(MLPrec == Teuchos::null) // first repetition we need to create a new preconditioner
        MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
    else if(isReusingHierarchy_m)
        MLPrec->ReComputePreconditioner();
    else if(isReusingPreconditioner_m) {
        // do nothing since we are reusing old preconditioner
    } else { // create a new preconditioner in every repetition
        delete MLPrec.get();//MLPrec now RCP => delete??? TODO
        MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
    }

    // Setup problem
    problem.setOperator(A);
    problem.setLHS(LHS);
    problem.setRHS(RHS);
    prec = rcp(new Belos::EpetraPrecOp(MLPrec));
    problem.setLeftPrec(prec);
    solver->setParameters(rcp(&belosList, false));
    solver->setProblem(rcp(&problem,false));
    if(!problem.isProblemSet()){
        if (problem.setProblem() == false) {
            std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        }
    }

    // Solve problem
    // Timer
    MPI_Barrier(MPI_COMM_WORLD);

    solver->solve();

    // Timer
    MPI_Barrier(MPI_COMM_WORLD);

    // Store new LHS in OldLHS
    OldLHS.push_front(*(LHS.get()));
    if(OldLHS.size() > nLHS_m) OldLHS.pop_back();
    std::cout<<"#OldLHS: "<<OldLHS.size()<<std::endl;
}

void 
Solver::printLoadBalanceStats() {

    //compute some load balance statistics
    size_t myNumPart = Map->NumMyElements();
    size_t NumPart = Map->NumGlobalElements() * 1.0/Comm_m.NumProc();
    double imbalance = 1.0;
    if(myNumPart >= NumPart)
        imbalance += (myNumPart-NumPart)/NumPart;
    else
        imbalance += (NumPart-myNumPart)/NumPart;

    double max=0.0, min=0.0, avg=0.0;
    int minn=0, maxn=0;
    MPI_Reduce(&imbalance, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&imbalance, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&imbalance, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myNumPart, &minn, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myNumPart, &maxn, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    avg /= Comm_m.NumProc();
    if(Comm_m.MyPID() == 0) cout << "LBAL min = " << min << ", max = " << max << ", avg = " << avg << endl;

}

void
Solver::construct_mask_and_imap()
{
    int nlevs = mask.size();
    BoxArray ba;

    int comp = 0;

    // Set the ghost cell values of the level = 0 mask to OUTSIDE
    // Set the interior   values of the level = 0 mask to INSIDE
    int lev = 0;
    for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
    {
        const Box& fab_bx = mask[lev][mfi].box();
        mask[lev][mfi].setVal(OUTSIDE,fab_bx,comp,1);

        const Box& bx = mfi.validbox();
        mask[lev][mfi].setVal(INSIDE,bx,comp,1);
    }

    // Set the ghost cell values of the level > 0 mask to COARSER
    // Set the interior   values of the level > 0 mask to INSIDE
    for (int lev = 1; lev < nlevs; lev++)
    {
        for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
        {
            const Box& fab_bx = mask[lev][mfi].box();
            mask[lev][mfi].setVal(COARSER,fab_bx,comp,1);

            const Box& bx = mfi.validbox();
            mask[lev][mfi].setVal(INSIDE,bx,comp,1);
        }
    }

    // Just in case ...
    // Fill the ghost cells of all IABS 
    for (int lev = 0; lev < nlevs; lev++)
    {
       mask[lev].FillBoundary();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "PROBLO " << prob_lo[0] << " " << prob_lo[1] << " " << prob_lo[2] << std::endl;
        std::cout << "DX     " <<      hr[0] << " " <<      hr[1] << " " <<      hr[2] << std::endl;
    }

    std::multimap < std::pair<int, int>, double >::iterator itr; 

    // At coarsest level we test for cells outside the boundary geometry
    lev = 0;
    BoxArray lev0_ba(mask[0].boxArray());
    Box prob_domain = lev0_ba.minimalBox();

    for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
    {
         const Box& bx = mfi.validbox();
         IArrayBox& fab_mask = mask[lev][mfi];

         // First find lo and hi intersections (xlo, xhi) at fixed i
         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; k++)
         for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; j++)
         {
      	     std::pair<int, int> coordyz(j,k);

       	     itr = xlo.find(coordyz);
             Real inter_xlo = itr->second;     
             
       	     itr = xhi.find(coordyz);
             Real inter_xhi= itr->second;     

      	     // Test if we are outside of domain in x-dir
             for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; i++)
             {
                IntVect iv(i,j,k);
                Real x = prob_lo[0] + (i+.5)*hr[0];
                if (x < inter_xlo               || x > inter_xhi || 
                    i < prob_domain.loVect()[0] || i > prob_domain.hiVect()[0])  
                {
                    fab_mask(iv,0) = OUTSIDE;
                }
             }
         }

         // Next find lo and hi intersections (ylo, yhi) at fixed j
         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; k++)
         for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; i++)
         {
      	     std::pair<int, int> coordxz(i,k);

       	     itr = ylo.find(coordxz);
             Real inter_ylo = itr->second;     

       	     itr = yhi.find(coordxz);
             Real inter_yhi = itr->second;     

      	     // Test if we are outside of domain in y-dir
             for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; j++)
             {
                IntVect iv(i,j,k);
                Real y = prob_lo[1] + (j+.5)*hr[1];
                if (y < inter_ylo               || y > inter_yhi || 
                    j < prob_domain.loVect()[1] || j > prob_domain.hiVect()[1])  
                {
                    fab_mask(iv,0) = OUTSIDE;
                }
             }
         }
    }

    // At all levels but the finest we test to see cells are covered by finer grids
    if (nlevs > 1)
       for (int lev = 0; lev < nlevs-1; lev++)
       {
          ba = mask[lev+1].boxArray();
          ba.coarsen(2);
          for (int b = 0; b < ba.size(); b++)
          {
             mask[lev].setVal(COVERED,ba[b],0,1);
          }
       }

    // Fill the ghost cells of all IABS at each level
    for (int lev = 0; lev < nlevs; lev++)
        mask[lev].FillBoundary();

    // Now define the imap -- first default to -1
    for (int lev = 0; lev < nlevs; lev++)
        imap[lev].setVal(-1);

    // We know there are the same number of processes as grids at level 0
    int nprocs = mask[0].boxArray().size();
    int cells_per_proc[nprocs];
    for (int i = 0; i < nprocs; i++) cells_per_proc[i] = 0;

    // First figure out how many valid points on each processor
    int ct = 0;
    for (int lev = 0; lev < nlevs; lev++)
    {
       for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
       {
         int grid_id = mfi.index();
         IArrayBox& fab_mask =  mask[lev][mfi];
         const Box& bx = mfi.validbox();

         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; k++)
         for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; j++)
         for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; i++)
         {
             IntVect iv(i,j,k);
             if (fab_mask(iv,0) == INSIDE) ct++;
         }
         cells_per_proc[grid_id] = ct;
       }
    }
    // Send cells_per_proc from every processor to every processor
    ParallelDescriptor::ReduceIntSum(cells_per_proc,nprocs);
    ParallelDescriptor::Bcast(cells_per_proc,nprocs);

    proc_offset = new int[nprocs];
    proc_offset[0] = 0;
    for (int i = 0; i < nprocs; i++)
        proc_offset[i] = proc_offset[i-1] + cells_per_proc[i-1];

    //
    // Send proc_offset of every grid to every processor
    //
    ParallelDescriptor::Bcast(proc_offset,nprocs);

    int my_offset = proc_offset[ParallelDescriptor::MyProc()];

    // Make sure to set idx to 0
    int local_idx = 0;
    
    // Here do just the level 0 grids -- these don't have the
    // special factor 2 indexing
    lev = 0;
    {
       for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
       {
         IArrayBox& fab_mask =  mask[lev][mfi];
         IArrayBox& fab_imap =  imap[lev][mfi];
         const Box& bx = mfi.validbox();
         std::cout << "LEVEL 0 BOX " << bx << std::endl;

         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; k++)
         for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; j++)
         for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; i++)
         {
             IntVect iv(i,j,k);
             if (fab_mask(iv,0) == INSIDE)
             {
                 fab_imap(iv,0) = my_offset + local_idx;
                 local_idx++;
             }
         }
       }
    }

    // Now do the level > 0 grids -- these use a special indexing
    // so that if we know just one index we know all indices in that
    // 2x2x2 box.
    for (int lev = 1; lev < nlevs; lev++)
    {
       for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
       {
         IArrayBox& fab_mask =  mask[lev][mfi];
         IArrayBox& fab_imap =  imap[lev][mfi];
         const Box& bx = mfi.validbox();
         std::cout << "LEVEL " << lev << "  BOX " << bx << std::endl;

         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; k+=2)
         for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; j+=2)
         for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; i+=2)
         {
             // DO NOT THE ORDER OF THESE LINES -- this is the order assumed
             //   in ComputeStencil -- if you change it here you must
             //   change the definitions for COVERED cells in ComputeSTencil
             IntVect iv_lll(i  ,j  ,k  );
             IntVect iv_hll(i+1,j  ,k  );
             IntVect iv_lhl(i  ,j+1,k  );
             IntVect iv_hhl(i+1,j+1,k  );
             IntVect iv_llh(i  ,j  ,k+1);
             IntVect iv_hlh(i+1,j  ,k+1);
             IntVect iv_lhh(i  ,j+1,k+1);
             IntVect iv_hhh(i+1,j+1,k+1);

             // Note that if iv_lll is INSIDE then all 8 cells are INSIDE
             if (fab_mask(iv_lll,0) == INSIDE)
             {
                 fab_imap(iv_lll,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_hll,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_lhl,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_hhl,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_llh,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_hlh,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_lhh,0) = my_offset + local_idx; local_idx++;
                 fab_imap(iv_hhh,0) = my_offset + local_idx; local_idx++;
             }
 
             if (i == 30 && j == 30 && k == 30) 
                 std::cout << "CONSTRUCTING IMAP(31,31,31) " << fab_imap(iv_hhh,0) << std::endl; 
         }
       }
    }
    // Replace the imap value of COVERED cells by the imap value of the iv_lll cell
    //         on the finer grid covering it
    for (int lev = 0; lev < nlevs-1; lev++)
        avgDownImap(imap[lev],imap[lev+1]);

    // Replace the imap value of COVERED cells by the imap value of the iv_lll cell
    //         on the finer grid covering it
    for (int lev = 1; lev < nlevs; lev++)
        FillImapGhostCells(imap[lev-1],imap[lev], mask[lev]);

    // Fill the ghost cells of all IABS at this level
    for (int lev = 0; lev < nlevs; lev++)
    {
       imap[lev].FillBoundary();
    }
}

void
Solver::avgDownImap(iMultiFab& crse, iMultiFab& fine)
{
    IntVect ratio(2,2,2);

    BoxArray crse_fine_BA(fine.boxArray().size());

    for (int i = 0; i < fine.boxArray().size(); ++i)
        crse_fine_BA.set(i,BoxLib::coarsen(fine.boxArray()[i],ratio));

    iMultiFab crse_fine(crse_fine_BA,1,0);

    for (MFIter mfi(fine); mfi.isValid(); ++mfi)
    {
        const int        i        = mfi.index();
        const Box&       ovlp     = crse_fine_BA[i];
        IArrayBox&       crse_fab = crse_fine[i];
        const IArrayBox& fine_fab = fine[i];

        // Loop over the coarse indices of the box where the fine and
        // coarse grids overlap
        for (int k = ovlp.loVect()[2]; k <= ovlp.hiVect()[2]; k++)
        for (int j = ovlp.loVect()[1]; j <= ovlp.hiVect()[1]; j++)
        for (int i = ovlp.loVect()[0]; i <= ovlp.hiVect()[0]; i++)
        {   
              IntVect iv_crse(  i,  j,  k);
              IntVect iv_lll (2*i,2*j,2*k);
              crse_fab(iv_crse,0) = fine_fab(iv_lll,0);
        }   
    }

    crse.copy(crse_fine);
}

void
Solver::FillImapGhostCells(iMultiFab& crse_imap, iMultiFab& fine_imap, iMultiFab& fine_mask)
{
    IntVect ratio(2,2,2);

    BoxArray crse_fine_BA(fine_imap.boxArray().size());

    // Create a coarse boxArray that holds the fine boxes PLUS ghost cells
    // This boxArray is at coarse resolution but uses the processor distribution
    // of the fine boxArray.
    for (int i = 0; i < fine_imap.boxArray().size(); ++i)
        crse_fine_BA.set(i,BoxLib::grow(BoxLib::coarsen(fine_imap.boxArray()[i],ratio),1));

    iMultiFab crse_fine(crse_fine_BA,1,0);

    // This boxArray is at coarse resolution but uses the processor distribution
    // We copy the crse imap values into this new MultiFab.
    crse_fine.copy(crse_imap);

    for (MFIter mfi(fine_imap); mfi.isValid(); ++mfi)
    {
        const int        i   = mfi.index();
        const Box&   fine_bx = mfi.validbox();

        const IArrayBox& crse_imap_fab = crse_fine[i];
              IArrayBox& fine_imap_fab = fine_imap[i];
        const IArrayBox& fine_mask_fab = fine_mask[i];

        // Loop over the fine indices including ghost cells -- but only fill
        // from the coarse grid if the mask says COARSER
        for (int k = fine_bx.loVect()[2]-1; k <= fine_bx.hiVect()[2]+1; k++)
        for (int j = fine_bx.loVect()[1]-1; j <= fine_bx.hiVect()[1]+1; j++)
        for (int i = fine_bx.loVect()[0]-1; i <= fine_bx.hiVect()[0]+1; i++)
        {   
              IntVect iv_fine(i,j,k);
              if (fine_mask_fab(iv_fine,0) == COARSER)
              {
                 IntVect iv_crse(i/2,j/2,k/2);
                 fine_imap_fab(iv_fine,0) = crse_imap_fab(iv_crse,0);
              }
        }
    }
}
