#include <ParallelDescriptor.H>
#include "TrilinosSolver.h"

void 
Solver::ComputeStencil()
{
    int NumMyElements = Map->NumMyElements();
    printf("NumMyElements = %d\n",NumMyElements);

    int* MyGlobalElements = Map->MyGlobalElements();

    std::vector<double> Values(9);
    std::vector<int> Indices(9);

    int nlevs = mask.size();

    // The stencil at a coarse cell (x_i) next to a fine cell on the right and coarse cell on the left is
    //    [ (32/21) x_fine + (24/21) x_{i-1} - (56/21) x_i ] / (h_crse * h_crse)
    // or [ (32/21) x_fine + ( 8/7 ) x_{i-1} - ( 8/3 ) x_i ] / (h_crse * h_crse)
    //
    // The stencil at a fine cell (x_i) next to a fine cell on the left and coarse cell on the right is
    //    [ (8/15) x_crse + (12/15) x_{i-1} - (20/15) x_i ]  / (h_fine * h_fine)
    // or [ (8/15) x_crse + ( 4/5 ) x_{i-1} - ( 4/3 ) x_i ]  / (h_fine * h_fine)
    // or [ (32/15) x_crse + (16/5 ) x_{i-1} - (16/3 ) x_i ]  / (h_crse * h_crse)

    for (int lev = 0; lev < nlevs; lev++) 
    {
       if (ParallelDescriptor::IOProcessor()) 
           std::cout << "... Starting to compute stencil at level " << lev << std::endl;

       for (MFIter mfi(mask[lev]); mfi.isValid(); ++mfi)
       {
          IArrayBox& fab_mask =  mask[lev][mfi];
          IArrayBox& fab_imap =  imap[lev][mfi];

          const Box& bx = mfi.validbox();

          for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; z++) 
          for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; y++) 
          for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; x++) 
          {
             IntVect iv(x,y,z);
             int global_idx = fab_imap(iv,0);

             if (lev == 1 && x == 31 && y == 31 && z == 31) std::cout << "IMAP IN COMPUTE_STENCIL " << fab_imap(iv,0) << std::endl;

             if (fab_mask(iv,0) == INSIDE)
             {
                 IntVect iv_w(x-1,y,z);
                 IntVect iv_e(x+1,y,z);
                 IntVect iv_s(x,y-1,z);
                 IntVect iv_n(x,y+1,z);
                 IntVect iv_f(x,y,z-1);
                 IntVect iv_b(x,y,z+1);

                 int NumEntries = 0;

                 Real CV = 0.;

                 // Sanity checks!
                 // If one diretion is not INSIDE then the opposite direction must be. 
                 if (fab_mask(iv_w) != INSIDE && fab_mask(iv_e) != INSIDE) 
                 {
                    std::cout << "AT IV " << iv << std::endl;
                    std::cout << "EAST MASK IS " << fab_mask(iv_e) << std::endl;
                    std::cout << "WEST MASK IS " << fab_mask(iv_w) << std::endl;
                    BoxLib::Error("East and West both not INSIDE!!");
                 }
                 if (fab_mask(iv_s) != INSIDE && fab_mask(iv_n) != INSIDE) 
                 {
                    std::cout << "AT IV " << iv << std::endl;
                    std::cout << "SOUTH MASK IS " << fab_mask(iv_s) << std::endl;
                    std::cout << "NORTH MASK IS " << fab_mask(iv_n) << std::endl;
                    BoxLib::Error("South and North both not INSIDE!!");
                 }
                 if (fab_mask(iv_f) != INSIDE && fab_mask(iv_b) != INSIDE) 
                 {
                    std::cout << "AT IV " << iv << std::endl;
                    std::cout << "FRONT MASK IS " << fab_mask(iv_f) << std::endl;
                    std::cout << "BACK  MASK IS " << fab_mask(iv_b) << std::endl;
                    BoxLib::Error("Front and Back both not INSIDE!!");
                 }

                 // If INSIDE or OUTSIDE in both directions add "usual" contribution
                 if ( (fab_mask(iv_e) == INSIDE || fab_mask(iv_e) == OUTSIDE) && 
                      (fab_mask(iv_w) == INSIDE || fab_mask(iv_w) == OUTSIDE) )
                     CV += 2/(hr[0]*hr[0]);
                 if ( (fab_mask(iv_s) == INSIDE || fab_mask(iv_s) == OUTSIDE) && 
                      (fab_mask(iv_n) == INSIDE || fab_mask(iv_n) == OUTSIDE) )
                     CV += 2/(hr[1]*hr[1]);
                 if ( (fab_mask(iv_f) == INSIDE || fab_mask(iv_f) == OUTSIDE) && 
                      (fab_mask(iv_b) == INSIDE || fab_mask(iv_b) == OUTSIDE) )
                     CV += 2/(hr[2]*hr[2]);
                 
                 // If COVERED in one direction add coarse/fine contribution
                 if (fab_mask(iv_e) == COVERED || fab_mask(iv_w) == COVERED) 
                     CV += (8./3.)/(hr[0]*hr[0]);
                 if (fab_mask(iv_s) == COVERED || fab_mask(iv_n) == COVERED) 
                     CV += (8./3.)/(hr[1]*hr[1]);
                 if (fab_mask(iv_f) == COVERED || fab_mask(iv_b) == COVERED) 
                     CV += (8./3.)/(hr[2]*hr[2]);
                 
                 // If COARSER in one direction add fine/coarse contribution
                 if (fab_mask(iv_e) == COARSER || fab_mask(iv_w) == COARSER) 
                     CV += (16./3.)/(hr[0]*hr[0]);
                 if (fab_mask(iv_s) == COARSER || fab_mask(iv_n) == COARSER) 
                     CV += (16./3.)/(hr[1]*hr[1]);
                 if (fab_mask(iv_f) == COARSER || fab_mask(iv_b) == COARSER) 
                     CV += (16./3.)/(hr[2]*hr[2]);

                 // **************************************************************************
                 // WEST
                 // **************************************************************************
                 // If the cell is adjacent at the same level
                 if (fab_mask(iv_w) == INSIDE) {

                     Indices[NumEntries] = fab_imap(iv_w,0);
                     if (fab_mask(iv_e) == COVERED) {
                        Values[NumEntries++] = -(8./7.)/(hr[0]*hr[0]);
                     } else if (fab_mask(iv_e) == COARSER) {
                        Values[NumEntries++] = -(16./5.)/(hr[0]*hr[0]);
                     } else {
                        Values[NumEntries++] = -1/(hr[0]*hr[0]);
                     }

                 // If the cell is outside the domain, then we do not use it in the stencil
                 //    because the value there is 0
                 } else if (fab_mask(iv_w) == OUTSIDE) {

                 // If the cell is covered by a finer grid -- note that fab_imap(iv_w)
                 //     is in fact the index of the LLL fine cell covering this coarse cell
                 } else if (fab_mask(iv_w) == COVERED) {
#if 1
                     Indices[NumEntries] = fab_imap(iv_w)+1;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);

                     Indices[NumEntries] = fab_imap(iv_w)+3;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);

                     Indices[NumEntries] = fab_imap(iv_w)+5;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);

                     Indices[NumEntries] = fab_imap(iv_w)+7;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);
#endif

                 //
                 // If the cell is on a coarser grid, then the ghost cell value of the fine
                 //     grid is in fact the index of the coarse cell there
                 } else if (fab_mask(iv_w) == COARSER) {
                     Indices[NumEntries] = fab_imap(iv_w);
                     Values[NumEntries++] = -(32./15.)/(hr[0]*hr[0]);
                 } 

                 // **************************************************************************
                 // EAST
                 // **************************************************************************
                 // If the cell is adjacent at the same level
                 if (fab_mask(iv_e) == INSIDE) {
                     Indices[NumEntries] = fab_imap(iv_e,0);
                     if (fab_mask(iv_w) == COVERED) {
                        Values[NumEntries++] = -(8./7.)/(hr[0]*hr[0]);
                     } else if (fab_mask(iv_w) == COARSER) {
                        Values[NumEntries++] = -(16./5.)/(hr[0]*hr[0]);
                     } else {
                        Values[NumEntries++] = -1/(hr[0]*hr[0]);
                     }

                 // If the cell is outside the domain, then we do not use it in the stencil
                 //    because the value there is 0
                 } else if (fab_mask(iv_e) == OUTSIDE) {

                 // If the cell is covered by a finer grid -- note that fab_imap(iv_e)
                 //     is in fact the index of the LLL fine cell covering this coarse cell
                 } else if (fab_mask(iv_e) == COVERED) {

                     Indices[NumEntries] = fab_imap(iv_e)+0;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);
#if 1
                     Indices[NumEntries] = fab_imap(iv_e)+2;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);

                     Indices[NumEntries] = fab_imap(iv_e)+4;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);

                     Indices[NumEntries] = fab_imap(iv_e)+6;
                     Values[NumEntries++] = -(8./21.)/(hr[0]*hr[0]);
#endif
                 
                 //
                 // If the cell is on a coarser grid, then the ghost cell value of the fine
                 //     grid is in fact the index of the coarse cell there
                 } else if (fab_mask(iv_e) == COARSER) {
                     Indices[NumEntries] = fab_imap(iv_e);
                     Values[NumEntries++] = -(32./15.)/(hr[0]*hr[0]);
                 } 
                 
                 // **************************************************************************
                 // SOUTH
                 // **************************************************************************
                 // If the cell is adjacent at the same level
                 if (fab_mask(iv_s) == INSIDE) {
                     Indices[NumEntries] = fab_imap(iv_s,0);
                     if (fab_mask(iv_n) == COVERED) {
                        Values[NumEntries++] = -(8./7.)/(hr[1]*hr[1]);
                     } else if (fab_mask(iv_n) == COARSER) {
                        Values[NumEntries++] = -(16./5.)/(hr[1]*hr[1]);
                     } else {
                        Values[NumEntries++] = -1/(hr[1]*hr[1]);
                     }

                 // If the cell is outside the domain, then we do not use it in the stencil
                 //    because the value there is 0
                 } else if (fab_mask(iv_s) == OUTSIDE) {

                 // If the cell is covered by a finer grid -- note that fab_imap(iv_s)
                 //     is in fact the index of the LLL fine cell covering this coarse cell
                 } else if (fab_mask(iv_s) == COVERED) {
#if 1
                     Indices[NumEntries] = fab_imap(iv_s)+2;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);

                     Indices[NumEntries] = fab_imap(iv_s)+3;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);

                     Indices[NumEntries] = fab_imap(iv_s)+6;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);

                     Indices[NumEntries] = fab_imap(iv_s)+7;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);
#endif
                 
                 //
                 // If the cell is on a coarser grid, then the ghost cell value of the fine
                 //     grid is in fact the index of the coarse cell there
                 } else if (fab_mask(iv_s) == COARSER) {
                     Indices[NumEntries] = fab_imap(iv_s);
                     Values[NumEntries++] = -(32./15.)/(hr[1]*hr[1]);
                 } 
                  
                 // **************************************************************************
                 // NORTH
                 // **************************************************************************
                 // If the cell is adjacent at the same level
                 if (fab_mask(iv_n) == INSIDE) {
                     Indices[NumEntries] = fab_imap(iv_n,0);
                     if (fab_mask(iv_s) == COVERED) {
                        Values[NumEntries++] = -(8./7.)/(hr[1]*hr[1]);
                     } else if (fab_mask(iv_s) == COARSER) {
                        Values[NumEntries++] = -(16./5.)/(hr[1]*hr[1]);
                     } else {
                        Values[NumEntries++] = -1/(hr[1]*hr[1]);
                     }

                 // If the cell is outside the domain, then we do not use it in the stencil
                 //    because the value there is 0
                 } else if (fab_mask(iv_n) == OUTSIDE) {

                 // If the cell is covered by a finer grid -- note that fab_imap(iv_n)
                 //     is in fact the index of the LLL fine cell covering this coarse cell
                 } else if (fab_mask(iv_n) == COVERED) {

                     Indices[NumEntries] = fab_imap(iv_n)+0;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);

#if 1
                     Indices[NumEntries] = fab_imap(iv_n)+1;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);

                     Indices[NumEntries] = fab_imap(iv_n)+4;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);

                     Indices[NumEntries] = fab_imap(iv_n)+5;
                     Values[NumEntries++] = -(8./21.)/(hr[1]*hr[1]);
#endif
                 
                 //
                 // If the cell is on a coarser grid, then the ghost cell value of the fine
                 //     grid is in fact the index of the coarse cell there
                 } else if (fab_mask(iv_n) == COARSER) {
                     Indices[NumEntries] = fab_imap(iv_n);
                     Values[NumEntries++] = -(32./15.)/(hr[1]*hr[1]);
                 } 
                  
                 // **************************************************************************
                 // FRONT
                 // **************************************************************************
                 // If the cell is adjacent at the same level
                 if (fab_mask(iv_f) == INSIDE) {
                     Indices[NumEntries] = fab_imap(iv_f,0);
                     if (fab_mask(iv_b) == COVERED) {
                        Values[NumEntries++] = -(8./7.)/(hr[2]*hr[2]);
                     } else if (fab_mask(iv_b) == COARSER) {
                        Values[NumEntries++] = -(16./5.)/(hr[2]*hr[2]);
                     } else {
                        Values[NumEntries++] = -1/(hr[2]*hr[2]);
                     }

                 // If the cell is outside the domain, then we do not use it in the stencil
                 //    because the value there is 0
                 } else if (fab_mask(iv_f) == OUTSIDE) {

                 // If the cell is covered by a finer grid -- note that fab_imap(iv_f)
                 //     is in fact the index of the LLL fine cell covering this coarse cell
                 } else if (fab_mask(iv_f) == COVERED) {

                     Indices[NumEntries] = fab_imap(iv_f)+0;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);
#if 1
                     Indices[NumEntries] = fab_imap(iv_f)+1;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);

                     Indices[NumEntries] = fab_imap(iv_f)+4;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);

                     Indices[NumEntries] = fab_imap(iv_f)+5;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);
#endif
                 
                 //
                 // If the cell is on a coarser grid, then the ghost cell value of the fine
                 //     grid is in fact the index of the coarse cell there
                 } else if (fab_mask(iv_f) == COARSER) {
                     Indices[NumEntries] = fab_imap(iv_f);
                     Values[NumEntries++] = -(32./15.)/(hr[2]*hr[2]);
                 } 
                 
                 // **************************************************************************
                 // BACK
                 // **************************************************************************
                 // If the cell is adjacent at the same level
                 if (fab_mask(iv_b) == INSIDE) {
                     Indices[NumEntries] = fab_imap(iv_b,0);
                     if (fab_mask(iv_f) == COVERED) {
                        Values[NumEntries++] = -(8./7.)/(hr[2]*hr[2]);
                     } else if (fab_mask(iv_f) == COARSER) {
                        Values[NumEntries++] = -(16./5.)/(hr[2]*hr[2]);
                     } else {
                        Values[NumEntries++] = -1/(hr[2]*hr[2]);
                     }

                 // If the cell is outside the domain, then we do not use it in the stencil
                 //    because the value there is 0
                 } else if (fab_mask(iv_b) == OUTSIDE) {

                 // If the cell is covered by a finer grid -- note that fab_imap(iv_b)
                 //     is in fact the index of the LLL fine cell covering this coarse cell
                 } else if (fab_mask(iv_b) == COVERED) {

                     Indices[NumEntries] = fab_imap(iv_b)+0;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);
#if 1
                     Indices[NumEntries] = fab_imap(iv_b)+1;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);

                     Indices[NumEntries] = fab_imap(iv_b)+4;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);

                     Indices[NumEntries] = fab_imap(iv_b)+5;
                     Values[NumEntries++] = -(8./21.)/(hr[2]*hr[2]);
#endif
                 
                 //
                 // If the cell is on a coarser grid, then the ghost cell value of the fine
                 //     grid is in fact the index of the coarse cell there
                 } else if (fab_mask(iv_b) == COARSER) {
                     Indices[NumEntries] = fab_imap(iv_b);
                     Values[NumEntries++] = -(32./15.)/(hr[2]*hr[2]);
                 } 

                 if (lev == nlevs-1 && NumEntries > 6)
                 {
                    std::cout << " AT CELL " << iv << " AT LEVEL " << lev << std::endl;
                    std::cout << " NUMENTRIES " << NumEntries << std::endl;
                    BoxLib::Error("Should only have 6 entries at the finest level in ComputeStencil!!!");
                 }

                 for (int ii = 0; ii < NumEntries; ii++) 
                 {
                     if (Indices[ii] < 0)
                     {
                         std::cout << "INDICES AT GLOBAL IDX " << global_idx << 
                                      " AND ENTRY " << ii << " IS " << Indices[ii] << std::endl; 
                         BoxLib::Error("Bad Indices in ComputeStencil!!! ");
                     }
                 }

                 if (lev == 1 && x == 31 && y == 31 && z == 31)
                 {
                     std::cout << "CV AT IDX " << global_idx << " " << CV << std::endl;
                     for (int ii = 0; ii < NumEntries; ii++) 
                         std::cout << "INDICES AT IDX " << global_idx << " AND ENTRY " << ii << " IS " << 
                                   Indices[ii] << " " << Values[ii] <<  std::endl; 
                 }

                 // put the off-diagonal entries
                 A->InsertGlobalValues(MyGlobalElements[global_idx], NumEntries, &Values[0], &Indices[0]);

                 // put in the diagonal entry
                 A->InsertGlobalValues(MyGlobalElements[global_idx], 1, &CV, &(MyGlobalElements[global_idx]));
             }
          }
       }
    }

    A->FillComplete();
    A->OptimizeStorage();
}
