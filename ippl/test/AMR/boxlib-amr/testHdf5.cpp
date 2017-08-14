/*  
 *  This example writes dataset sing chunking. Each process writes
 *  exactly one chunk.
 *             - | 
 *             * V 
 *  Number of processes is assumed to be 4.
 */
 
#include "Ippl.h"
#include "hdf5.h"
#include "stdlib.h"

#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <cassert>
#include <memory>

#ifdef UNIQUE_PTR
typedef Array<std::unique_ptr<MultiFab> > container_t;
#else
#include <PArray.H>
typedef PArray<MultiFab> container_t;
#endif

#define FILENAME    "extend.h5"
#define DATASETNAME "ExtendibleArray"

int main (int argc, char **argv) {
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc, argv, false);
    
    assert(BL_SPACEDIM == 3);
    
    int nr = std::atoi(argv[1]);
    int size = std::atoi(argv[2]);
    
    std::cout << nr << "x" << nr << "x" << nr << std::endl
              << "max. box size: " << size << std::endl;
    
    IntVect low(0, 0, 0);
    IntVect high(nr - 1, nr - 1, nr - 1);    
    Box bx(low, high);
    
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    Array<Geometry> geom;
    geom.resize(1);
    
    // level 0 describes physical domain
    geom[0].define(bx, &domain, 0, bc);
    
    Array<BoxArray> ba;
    
    ba.resize(1);
    
    ba[0].define(bx);
    ba[0].maxSize(size);
    
    Array<DistributionMapping> dmap;
    dmap.resize(1);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    container_t rhs(1);
    rhs.set(0, new MultiFab(ba[0], 1, 0));
    rhs[0].setVal(1.0);
    
    /* 
     * Set up file access property list with parallel I/O access
     */
     hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
     MPI_Comm comm = Ippl::getComm();
     MPI_Info info  = MPI_INFO_NULL;
     H5Pset_fapl_mpio(plist_id, comm, info);
    
    hsize_t      dims[3]  = {0, 0, 0};
    hsize_t      maxdims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED};
    hid_t dataspace = H5Screate_simple (3 /*number of dimensions*/, dims, maxdims); 
    hid_t file = H5Fcreate (FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    
    hsize_t      chunk_dims[3] = {1, 1, 1};
    hid_t prop = H5Pcreate (H5P_DATASET_CREATE);
    herr_t status = H5Pset_chunk (prop, 3 /*number of dimensions*/, chunk_dims);
    
    hid_t dataset = H5Dcreate2 (file, DATASETNAME, H5T_NATIVE_INT, dataspace,
                                H5P_DEFAULT, prop, H5P_DEFAULT);
    status = H5Pclose (prop);
    
    int count = 0;
    
    // minimal box extent per level
    Box minBox = rhs[0].boxArray().minimalBox();
    IntVect minsm = minBox.smallEnd();
    IntVect minbg = minBox.bigEnd();
    int xmax = minbg[0] - minsm[0] + 1;
    int ymax = minbg[1] - minsm[1] + 1;
    int zmax = minbg[2] - minsm[2] + 1;
    
    
    std::cout << "Minimal box: " << rhs[0].boxArray().minimalBox() << std::endl;
    
    for (MFIter mfi(rhs[0]); mfi.isValid(); ++mfi) {
        const FArrayBox& field = (rhs[0])[mfi];
        const Box& bx = field.box();
        IntVect sm = bx.smallEnd();
        IntVect bg = bx.bigEnd();
        
        hsize_t      dimsext[3] = {
            hsize_t(bg[0] - sm[0] + 1),
            hsize_t(bg[1] - sm[1] + 1),
            hsize_t(bg[2] - sm[2] + 1)
        };
        
        std::vector< int > dataext(dimsext[0] * dimsext[1] * dimsext[2]);
        int i = 0;
        for(IntVect p(sm); p <= bg; bx.next(p)) {
            for(int k(0); k < 1/*num_comp*/; ++k) {
                std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
                dataext[i++/*p[0] % dataext.size()][p[1] % dataext[0].size()][p[2] % dataext[1].size()*/] = field(p,k+0/*comp*/);
            }
        }
        
//         xmax = ( xmax > bg[0] ) ? xmax : bg[0] + 1;
//         ymax = ( ymax > bg[1] ) ? ymax : bg[1] + 1;
//         zmax = ( zmax > bg[2] ) ? zmax : bg[2] + 1;
        
//         std::cout << "max: " << xmax << " " << ymax << " " << zmax << std::endl;
        
        hsize_t size[3] = {
            /*4, */ hsize_t(xmax),
            /*4, */ hsize_t(ymax),
            /*4  */ hsize_t(zmax)
        };
        
        status = H5Dset_extent (dataset, size);
        
        hsize_t offset[3] = {
            hsize_t(sm[0]),
            hsize_t(sm[1]),
            hsize_t(sm[2])
        };
        
        hid_t filespace = H5Dget_space(dataset);
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     dimsext, NULL);
        hid_t memspace = H5Screate_simple(3 /*number of dimensions*/, dimsext, NULL); 
        
        status = H5Dwrite (dataset, H5T_NATIVE_INT, memspace, filespace,
                       H5P_DEFAULT, &dataext[0]);
            
        status = H5Sclose (memspace);
        
//         std::cout << "Status memspace: " << status << std::endl;
        
        status = H5Sclose (filespace);
        
//         std::cout << "Status filespace: " << status << std::endl;
    }
    
    rhs.clear();
    
    status = H5Sclose (dataspace);
//     std::cout << "Status dataspace: " << status << std::endl;
    status = H5Dclose (dataset);
//     std::cout << "Status dataset: " << status << std::endl;
    status = H5Pclose(plist_id);
//     std::cout << "Status plist_id: " << status << std::endl;
    status = H5Fclose (file);
    
    std::cout << "Status file: " << status << std::endl;
    
    return 0;
}     
