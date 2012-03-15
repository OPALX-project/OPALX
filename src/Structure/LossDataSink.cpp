#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <Structure/LossDataSink.h>


void LossDataSink::openH5(string element) {
    
    if ((!hdf5FileIsOpen_m) && doHdf5Save_m) {
	// open h5 file
	fn_m = element + string(".h5");
	/// Open H5 file. Check that it opens correctly.
#ifdef PARALLEL_IO
	H5file_m=H5PartOpenFileParallel(fn_m.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
	H5file_m=H5PartOpenFile(fn_m.c_str(),H5PART_WRITE);
#endif
	
	if(!H5file_m) {
	    ERRORMSG("h5 file open failed: exiting!" << endl);
	    exit(0);
	}
	
	writeH5FileAttributes() ;
	
	hdf5FileIsOpen_m = true;
	INFOMSG("H5 file of proble " << element << " is open" << endl); 
    }	
}
               
void LossDataSink::writeH5FileAttributes() {

    // Write file attributes to describe phase space to H5 file.
    stringstream OPAL_version;
    OPAL_version << PACKAGE_STRING << " rev. " << SVN_VERSION;
    H5PartWriteFileAttribString(H5file_m,"OPAL_version",OPAL_version.str().c_str());

    H5PartWriteFileAttribString(H5file_m,"tUnit","s");
    H5PartWriteFileAttribString(H5file_m,"xUnit","mm");
    H5PartWriteFileAttribString(H5file_m,"yUnit","mm");
    H5PartWriteFileAttribString(H5file_m,"zUnit","mm");
    H5PartWriteFileAttribString(H5file_m,"pxUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"pyUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"pzUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"idUnit","1");
    H5PartWriteFileAttribString(H5file_m,"SPOSUnit","mm");
    H5PartWriteFileAttribString(H5file_m,"TIMEUnit","s");

    H5PartWriteFileAttribString(H5file_m,"mpart","GeV");
    H5PartWriteFileAttribString(H5file_m,"qi","C");

    stringstream inputFileContent;
    hsize_t ContentLength;
    hsize_t write_length;
    hsize_t length;
    hsize_t start = 0;
    hsize_t dmax = H5S_UNLIMITED;

    herr_t herr;

    hid_t group_id;
    hid_t shape;
    hid_t dataset_id;
    hid_t diskshape;
    hid_t memshape;

    char group_name[] = "INPUT";
    char dataset_name[] = "InputFile";

    char *FileContent = NULL;

    if ( H5file_m->timegroup >= 0 ) {
        herr = H5Gclose ( H5file_m->timegroup );
        H5file_m->timegroup = -1;
    }

    if (Ippl::myNode() == 0) {
        struct stat st;
        off_t fsize;
        if (stat(OPAL.getInputFn().c_str(), &st) == 0) {
            fsize = st.st_size;
        }
        ContentLength = fsize / sizeof(char);
        FileContent = new char[ContentLength];

        filebuf inputFileBuffer;
        inputFileBuffer.open(OPAL.getInputFn().c_str(), ios::in);
        istream inputFile(&inputFileBuffer);
        
        inputFile.get(FileContent, ContentLength, '\0');

        inputFileBuffer.close();
        write_length = ContentLength;

    } else {
        FileContent = new char[1];
        //        FileContent[0] = '.';
        write_length = 0;
    }

    MPI_Bcast ( &ContentLength, 
                1, 
                MPI_LONG_LONG_INT, 
                0, 
                MPI_COMM_WORLD);
    
    herr = H5Gget_objinfo( H5file_m->file, group_name, 1, NULL );
    if (herr >= 0) { // there exists a group 'INPUT'
        delete[] FileContent;
        return;
    }
    
    group_id = H5Gcreate ( H5file_m->file, group_name, 0 );

    shape = H5Screate_simple(1, &ContentLength, &ContentLength);
    dataset_id = H5Dcreate( group_id,
                            dataset_name,
                            H5T_NATIVE_CHAR,
                            shape,
                            H5P_DEFAULT );
    H5Sclose ( shape );

    diskshape = H5Dget_space ( dataset_id );
    H5Sselect_hyperslab ( diskshape, 
                          H5S_SELECT_SET, 
                          &start, 
                          NULL, 
                          &write_length, 
                          NULL );

    memshape = H5Screate_simple(1, &write_length, &dmax);

    herr = H5Dwrite ( dataset_id,
                      H5T_NATIVE_CHAR,
                      memshape,
                      diskshape,
                      H5file_m->xfer_prop,
                      FileContent );

    H5Sclose ( memshape );
    H5Dclose ( dataset_id );
    H5Sclose ( diskshape );
    H5Gclose ( group_id );

    delete[] FileContent;

    INFOMSG("H5 Header of is written" << endl); 

}

void LossDataSink::addParticle(Vektor<double,3> x, Vektor<double,3> p, size_t  id) {
    x_m.push_back(x(0));
    y_m.push_back(x(1));
    z_m.push_back(x(2));
    px_m.push_back(p(0));
    py_m.push_back(p(1));
    pz_m.push_back(p(2));
    id_m.push_back(id);
}

void LossDataSink::save(string element)	{

    element_m=element;

    if (hdf5FileIsOpen_m) {
	size_t nLoc = x_m.size();
	void *varray = malloc(nLoc*sizeof(double));
	double *farray = (double*)varray;
	h5part_int64_t *larray = (h5part_int64_t *)varray;
	
	///Get the particle decomposition from all the compute nodes.
	size_t *locN = (size_t *) malloc(Ippl::getNodes()*sizeof(size_t));
	size_t  *globN = (size_t*) malloc(Ippl::getNodes()*sizeof(size_t));

	for(int i=0; i<Ippl::getNodes(); i++) {
	    globN[i] = locN[i]=0;
	}
	locN[Ippl::myNode()] = nLoc;
	reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());

	/// Set current record/time step.
	H5PartSetStep(H5file_m,H5call_m);
	H5PartSetNumParticles(H5file_m,nLoc);	

	// Write all data
	for (size_t i=0; i<nLoc;i++)
	    farray[i] = x_m[i];
	H5PartWriteDataFloat64(H5file_m,"x",farray);

	for (size_t i=0; i<nLoc;i++)
	    farray[i] = y_m[i];
	H5PartWriteDataFloat64(H5file_m,"y",farray);

	for (size_t i=0; i<nLoc;i++)
	    farray[i] = z_m[i];
	H5PartWriteDataFloat64(H5file_m,"z",farray);

	for (size_t i=0; i<nLoc;i++)
	    farray[i] = px_m[i]; 
	H5PartWriteDataFloat64(H5file_m,"px",farray);

	for (size_t i=0; i<nLoc;i++)
	    farray[i] = py_m[i];
	H5PartWriteDataFloat64(H5file_m,"py",farray);

	for (size_t i=0; i<nLoc;i++)
	    farray[i] = pz_m[i];
	H5PartWriteDataFloat64(H5file_m,"pz",farray);

	/// Write particle id numbers.
	for (size_t i = 0; i < nLoc; i++)
	    larray[i] =  id_m[i];
	H5PartWriteDataInt64(H5file_m,"id",larray);
	
	H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

	/// Step record/time step index.
	H5call_m++;

	if(varray)
	    free(varray);
    }
    else {
	int tag = Ippl::Comm->next_tag(IPPL_APP_TAG3, IPPL_APP_CYCLE);
	if( Ippl::Comm->myNode() == 0 ) { 
	    const unsigned partCount = x_m.size();
	    for(unsigned i=0; i<partCount; i++) {
		os_m << element_m   << "   ";
		os_m << x_m[i] << "   ";
		os_m << y_m[i] << "   ";
		os_m << z_m[i] << "   ";
		os_m << px_m[i] << "   ";
		os_m << py_m[i] << "   ";
		os_m << pz_m[i] << "   ";
		os_m << id_m[i]   << "   "<<endl;

	    }
	    int notReceived =  Ippl::getNodes() - 1;  
	    while (notReceived > 0) {
		unsigned dataBlocks=0;
		int node = COMM_ANY_NODE;
		Message* rmsg =  Ippl::Comm->receive_block(node, tag);
		if (rmsg == 0) {
		    ERRORMSG("LossDataSink: Could not receive from client nodes output." << endl);
		}
		notReceived--;
		rmsg->get(&dataBlocks);
		for (unsigned i=0; i < dataBlocks; i++) {
		    long id;
		    double rx,ry,rz,px,py,pz;			  
		    rmsg->get(&id);
		    rmsg->get(&rx);
		    rmsg->get(&ry);
		    rmsg->get(&rz);
		    rmsg->get(&px);
		    rmsg->get(&py);
		    rmsg->get(&pz);
		    os_m << element_m <<"   ";	
		    os_m << rx << "   ";
		    os_m << ry << "   ";
		    os_m << rz << "   ";
		    os_m << px << "   ";
		    os_m << py << "   ";
		    os_m << pz << "   ";
		    os_m << id <<" "<< endl;

		}
		delete rmsg;
	    }
	}
	else {
	    Message* smsg = new Message();
	    const unsigned msgsize = x_m.size();
	    smsg->put(msgsize);
	    for(unsigned i=0; i<msgsize; i++) {
		smsg->put(id_m[i]);
		smsg->put(x_m[i]);
		smsg->put(y_m[i]);
		smsg->put(z_m[i]);
		smsg->put(px_m[i]);
		smsg->put(py_m[i]);
		smsg->put(pz_m[i]);
	    }
	    bool res = Ippl::Comm->send(smsg, 0, tag);
	    if (! res) 
		ERRORMSG("LossDataSink Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
	}
    }
    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
}
