#include <hdf5.h>
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>


/**
 *\file DarkcurrentAlgo.cpp
 */

/**
 *\typedef Vector_t is an Opal Vektor<double,3> type.
 */
typedef Vektor<double,3> Vector_t;
/**
 * OPAL BoundaryGeometry Class Test Program.
 * 
 * Compile command:
 * CXX=mpicxx make.
 * Run command:
 * mpirun -np 1 ./DarkCurrentAlgo --commlib mpi.
 * 
 * At present,the BoundaryGeometry class read 
 * in .h5 created by Heronion. From .h5 file,
 * we store all points of the solid geometry
 * model and all surface triangles to obtain
 * the boundary geometry. This class provide
 * a simple stratergy to determine whether a
 * particle hits the boundary and to count the 
 * num of particles dumped in each surface 
 * triangle.
 */
class BoundaryGeometry
{
public:
    /**
     * Exemplar Constructor 
     */

    BoundaryGeometry() { 
	string dir("vtk");
	mkdir (dir.c_str(), O_CREAT);
	chmod (dir.c_str(), S_IRWXU);

	TMain_m = IpplTimings::getTimer("Main");	
	TPreProc_m = IpplTimings::getTimer("Pre Processing");
	TRayTrace_m = IpplTimings::getTimer("Ray tracing");
	TStep_m = IpplTimings::getTimer("Push particles");

	IpplTimings::startTimer(TMain_m);

    }
    /**
     * Destructor.
     *
     * Delete the previous defined member arrays 
     */
    ~BoundaryGeometry() { 
	if (allbfaces_m)
	    delete allbfaces_m;
	if (bfaces_idx_m)
	    delete bfaces_idx_m;
	if (TriPartloss_m)
	    delete TriPartloss_m;
	if (Tribarycent_m)
	    delete Tribarycent_m;
	if (TriPartlossZ_m)
	    delete TriPartlossZ_m;

	IpplTimings::stopTimer(TMain_m);
	IpplTimings::print();
	IpplTimings::print(string("timing.dat"));

    }
    /**
     * Initialize the BoundaryGeometry.
     *
     * @param fn specifies the name of .h5 file contains the geometry
     * 
     */

    void initialize(string fn) {

	IpplTimings::startTimer(TPreProc_m);
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS); //Property list identifier
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	// Create a new file collectively and release property list identifier.
	hid_t file_id = H5Fopen(fn.c_str(), H5F_ACC_RDONLY, plist_id);
	assert(file_id >= 0);
	H5Pclose(plist_id);

	///////////////////////////////////////////// 
	//   Read dataset "surface" from .h5 file
	////////////////////////////////////////////

	hsize_t dimsf[2];//dataset dimensions 
	herr_t  status;

	hid_t dset_id = H5Dopen(file_id, "/surface");
	assert(dset_id >= 0);
	// Create the dataspace for the dataset.
	hid_t x = H5Dget_space(dset_id);
	H5Sget_simple_extent_dims(x, dimsf, NULL);
	hid_t filespace = H5Screate_simple(2, dimsf, NULL);

	numbfaces_global_m = dimsf[0];
	INFOMSG("total_num= "<< numbfaces_global_m<< endl);

	allbfaces_m = new int [numbfaces_global_m * 4];

	TriPartloss_m = new int [numbfaces_global_m];

	bfaces_idx_m = new int [numbfaces_global_m];

	Tribarycent_m = new Vector_t[numbfaces_global_m];
	// Create property list for collective dataset write.
	hid_t plist_id2 = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id2, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace, plist_id2,allbfaces_m);
	assert(status >= 0);
	// Store local boundary faces, discard others
	nof_sym_m = 0; //Count number of symmetry planes.
	int const nogeosym_flag = 0; //heronion outputs a index, case 0 indicates no symmetry plane, 1 for (x,y) sym, 2 for (y,z), 3 for (z,x),none 0 means not a surface.
	for (int i = 0; i < numbfaces_global_m; i ++) {
	    bfaces_idx_m[i]=i;
	    if (allbfaces_m[4 * i] > nogeosym_flag) {
		nof_sym_m += 1;
		if (i<numbfaces_global_m-1){
		    for (int j = 0; j<numbfaces_global_m-i; j ++) {
			allbfaces_m[4*(i+j)]= allbfaces_m[4*(i+j+1)];
			allbfaces_m[4*(i+j)+1]= allbfaces_m[4*(i+j+1)+1];
			allbfaces_m[4*(i+j)+2]= allbfaces_m[4*(i+j+1)+2];
			allbfaces_m[4*(i+j)+3]= allbfaces_m[4*(i+j+1)+3];
		    }
		}
		else  
		    numbfaces_global_m=numbfaces_global_m-1;
	    }
	}
	H5Dclose(dset_id);
	H5Sclose(filespace);
	////////////////////////////////
	//  Also read dataset "coords"   
	////////////////////////////////
	hsize_t dimsf_c[2];
	herr_t status_c;
	hid_t dset_id_c = H5Dopen(file_id, "/coords");
	assert(dset_id_c >= 0);

	// Create the dataspace for the dataset.
	hid_t cp_space = H5Dget_space(dset_id_c);
	H5Sget_simple_extent_dims(cp_space, dimsf_c, NULL);
	hid_t filespace_c = H5Screate_simple(2, dimsf_c, NULL);

	numpoints_global_m = dimsf_c[0];
	double* point_coords = new double[3 * numpoints_global_m];
	hid_t plist_id3 = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id3, H5FD_MPIO_COLLECTIVE);
	status_c = H5Dread(dset_id_c, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_c, plist_id3, point_coords);
	assert(status_c >= 0);
	Vector_t geo3d_tmp;
	for (int i = 0;i<numpoints_global_m;i++) {
	    geo3d_tmp[0]=point_coords[3*i];
	    geo3d_tmp[1]=point_coords[3*i+1];
	    geo3d_tmp[2]=point_coords[3*i+2];

	    geo3Dcoords_m.push_back(geo3d_tmp);

	}
	INFOMSG("Points_tot_num= "<<numpoints_global_m<<" "<<geo3Dcoords_m.size()<< endl);

	for (int i = 0; i < numbfaces_global_m; i++) {

	    Tribarycent_m[i] = (geo3Dcoords_m[allbfaces_m[4*i+1]]+geo3Dcoords_m[allbfaces_m[4*i+2]]+geo3Dcoords_m[allbfaces_m[4*i+3]])/3.0 ;

	}
	if (point_coords)
	    delete point_coords;
	// Close HDF5 stuff
	H5Dclose(dset_id_c);
	H5Sclose(filespace_c);
	INFOMSG("Points_tot_num= "<<numpoints_global_m<<" "<<geo3Dcoords_m.size()<< endl);
	mincoords_m=getMinExtend();

	maxcoords_m=getMaxExtend();

	// fixme

	len_m = maxcoords_m - mincoords_m;  // have to fix !

	nr_m(0) = (int)floor(len_m(0)/getMaxDimenssion());
	nr_m(1) = (int)floor(len_m(1)/getMaxDimenssion()); 
	nr_m(2) = (int)floor(len_m(2)/getMaxDimenssion());

	hr_m = len_m / nr_m;
	INFOMSG("Points_tot_num= "<<numpoints_global_m<<" "<<geo3Dcoords_m.size()<< endl);
	TriPartlossZ_m = new int [nr_m(2)];


	INFOMSG("maxExtend= " << getMaxExtend() << " minExtend= " << getMinExtend()  << endl);
	INFOMSG("len= " << len_m << endl);
	INFOMSG("nr= " << nr_m << endl);
	INFOMSG("hr= " << hr_m << endl);
    }
    /**
     * Write geometry points and surface triangles to vtk file
     *
     * @param fn specifies the name of vtk file contains the geometry
     * 
     */

    void writeGeomToVtk(string fn) {
	std::ofstream of;
	of.open(fn.c_str());
	assert(of.is_open());
	of.precision (6);
	of << "# vtk DataFile Version 2.0" << endl;                      // first line of a vtk file
	of << "generated using BoundaryGeometry::writeToVtk" << endl;   // header
	of << "ASCII" << endl << endl;                                   // file format
	of << "DATASET UNSTRUCTURED_GRID" << endl;                 // dataset structrue
	of << "POINTS " << numpoints_global_m << " float" << endl;   // data num and type
	for (int i = 0; i < numpoints_global_m ; i ++)
	    of << geo3Dcoords_m[i](0) << " " << geo3Dcoords_m[i](1) << " " << geo3Dcoords_m[i](2) << endl;
	of << endl;

	of << "CELLS " << numbfaces_global_m << " " << 4*numbfaces_global_m << endl;
	for ( int i = 0; i < numbfaces_global_m ; i ++)
	    of << "3 " << allbfaces_m[4*i+1] << " " << allbfaces_m[4*i+2] << " " << allbfaces_m[4*i+3] << endl;
	of << "CELL_TYPES " << numbfaces_global_m << endl;
	for ( int i = 0; i < numbfaces_global_m ; i ++)
	    of << "5" << endl;
	of << "CELL_DATA " << numbfaces_global_m << endl;
	of << "SCALARS " << "cell_attribute_data" << " float " << "1" << endl;
	of << "LOOKUP_TABLE " << "default" << endl ;
	for (int i = 0; i < numbfaces_global_m ; i ++)
	    of <<(float)(i) << endl;
	of << endl;
	of << "COLOR_SCALARS " << "TriangleColor "<<4 << endl ;
	for (int i = 0; i < numbfaces_global_m ; i ++){
	    //     if(TriPartloss_m[i]!=0)
	    //	of <<"1.0"<<" 1.0 "<< endl;
	    //   else
	    of <<"0"<<" 1.0 "<<"1.0"<<" 1.0"<< endl;
	}

	of << endl;

    }
    /**
     * Write particle loss data to an ASCII fille for histogram
     * @param fn specifies the name of ASCII file 
     */

    void writePartlossToASCII(string fn) {
	std::ofstream of;
	of.open(fn.c_str());
	assert(of.is_open());
	of.precision (6);

	for (int i = 0; i < nr_m(2) ; i++ ){
	    TriPartlossZ_m[i] = 0;
	    for (int j = 0; j <numbfaces_global_m; j++ ){
		if(((mincoords_m[2]+hr_m(2)*i)<Tribarycent_m[j](2))&&(Tribarycent_m[j](2)<(hr_m(2)*i+hr_m(2)+mincoords_m[2]))){
		    TriPartlossZ_m[i]+=TriPartloss_m[j];
		}
	    }

	}
	of << "Delta_Z" <<std::setw(15)<< "Zcoordinates" <<std::setw(15)<<"numParticles" << endl ; 
	for (int i = 0; i < nr_m(2) ; i++ ){

	    of <<hr_m(2)<<std::setw(15)<<mincoords_m[2]+hr_m(2)*i<<std::setw(15)<<TriPartlossZ_m[i]<<endl;
	}

    }
    /**
     * Write particles to vtk file for visualization
     * @param fn specifies the name of vtk file 
     */

    void writePartToVtk(string fn) {
	std::ofstream of;
	of.open(fn.c_str());
	assert(of.is_open());
	of.precision (6);

	const int numParticles = partsr_m.size();

	of  << setprecision(5)
	    << "# vtk DataFile Version 2.0" << endl                
	    << "unstructured grid and vector field on the nodes" << endl    
	    << "ASCII" << endl                
	    << "DATASET UNSTRUCTURED_GRID" << endl                
	    << "POINTS " << numParticles << " float" << endl;            

	// Particle positions            
	for (int i=0; i < numParticles; i++)
	    // for (int i=0; i < partsr_m.size(); i++)
	    of << partsr_m[i](0) << "  " << partsr_m[i](1) << "  " << partsr_m[i](2) << endl; 

	of << endl;            // defining VTK_poly_vertex            

	of << "CELLS " << numParticles << " " << 2*numParticles << endl;            
	for (int i=0; i < numParticles; i++)                

	    of << "1 " << i << endl;           

	of << endl;            
	// defining Cell_types            
	of << "CELL_TYPES " << numParticles << endl;            
	for (int i=0; i < numParticles; i++)                

	    of << "2" << endl;            
    }

    /**
     * Only for code debuging
     */

    Vector_t getCoordinate() {
	Inform msg("DarkCurrentAlgo ");  
	double x = 0.0,y=0.0,z=0.0;

	msg << "Please specify x-coordinate " << endl;
	cin >> x;
	msg << "Please specify y-coordinate " << endl;
	cin >> y;
	msg << "Please specify z-coordinate " << endl;
	cin >> z;

	Vector_t xv(x,y,z);

	msg << " ---> Read in " << xv << endl;
	return xv;
    }
    /**
     * Used to determine whether a particle hits boundary.
     * @param x stands for the coordinates of the testing particle
     */
    // Not sufficient for particle hitting the triangle area.
    bool isInGeometry(Vector_t x) {
	return boundary_ids_m.find(f(x)) != boundary_ids_m.end();
    }

    /**
     * Calculate the maximum of coordinates of geometry,i.e the maximum of X,Y,Z
     */
    Vector_t getMaxExtend() {
	const Vector_t x = *max_element(geo3Dcoords_m.begin(),geo3Dcoords_m.end(), myLessx);
	const Vector_t y = *max_element(geo3Dcoords_m.begin(),geo3Dcoords_m.end(), myLessy);	
	const Vector_t z = *max_element(geo3Dcoords_m.begin(),geo3Dcoords_m.end(), myLessz);
	return Vector_t(x(0),y(1),z(2));
    }
    /**
     * Calculate the minimum of coordinates of geometry,i.e the minimum of X,Y,Z
     */
    Vector_t getMinExtend() {
	const Vector_t x = *min_element(geo3Dcoords_m.begin(),geo3Dcoords_m.end(), myLessx);
	const Vector_t y = *min_element(geo3Dcoords_m.begin(),geo3Dcoords_m.end(), myLessy);	
	const Vector_t z = *min_element(geo3Dcoords_m.begin(),geo3Dcoords_m.end(), myLessz);
	return Vector_t(x(0),y(1),z(2));
    }
    /**
     * Calculate the maximum dimention of triangles. This value will be used to define the cubic box size 
     */
    double getMaxDimenssion() {
	double maximum=0,min=0.1;
	for (int i = 0; i<numbfaces_global_m; i++) {
	    maximum = maximum>=Trianglelen(i)?maximum:Trianglelen(i);
	    min=min<Trianglelen(i)?min:Trianglelen(i);}
	INFOMSG("size of Triangle " << maximum <<" "<<min<< endl);
	return maximum ;
    }
    /**
     * Initialize some particles near the surface with inward momenta.
     */
    void createParticles() {
	INFOMSG("creat Particle called"<< endl);
	vector<Vector_t>::iterator myit;
	for (int i = 0;i< numbfaces_global_m;i++){
	    Vector_t  temp = Tribarycent_m[i];
	    partsr_m.push_back(temp);
	    temp =TriNormal_m[i]*hr_m[0];
	    partsp_m.push_back(temp);
	}
	assert(!partsp_m.empty());
	INFOMSG("Initial_Particle_num= "<<partsr_m.size()<< endl);
    }
    /**
     * Make the boundary set by using triangle vertex .
     */
    void makeBoundaryIndexSet() {
   
	for(int i = 0; i< numbfaces_global_m; i++){
	    vector<Vector_t> discreteTri;
	    vector<Vector_t>::iterator disIt;
	    for (int j = 0; j<20;j++) {
		discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]]+0.05*(j+1)*(0.5*(geo3Dcoords_m[allbfaces_m[4*i+2]]+geo3Dcoords_m[allbfaces_m[4*i+3]])-geo3Dcoords_m[allbfaces_m[4*i+1]]));
		discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]]+0.05*(j+1)*(0.5*(geo3Dcoords_m[allbfaces_m[4*i+3]]+geo3Dcoords_m[allbfaces_m[4*i+1]])-geo3Dcoords_m[allbfaces_m[4*i+2]]));
		discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]]+0.05*(j+1)*(0.5*(geo3Dcoords_m[allbfaces_m[4*i+1]]+geo3Dcoords_m[allbfaces_m[4*i+2]])-geo3Dcoords_m[allbfaces_m[4*i+3]]));
	    }
	    discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+1]]);
	    discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+2]]);
	    discreteTri.push_back(geo3Dcoords_m[allbfaces_m[4*i+3]]);
	    discreteTri.push_back(Tribarycent_m[i]);
	    for (disIt=discreteTri.begin();disIt!=discreteTri.end();disIt++) {
		size_t id = f(*disIt);
		assert(id > 0);
		boundary_ids_m.insert(id);
		vector<size_t> tmp;
		map< size_t,vector<size_t> >::iterator It;
		It =  CubicLookupTable_m.find(id);
		if (It==CubicLookupTable_m.end()) {
		    tmp.push_back(i);
		    CubicLookupTable_m.insert( pair <size_t,vector<size_t> > (id, tmp));
		}
		else 
		    (*It).second.push_back(i);
	    }
     
	}
	INFOMSG("size of boundary index set " << boundary_ids_m.size() <<" "<<b_cubcentre_m.size()<< endl);
	INFOMSG("size of all index " << nr_m(0)*nr_m(1)*nr_m(2) << endl);
	INFOMSG("size of geo3Dcoords " << geo3Dcoords_m.size() << endl);
	IpplTimings::stopTimer(TPreProc_m);
    }
  
    /**
     * Make the boundary set by using triangle barycentric point.
     */
    void makeTriangleIndexSet() {

	vector<Vector_t>::iterator myit;

	for (int i = 0;i< numbfaces_global_m; i++){
	    size_t id = f(Tribarycent_m[i]);
	    assert(id > 0);
	    boundary_ids_m.insert(id);        
	}
	INFOMSG("size of boundary index set " << boundary_ids_m.size() << endl);
	INFOMSG("size of geo3Dcoords " << geo3Dcoords_m.size() << endl);
    }
    /**
     * Test program for particle movement. 
     *
     * @param n define the tracking time of particle movement.
     * @param dt define the time interval of particle movement integration.
     */

    void integPartcleMove(int n,double dt) {

	vector<Vector_t>::iterator r,v;
	Vector_t temp,tempold;
	vector<Vector_t> escape;
	double maxnt;
	IpplTimings::startTimer( TStep_m);
	for(int i = 0; i<n ; i++){
	    for ( r = partsr_m.begin(), v = partsp_m.begin(); (r != partsr_m.end()); ){

		Vector_t temp1 = *r;
		Vector_t temp = *r+(*v*dt);
		if (isInGeometry(temp1)){
		    size_t id = f(temp1);
		    map< size_t,vector<size_t> >::iterator It;
		    vector<size_t> ::iterator myIt;
		    It=CubicLookupTable_m.find(id);//Lookup table (stl map) contains the id of cubic box as key and triangle ids inside the id th cubic box as value
		    int flg =0;
		    for(myIt=(*It).second.begin();myIt!=(*It).second.end();myIt++ ) {//For each particle inside the cubic box, test if this particle have intersect with each triangle inside the same cubic box 
	    
			if (Dotproduct(*v,TriNormal_m[*myIt])<0.0) {
			    IpplTimings::startTimer(TRayTrace_m);
			    //out
			    Vector_t tmp=LineInsTri(temp,temp1,*myIt);
			    if((tmp!=temp1)){
				lostedPart_m.push_back(*r);
				partsr_m.erase(r);
				partsp_m.erase(v);
				flg++;
				break;
			    }
			    else {//In some accassion,  the boundary box lies outward the boundary.
				if (Dotproduct(temp-Tribarycent_m[(*myIt)],TriNormal_m[(*myIt)])<0.0) {
				    lostedPart_m.push_back(*r);
				    partsr_m.erase(r);
				    partsp_m.erase(v);
				    flg++;
				    break;
				}
		
			    }
			    IpplTimings::stopTimer(TRayTrace_m); 
			}
	    
		    }
		    if (flg ==0) {
	    
			r++;
			v++;
		    }
	 
		}
		else {
		    r++;
		    v++;
		}
	    }

	    IpplTimings::stopTimer( TStep_m);
	    INFOMSG("time"<<i<<"_Particle_num= "<<partsr_m.size()<< endl);
	    INFOMSG("Lost Particle "<<lostedPart_m.size()<< endl);

	    writeGeomToVtk(string("vtk/testGeometry-") + convert2Int(i+1) + string(".vtk"));
	    writePartToVtk(string("vtk/testParticle-") + convert2Int(i+1) + string(".vtk"));
	    writePartlossToASCII(string("vtk/testPartloss-") + convert2Int(i+1) + string(".dat"));
	    for ( r = partsr_m.begin(), v = partsp_m.begin(); (r != partsr_m.end())&&(v!=partsp_m.end());r++,v++ ){
		*r += (*v*dt);
	    }
	}

	Triangle_lossmax_m=0;
	for (int i = 0;i< numbfaces_global_m; i++) {
	    Triangle_lossmax_m =( Triangle_lossmax_m>=TriPartloss_m[i])?Triangle_lossmax_m:TriPartloss_m[i]; 
	}
	INFOMSG(" Triangle maximum particle loss "<<Triangle_lossmax_m<<" "<<partsr_m.size()<< endl);
    }
    /**
     * Calculate the number of intersects between a line segment and the geometry boundary,make sure x1 is outside the geometry.
     * @param
     */

    vector<Vector_t>  PartBoundaryInteNum(Vector_t x0, Vector_t x1) {

	using namespace std;

	vector<Vector_t> SegDiscrete;
	vector<Vector_t> Isp;
	vector<Vector_t>::iterator TriIscIt;

	vector<int> TriId;
	vector<int>::iterator TriIdIt;
	double len =sqrt((x0[0]-x1[0])*(x0[0]-x1[0])+(x0[1]-x1[1])*(x0[1]-x1[1])+(x0[2]-x1[2])*(x0[2]-x1[2]));

	int Seglen =(((int)floor(len/hr_m[0]))+1);
	int count = 0;
	for (int i = 0; i<Seglen ;i++) {
	    SegDiscrete.push_back(x0+hr_m[0]*i*(x1-x0)/len);
	}
	SegDiscrete.push_back(x1);
	vector<Vector_t>::iterator myit;
	for (myit = SegDiscrete.begin() ; myit != SegDiscrete.end() ; myit++) {
	    if (!isInGeometry(*myit)) {//segment do not cross the boundary cubic;
		continue;
	    }
	    else {//segment cross the boundary cubic;
		int id = f(*myit);
		for (int i = 0;i< numbfaces_global_m; i++) {
		    if (((f(geo3Dcoords_m[allbfaces_m[4*i+1]])==id)||(f(geo3Dcoords_m[allbfaces_m[4*i+2]])==id)||(f(geo3Dcoords_m[allbfaces_m[4*i+3]])==id))){
			Vector_t tmp =  LineInsTri(x0, x1, i);
			TriIdIt = std::find(TriId.begin(),TriId.end(),i);
			if (((tmp!=x1)&&((TriIdIt==TriId.end())||TriId.size()==0))) {
			    TriIscIt = std::find(Isp.begin(),Isp.end(),tmp);
			    if((TriIscIt==Isp.end())||(Isp.size()==0)) {
				TriId.push_back(i);
				Isp.push_back(tmp);
				++count;
			    }
			}
			else 
			    continue;
		    }
		}
	    }
	}
	if (count==0)
	    Isp.push_back(x1);
	return Isp;
    }



    /**
     * Determine if a point is outside (return 0), inside (return 1) or just on (return 2) the boundary.
     * @param x stands for coordinates of the test point.
     * The basic idea is if a line segment starting from the test point has odd intersects with a closed boundary, 
     * then the test point is inside the geometry; if the intersects have even number, then the test points is outside the geometry;
     * if the test point is amoung the intersects, then the test point is just on the boundary. Make sure the end point of the line 
     * segment is outside the geometry boundary.
     */   
    int isInsideGeometry(Vector_t x) {
	Vector_t x0=x;
	Vector_t x1;
	x1[0] = x0[0];
	x1[1] = x0[1];
	x1[2] = maxcoords_m[2]*(1+IpplRandom());//Random number could avoid some specific situation,like line parallel to boundary......
	// x1 could be any point outside the boundary like :x1[0] = hr_m[0]+maxcoords_m[0]; x1[1] = hr_m[1]+maxcoords_m[1]  x1[2] = x0[2] ;
	vector<Vector_t> IntesecNum= PartBoundaryInteNum( x0, x1);
	vector<Vector_t>::iterator InIt;
	for (InIt = IntesecNum.begin();InIt != IntesecNum.end();InIt++){
	    //     INFOMSG("IntesecNum "<< *InIt<<IntesecNum[0]-x0<<" x0= "<<x0 <<" Begin "<<endl);
	}
	Vector_t temp = IntesecNum[0]-x0;
	// INFOMSG("temp "<< temp[0]<<" "<<temp[1] <<" "<<temp[2]<<endl);
	if (((temp[2])<1.0e-16)&&((temp[2])>-1.0e-16)&&((temp[1])<1.0e-16)&&((temp[1])>-1.0e-16)&&((temp[0])<1.0e-16)&&((temp[0])>-1.0e-16)) {
	    //    INFOMSG("IntesecNum[0] "<< IntesecNum[0]<<endl);
	    return 2;//x0 is just on the boundary;
	}
	else {
	    if (((IntesecNum.size()%2)==0)||(*IntesecNum.begin()==x1)) {
		return 0;//x0 is  outside the boundary; 
	    }
	    else
		return 1;//x0 is inside the boundary;
	}
    }
    //(we have to make sure that x0 will not jump outside in one time step and the initial particle will not outside the geometry)



    void testNp() {
	Vector_t temp,tp;
	makeTriNormal();
    }

private:
    
    IpplTimings::TimerRef TMain_m;
    IpplTimings::TimerRef TPreProc_m;
    IpplTimings::TimerRef TRayTrace_m;
    IpplTimings::TimerRef TGuessIsInGeom_m;
    IpplTimings::TimerRef TStep_m;



    map< size_t,vector<size_t> >  CubicLookupTable_m;          
    set<size_t> isOriented_m;
    map< size_t, vector<size_t> > triangleLookupTable_m;
    void orientAllTriangles(size_t caller, size_t idx) {
	// cout << "." << flush;
	orientTriangle(idx, caller);
	isOriented_m.insert(idx);
	vector<size_t> neighbours = findNeighbours(idx);
	for(int i=0; i<neighbours.size(); i++) { 
	    if(isOriented_m.find(neighbours[i]) == isOriented_m.end())
		orientAllTriangles(idx, neighbours[i]);
	}
    }


    vector<size_t> findNeighbours(size_t idx) {
	vector<size_t> ret,ret1;
	set<size_t>temp;
	vector<size_t>::iterator retIt;
	map< size_t, vector<size_t> >::iterator it;
	for(int i=1; i<=3; i++) {
	    it = triangleLookupTable_m.find(triangleCorner(idx, i));
	    ret.insert(ret.end(), (*it).second.begin(), (*it).second.end());
	}
	for (retIt=ret.begin();retIt!=ret.end();retIt++) {
	    if(*retIt!=idx){
		if (temp.find(*retIt)==temp.end()){
		    temp.insert(*retIt);
		}
		else
		    ret1.push_back(*retIt);
	    }
	}
	assert(ret1.size()!=0);
	return ret1;
    }
  
    size_t triangleCorner(size_t idx, size_t corner) const {
	return allbfaces_m[4*idx + corner];
    }
    vector<size_t> alignedT_m,NaliT_m;
    void orientTriangle(size_t idx, size_t caller) {
	vector<size_t> id,ic;
	for(int i=1; i<=3; i++) {
	    for(int j=1; j<=3; j++) {
		if( triangleCorner(idx, j) == triangleCorner(caller, i)) {
		    id.push_back(j);
		    ic.push_back(i);
		}
	    }
	}
	if((ic[1]-ic[0])==1){
	    int idtmp=id[1]-id[0];
	    if (idtmp == 1) {
		alignedT_m.push_back(idx);
		allbfaces_m[4*idx+id[0]]=allbfaces_m[4*caller+ic[1]];
		allbfaces_m[4*idx+id[1]]=allbfaces_m[4*caller+ic[0]];
	    }
	    if (idtmp == -1) {
		NaliT_m.push_back(idx);
	    }
	    if (idtmp == 2) {
		NaliT_m.push_back(idx);
	    }
	    if (idtmp == -2){
		alignedT_m.push_back(idx);
		allbfaces_m[4*idx+id[0]]=allbfaces_m[4*caller+ic[1]];
		allbfaces_m[4*idx+id[1]]=allbfaces_m[4*caller+ic[0]];
	    }
	}
	if((ic[1]-ic[0])==2){
	    int idtmp=id[1]-id[0];
	    if (idtmp==-1) {
		alignedT_m.push_back(idx);
		allbfaces_m[4*idx+id[0]]=allbfaces_m[4*caller+ic[1]];
		allbfaces_m[4*idx+id[1]]=allbfaces_m[4*caller+ic[0]];
	    }
	    if (idtmp==1) {
		NaliT_m.push_back(idx);
	    }
	    if (idtmp==-2) {
		NaliT_m.push_back(idx);
	    }
	    if (idtmp==2) {
		alignedT_m.push_back(idx);
		allbfaces_m[4*idx+id[0]]=allbfaces_m[4*caller+ic[1]];
		allbfaces_m[4*idx+id[1]]=allbfaces_m[4*caller+ic[0]];
	    }
	
	}
 
    }

    void makeTriNormal() {

	Vector_t t0,t1,t2;
	t0 = geo3Dcoords_m[allbfaces_m[1]];
	t1 = geo3Dcoords_m[allbfaces_m[2]];
	t2 = geo3Dcoords_m[allbfaces_m[3]];
	Vector_t u = t1-t0;
	Vector_t v = t2-t0;
	Vector_t n;
	n[0] = u[1]*v[2]-v[1]*u[2];
	n[1] = u[2]*v[0]-v[2]*u[0];
	n[2] = u[0]*v[1]-v[0]*u[1];
	double nomtmp=sqrt(n(0)*n(0)+n(1)*n(1)+n(2)*n(2));
	n=n/nomtmp;
	TriNormal_m.push_back(n);

	Vector_t mytemp;
	mytemp[0]= t0[0]+ TriNormal_m[0](0)*0.00000005;
	mytemp[1]= t0[1]+ TriNormal_m[0](1)*0.00000005;
	mytemp[2]= t0[2]+ TriNormal_m[0](2)*0.00000005;
	if(isInsideGeometry(mytemp)==0){
	    INFOMSG("Inside normal outward "<< endl);
	    if(((mytemp[0]-t0[0])*TriNormal_m[0](0)+(mytemp[1]-t0[1])*TriNormal_m[0](1)+(mytemp[2]-t0[2])*TriNormal_m[0](2))<0){
		// TriNormal_m[i]=TriNormal_m[i];
	    }
	    else {
		TriNormal_m[0] = -TriNormal_m[0];
		int temp = allbfaces_m[2];
		allbfaces_m[2] = allbfaces_m[3];
		allbfaces_m[3] = temp; 
	    }
	}
	else {
	    if(((mytemp[0]-t0[0])*TriNormal_m[0](0)+(mytemp[1]-t0[1])*TriNormal_m[0](1)+(mytemp[2]-t0[2])*TriNormal_m[0](2))<0){
		TriNormal_m[0] = -TriNormal_m[0];
		int temp = allbfaces_m[2];
		allbfaces_m[2] = allbfaces_m[3];
		allbfaces_m[3] = temp; 
	    }
	    else {
		TriNormal_m[0] = TriNormal_m[0];
	    }
	}
	for (size_t i=0; i<numbfaces_global_m; i++) {//for every triangle find adjacent triangles for each vertex ;
	    for(int j=1; j<=3; j++) {
		map< size_t, vector<size_t> >::iterator it;
		it = triangleLookupTable_m.find(allbfaces_m[4*i + j]);
		if(it == triangleLookupTable_m.end()) {
		    vector <size_t> tmp;
		    tmp.push_back(i);
		    triangleLookupTable_m.insert( pair<size_t, vector<size_t> > (allbfaces_m[4*i + j], tmp));
		} else
		    (*it).second.push_back(i);
	    }
	}
	int temp = allbfaces_m[5];
	allbfaces_m[5] = allbfaces_m[6];
	allbfaces_m[6] = temp;  
	isOriented_m.insert(0);
	vector<size_t> neighbours = findNeighbours(0);

	for(int i=0; i<neighbours.size(); i++) { 
	    if(isOriented_m.find(neighbours[i]) == isOriented_m.end())
		orientAllTriangles(0, neighbours[i]);
	}
	cout << endl;


	INFOMSG("oriented " << isOriented_m.size() << " triangles (of " << numbfaces_global_m << ") delta = " << numbfaces_global_m - isOriented_m.size() <<""<<alignedT_m.size()<<" "<<NaliT_m.size()<< endl);
            
	/*=====================================================*/

	std::ofstream of;
	of.open(string("vtk/testnTri.vtk").c_str());
	assert(of.is_open());
	of.precision (6);
	of << "# vtk DataFile Version 2.0" << endl;                      // first line of a vtk file
	of << "generated using BoundaryGeometry::writeToVtk" << endl;   // header
	of << "ASCII" << endl << endl;                                   // file format
	of << "DATASET UNSTRUCTURED_GRID" << endl;                 // dataset structrue
	of << "POINTS " << numpoints_global_m << " float" << endl;   // data num and type
	for (int i = 0; i < numpoints_global_m ; i ++)
	    of << geo3Dcoords_m[i](0) << " " << geo3Dcoords_m[i](1) << " " << geo3Dcoords_m[i](2) << endl;
	of << endl;
	of << "CELLS " << numbfaces_global_m << " " << 4*numbfaces_global_m << endl;
	for ( int i = 0; i < numbfaces_global_m ; i ++)
	    of << "3 " << allbfaces_m[4*i+1] << " " << allbfaces_m[4*i+2] << " " << allbfaces_m[4*i+3] << endl;
	of << "CELL_TYPES " << numbfaces_global_m << endl;
	for ( int i = 0; i < numbfaces_global_m ; i ++)
	    of << "5" << endl;
	of << "CELL_DATA " << numbfaces_global_m << endl;
	of << "SCALARS " << "cell_attribute_data" << " float " << "1" << endl;
	of << "LOOKUP_TABLE " << "default" << endl ;
	for (int i = 0; i < numbfaces_global_m ; i ++)
	    of <<(float)(i) << endl;
	of << endl;
	of << "COLOR_SCALARS " << "TriangleColor "<<4 << endl ;
	for (int i = 0; i < numbfaces_global_m ; i++){
	    if(isOriented_m.find(i) == isOriented_m.end()) 
		of <<"0.5"<<" 1.0 "<<"1.0"<<" 1.0"<< endl;
	    else 
		of <<"0.0"<<" 1.0 "<<"0.5"<<" 1.0"<< endl;
	}
	INFOMSG("write GeonTri  "<< endl);

	of << endl;

	/*====================================================*/

	for (int i = 1; i < numbfaces_global_m; i ++) {
	    Vector_t t0,t1,t2;
	    t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
	    t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
	    t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
	    Vector_t u = t1-t0;
	    Vector_t v = t2-t0;
	    Vector_t n;
	    n[0] = u[1]*v[2]-v[1]*u[2];
	    n[1] = u[2]*v[0]-v[2]*u[0];
	    n[2] = u[0]*v[1]-v[0]*u[1];
	    double nomtmp=sqrt(n(0)*n(0)+n(1)*n(1)+n(2)*n(2));
	    if (nomtmp !=0){
		n=n/nomtmp;
	    }
	    TriNormal_m.push_back(n);
	}


	INFOMSG("Outside make normal "<< endl);
	/*=======================================================*/
	std::ofstream mof;
	mof.open(string("vtk/testNormal.vtk").c_str());
	assert(mof.is_open());
	mof.precision (6);
	mof  << setprecision(5)
	     << "# vtk DataFile Version 2.0" << endl                
	     << "unstructured grid and vector field on the nodes" << endl    
	     << "ASCII" << endl                
	     << "DATASET UNSTRUCTURED_GRID" << endl                
	     << "POINTS " << (int)floor(numbfaces_global_m/1) << " float" << endl;            
	//	<< "POINTS " <<numbfaces_global_m << " float" << endl;            
	// Particle positions            
	for (int i=0; i <numbfaces_global_m; i++)
	    // for (int i=0; i < partsr_m.size(); i++)
	    if ((i%1)==0){
		mof << Tribarycent_m[i](0) << "  " << Tribarycent_m[i](1) << "  " <<Tribarycent_m[i](2) << endl; 
	    }

	mof << endl;            // defining VTK_poly_vertex            

	mof << "CELLS " <<(int)floor(numbfaces_global_m/1) << " " << 2*((int)floor(numbfaces_global_m/1)) << endl;            
	for (int i=0; i <(int)floor(numbfaces_global_m/1); i++) 
	    // mof << "CELLS " <<numbfaces_global_m << " " << 2*numbfaces_global_m << endl;            
	    // for (int i=0; i <numbfaces_global_m; i++)                               

	    mof << "1 " << i << endl;           

	mof << endl;            
	// defining Cell_types            
	mof << "CELL_TYPES " <<(int)floor(numbfaces_global_m/1) << endl;            
	for (int i=0; i <(int)floor(numbfaces_global_m/1); i++)                
	    // mof << "CELL_TYPES " <<numbfaces_global_m << endl;            
	    //for (int i=0; i <numbfaces_global_m; i++)              
	    mof << "2" << endl;   
	mof << "POINT_DATA " <<(int)floor(numbfaces_global_m/1) << endl;
	// mof << "POINT_DATA " <<numbfaces_global_m << endl;
	mof << "SCALARS " << "point_attribute_data" << " float " << "1" << endl;
	mof << "LOOKUP_TABLE " << "default" << endl ;
	for (int i = 0; i < (int)floor(numbfaces_global_m/1) ; i ++)
	    // for (int i = 0; i < numbfaces_global_m; i ++)
	    mof <<(float)(i) << endl;
	mof << endl;
	mof << "VECTORS " << "TriangleNormals "<<"float" << endl ;
	for (int i = 0; i < numbfaces_global_m ; i ++){
	    //     if(TriPartloss_m[i]!=0)
	    //	mof <<"1.0"<<" 1.0 "<< endl;
	    //   else
	    if ((i%1)==0){
		mof <<TriNormal_m[i](0)<<" "<<TriNormal_m[i](1)<<" "<<TriNormal_m[i](2)<< endl;
	    }
	}
	mof << endl;
	/*================================================================*/
	map<size_t,Vector_t>::iterator outIt;
	for(outIt =  b_cubcentre_m.begin();outIt !=b_cubcentre_m.end();outIt++){
	    // INFOMSG("Cubic center: "<<(*outIt).first<<" "<<(*outIt).second<<" "<<b_cubcentre_m.size()<< endl); 
	}
	INFOMSG("Size of Triangle normal vector: "<<TriNormal_m.size()<< endl);
	//  INFOMSG("Inside normal outward= "<<counter<< endl);
    }

    double Dotproduct(Vector_t x0, Vector_t x1) {
	return (x0[0]*x1[0]+x0[1]*x1[1]+x0[2]*x1[2]);
    }



    /** 
     * Map a 3D coordinate to a cubic box with a unique Id
     * 
     *  @param x is an OPAL Vector_t type, stands for the 3D coordinates need to be mapped.
     */
    size_t f(Vector_t x) {
	/**  
	     We know: lenght of the structure len_m
	     Mesh size: hr_m
	     Number of mesh points nr_m
	*/	  
	size_t ret = 0;
	size_t id_tx,id_ty,id_tz;
	id_tx =(size_t)floor((x[0]-mincoords_m[0])/hr_m[0]);
	id_ty =(size_t)floor((x[1]-mincoords_m[1])/hr_m[1]);
	id_tz =(size_t)floor((x[2]-mincoords_m[2])/hr_m[2]);
	ret = 1 + id_tz*nr_m[0]*nr_m[1]+id_ty*nr_m[0]+id_tx;
	return ret;
    }
    /**
     * Determine if a line segment  has intersects with a triangle.
     * Algorithm :http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
     * @param x0 and @param x1 define the line segment. 
     * @param i is the id of the triangle.
     */

    Vector_t LineInsTri(Vector_t x0, Vector_t x1, size_t i) {
	//       INFOMSG(" LineInsTri called "<< endl);
	Vector_t lseg = x1-x0;
	Vector_t t0,t1,t2;
	t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
	t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
	t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
	Vector_t u = t1-t0;
	Vector_t v = t2-t0;
	Vector_t lt = t0-x0;
	Vector_t n;
	n[0] = u[1]*v[2]-v[1]*u[2];
	n[1] = u[2]*v[0]-v[2]*u[0];
	n[2] = u[0]*v[1]-v[0]*u[1];

	if ((n[0]*lseg[0]+n[1]*lseg[1]+n[2]*lseg[2])==0) {
	    //          INFOMSG("Triangle parallell to line segment, return x1 "<< endl);
	    return x1;
	}
	else {

	    double rI = (n[0]*lt[0]+n[1]*lt[1]+n[2]*lt[2])/(n[0]*lseg[0]+n[1]*lseg[1]+n[2]*lseg[2]);
	    Vector_t ItSec = x0+rI*lseg;
	    Vector_t w = ItSec-t0;
	    if ((rI<0)||(rI>1)) {
		//	INFOMSG("Intesect is on the extended line, return x1 "<< endl);
		return x1;
	    }
	    else {
		double tmp1,tmp2,tmp3,tmp4,tmp5;
		tmp1 = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
		tmp2 = w[0]*v[0]+w[1]*v[1]+w[2]*v[2];
		tmp3 = u[0]*w[0]+u[1]*w[1]+u[2]*w[2];
		tmp4 = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
		tmp5 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
		double sI = (tmp1*tmp2-tmp5*tmp3)/(tmp1*tmp1-tmp4*tmp5);
		double tI = (tmp1*tmp3-tmp4*tmp2)/(tmp1*tmp1-tmp4*tmp5);
		if ((sI>=0)&&(tI>=0)&&((sI+tI)<=1)) {
		    return ItSec;
		}
		else {
		    //		  INFOMSG("Intesect is on the extended plane, return x1 "<< endl);
		    return x1;
		}
	    }
	}
    }
    struct Intersection{
	double rI;
	Vector_t Isc;
    };


    Intersection FindIntersection(Vector_t x0, Vector_t x1, size_t i) {
	//       INFOMSG(" LineInsTri called "<< endl);
	Vector_t lseg = x1-x0;
	Vector_t t0,t1,t2;
	t0 = geo3Dcoords_m[allbfaces_m[4*i+1]];
	t1 = geo3Dcoords_m[allbfaces_m[4*i+2]];
	t2 = geo3Dcoords_m[allbfaces_m[4*i+3]];
	Vector_t u = t1-t0;
	Vector_t v = t2-t0;
	Vector_t lt = t0-x0;
	Vector_t n;
	n[0] = u[1]*v[2]-v[1]*u[2];
	n[1] = u[2]*v[0]-v[2]*u[0];
	n[2] = u[0]*v[1]-v[0]*u[1];
	Intersection myInter;
	myInter.rI=0;
	myInter.Isc=x1;
	if ((n[0]*lseg[0]+n[1]*lseg[1]+n[2]*lseg[2])==0) {
	    //          INFOMSG("Triangle parallell to line segment, return x1 "<< endl);
	    return myInter;
	}
	else {

	    myInter.rI = (n[0]*lt[0]+n[1]*lt[1]+n[2]*lt[2])/(n[0]*lseg[0]+n[1]*lseg[1]+n[2]*lseg[2]);
	    Vector_t ItSec = x0+myInter.rI*lseg;
	    Vector_t w = ItSec-t0;
	    /* if ((rI<0)||(rI>1)) {
	    //	INFOMSG("Intesect is on the extended line, return x1 "<< endl);
	    return x1;
	    }*/
    
	    double tmp1,tmp2,tmp3,tmp4,tmp5;
	    tmp1 = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	    tmp2 = w[0]*v[0]+w[1]*v[1]+w[2]*v[2];
	    tmp3 = u[0]*w[0]+u[1]*w[1]+u[2]*w[2];
	    tmp4 = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
	    tmp5 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	    double sI = (tmp1*tmp2-tmp5*tmp3)/(tmp1*tmp1-tmp4*tmp5);
	    double tI = (tmp1*tmp3-tmp4*tmp2)/(tmp1*tmp1-tmp4*tmp5);
	    if ((sI>=0)&&(tI>=0)&&((sI+tI)<=1)) {
		myInter.Isc=ItSec;
		return myInter;
	    }
	    else {
		//		  INFOMSG("Intesect is on the extended plane, return x1 "<< endl);
		return myInter;
	    }
      
	}
    }


    /**
     * Return the larger one.
     */

    bool gT(double x,double y) {
	return (x>=y);//Used for initializing the particles 
    }
    /**
     * Leading zeros.
     * @param number stands for specified bits number 
     */ 

    string convert2Int(int number) {

	stringstream ss;//create a stringstream
	ss <<   setw(5) << setfill('0') <<  number; //add number to the stream
	return ss.str();//return a string with the contents of the stream
    }

    /**
     * A struct.
     * 
     * Used as a compare function object for STL max_element
     */

    struct myLessclassx {
	bool operator() (Vector_t x1, Vector_t x2) { return x1(0) < x2(0); }
    } myLessx;
    /**
     * A struct.
     *
     * Used as a compare function object for STL max_element
     */
    struct myLessclassy {
	bool operator() (Vector_t x1, Vector_t x2) { return x1(1) < x2(1); }
    } myLessy;
    /**
     * A struct.
     *
     * Used as a compare function object for STL max_element
     */
    struct myLessclassz {
	bool operator() (Vector_t x1, Vector_t x2) { return x1(2) < x2(2); }
    } myLessz;
    /**
     * Calculate the maximum of a triangle length
     * @param i stands for the triangle Id.
     */
    double Trianglelen(int i) {
	Vector_t x1,x2,x3;
	int id1=allbfaces_m[4*i+1],id2=allbfaces_m[4*i+2],id3=allbfaces_m[4*i+3];    
	x1=geo3Dcoords_m[id1];
	x2=geo3Dcoords_m[id2];
	x3=geo3Dcoords_m[id3];
	double max =  sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]))>sqrt((x3[0]-x2[0])*(x3[0]-x2[0])+(x3[1]-x2[1])*(x3[1]-x2[1])+(x3[2]-x2[2])*(x3[2]-x2[2]))?sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2])):sqrt((x3[0]-x2[0])*(x3[0]-x2[0])+(x3[1]-x2[1])*(x3[1]-x2[1])+(x3[2]-x2[2])*(x3[2]-x2[2]));
	max = max>sqrt((x3[0]-x1[0])*(x3[0]-x1[0])+(x3[1]-x1[1])*(x3[1]-x1[1])+(x3[2]-x1[2])*(x3[2]-x1[2]))?max:sqrt((x3[0]-x1[0])*(x3[0]-x1[0])+(x3[1]-x1[1])*(x3[1]-x1[1])+(x3[2]-x1[2])*(x3[2]-x1[2]));
	return ( max );
    } 
    /**
     * Return the distance between triangle barycentric point and a arbitrary point
     * @param i stands for triangle Id,@param x1 stands for particle coordinates
     */
    double Barycentlength(int i,Vector_t x1) {
	Vector_t x2;
	x2=Tribarycent_m[i];
	double length =  sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]));
	return ( length );
    }  
    /**
     * @param numbfaces_global_m stores the number of triangles used to represent the geometry
     */
    int     numbfaces_global_m;
    /**
     * @param allbfaces_m stores the IDs of triangle vertex, which can be used to reference the  coordinates of vertex 
     */
    int*    allbfaces_m;
    /**
     * @param bfaces_idx_m store the Id number of a triangle,not in use at present.
     */
    int*    bfaces_idx_m;
    /**
     * @param  nof_sym_m stores the count number of symmetry planes.
     */
    int     nof_sym_m    ; 
    /**
     * @param numpoints_global_m stores the number of geometry points to be used to represent the geometry
     * caution: not only the surface points are included
     */
    int     numpoints_global_m;
    /**
     * @param  Tribarycent_m store the coordinates of barycentric points of triangles, The Id number is the same with triangle Id.
     */
    Vector_t*    Tribarycent_m;
    /**
     * @param TriNormal_m store the inward normal vector of triangles, The Id number is the same with triangle Id.
     */
    vector<Vector_t>  TriNormal_m;
    /**
     * @param TriPartloss_m store the number of particles hitting the Id th triangle. The Id number is the same with triangle Id(not vertex ID).
     */
    int*    TriPartloss_m;
    /**
     * @param TriPartlossZ_m is a counter of lost particle number in each Z intervals.
     */
    int*    TriPartlossZ_m;//counter along Z for histogram
    /**
     * @param Triangle_lossmax_m store the maximum number dumped in each triangle.
     */
    int     Triangle_lossmax_m;
    /**
     * @param geo3Dcoords_m store the geometry coordinates in a STL Vector
     * The ID of geo3Dcoords_m is equal to points' ID, and can be referenced by triangle to get vertex coordinates 
     */
    vector<Vector_t> geo3Dcoords_m;
    /**
     * @param boundary_ids_m stores the Ids of triangles which form the boundary in a STL set
     */
    set<size_t> boundary_ids_m;
    /**
     * @param  b_cubcentre_m stores the Ids of triangle convex in cubic and cubic center coordinates
     */
    map<size_t,Vector_t> b_cubcentre_m;
 
    /**
     * @param len_m is a OPAL Vector_t type parameter stores the length of geometry in 3D Cartesian coordinates.
     */
    Vector_t      len_m;
    /**
     * @param hr_m is a OPAL Vector_t type parameter stores the length of cubic box.
     */
    Vector_t      hr_m;
    /**
     * @param nr_m is a OPAL Vector_t type parameter stores the number of intervals of geometry in 3D Cartesian coordinates.
     */
    Vektor<int,3> nr_m;
    /**
     * @param mincoords_m stores minimum of geometry coordinate.
     */
    Vector_t      mincoords_m;
    /**
     * @param maxcoords_m stores maximum of geometry coordinate.
     */
    Vector_t      maxcoords_m;
    /**
     * @param partsr_m stores coordinates of particles using STL vector.
     */
    vector<Vector_t> partsr_m;  // add particle positions 
    /**
     * @param partsp_m stores momenta of particles using STL vector.
     */
    vector<Vector_t> partsp_m;  // add particle momenta
    /**
     * @param lostedPart_m stores coordinates of losted particles
     */
    vector<Vector_t> lostedPart_m;  // store particles hitting the surface
};


int main(int argc,char *argv[]){
    Ippl ippl(argc, argv);
    Inform msg("DarkCurrentAlgoTest ");  
    Inform msg2all("DarkCurrentAlgoTest ",INFORM_ALL_NODES);  
    BoundaryGeometry myBg;
    myBg.initialize(string("ctf3.h5"));
    myBg.makeBoundaryIndexSet();
    myBg.writeGeomToVtk(string("vtk/testGeometry-00000.vtk"));

    //  myBg.makeTriangleIndexSet();
    myBg.testNp();
    INFOMSG("After myBg.testNp(); "<< endl);
    myBg.createParticles();
    INFOMSG("After myBg.createParticles(); "<< endl);
    myBg.writePartToVtk(string("vtk/testParticle-00000.vtk"));
    // myBg.writePartlossToASCII(string("vtk/testPartloss-00000.dat"));

    myBg.integPartcleMove(150,0.6);
    Ippl::Comm->barrier();
    return 0;
}
