#ifndef mesh_hpp
#define mesh_hpp

#include <math.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <set>
#include <map>

#include "boost/bind.hpp"
#include "boost/function.hpp"

//for mesh writeout
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <cassert>
#include <stdio.h>

namespace Mesh {
  struct Output{
  public:
    //Control front/manifold output modes
    enum MODE {DEFAULT, FRONT};
    static MODE mode;
  };
  Output::MODE Output::mode;

  //Placeholder for a standard interface
  struct Interface {	
    virtual void Generate() {};
    virtual void Refresh() {};
    virtual void Print() =0;
    virtual void WriteOut( const char* ) =0;
  };

  typedef std::vector< double > point_t;
  typedef int ind_t;
  typedef std::vector< ind_t > coord_t;

  //Basic function for outputting vector to a stream
  template< typename T>
  
  
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec){
    unsigned i;      
    if( 1 ) {
      //if( Output::mode == Output::DEFAULT) {
      os << "[ " ;
      for(i=0;i<vec.size()-1;i++){ 
	os << vec[i];
	os << ", " ;	
      }
      os << vec[i] << " ]" ;
    }
    else{
      for(i=0;i<vec.size();i++){ 
	os << vec[i] ;
	os << std::endl ;	
      }
    }
    return os;
  }

  //Dis/Enable synchronous mesh updates
#define SYNCH

  class MeshBase : public Interface {
  public:
    class MeshPoint {
    public:
      //Maintain a global ID
      ind_t ID;
      //The three basic pieces of each meshpoint
      std::vector< double > DesignCoords;
      std::vector< double > ObjectiveCoords;
      std::vector< double > LambdaCoords;

      //Maintain neighbor references
      std::vector< MeshPoint* > NeighborP;

      //By index in the global array
      std::vector< int > Neighbors;
      
      MeshPoint(ind_t ID,
		std::vector< double > &Design, 
		std::vector< double > &Objective, 
		std::vector< double > &Lambdas ) :
	ID( ID ),
	DesignCoords( Design ), 
	ObjectiveCoords( Objective ), 
	LambdaCoords( Lambdas ) 
      { }
      
      //Simple print function that dumps all values
      void Print() {
	printf("ID: ");
	printf("%d ", this->ID );
	printf("Design: ( ");
	for(unsigned i=0;i<DesignCoords.size();i++)		
	  printf("%.15lf ", DesignCoords[i]);
	printf(") Obj: ( ");
	for(unsigned i=0;i<ObjectiveCoords.size();i++)	
	  printf("%.15lf ", ObjectiveCoords[i]);
	printf(") Lam: ( ");
	for(unsigned i=0;i<LambdaCoords.size();i++) 
	  printf("%lf ", LambdaCoords[i]);
	printf(") Links: ( ");
	for(unsigned i=0;i<Neighbors.size();i++)			
	  printf("%d ", Neighbors[i]);
	/*
	  printf(") Addresses: ( ");
	  for(i=0;i<NeighborP.size();i++)			printf("%p ", NeighborP[i]);
	  printf(") this: ( ");
	  printf("%p ", this);
	*/
	printf(")\n");
      }
      
      //For outputting a meshpoint to the front file
      friend std::ostream& operator<<(std::ostream& os, const MeshPoint& mp){
	//if( Output::mode == Output::MANIFOLD ) return os << mp.DesignCoords;
	//else return os << mp.ObjectiveCoords;
	os << "{" << std::endl;
	os << "\"ID\": " << mp.ID << " , " << std::endl;
	unsigned i;
	for(i=0;i<mp.ObjectiveCoords.size();i++){
	  os << "\"o" << i << "\": " << mp.ObjectiveCoords[i] << " , " << std::endl;
	}
	for(i=0;i<mp.DesignCoords.size()-1;i++){
	  os << "\"d" << i << "\": " << mp.DesignCoords[i] << " , " << std::endl;
	}
	os << "\"d" << i << "\": " << mp.DesignCoords[i] << std::endl;
	os << "}" << std::endl;
	return os;
      }

    };
    
    bool isValid;
    std::vector< MeshPoint > Corners;
    std::vector< MeshPoint > Points;
    
#ifdef SYNCH
  private:
    std::vector< MeshPoint > ShadowPoints;
#endif
    
  public:
    
    unsigned DesignDim;
    unsigned ObjectiveDim;
    unsigned MeshDim;
    
    //Verbose Constructor
    MeshBase(	std::vector< std::vector< double > > DesignSpace,
		std::vector< std::vector< double > > ObjectiveSpace,
		std::vector< std::vector< double > > Lambdas ) : 
      DesignDim( DesignSpace[0].size() ), 
      ObjectiveDim( ObjectiveSpace[0].size() ),
      MeshDim( ObjectiveDim - 1 )
    {
      
      //verify that we have the right number of points
      assert( DesignSpace.size() == ObjectiveSpace.size() );
      assert( ObjectiveSpace.size() == Lambdas.size() );
      
      unsigned i;
      
      for(i=0; i<DesignSpace.size(); i++) {
	
	//verify the dimensionality of each point
	assert( DesignDim == DesignSpace[i].size() );
	assert( ObjectiveDim == ObjectiveSpace[i].size() );
	assert( ObjectiveDim == Lambdas[i].size() );
	
	MeshPoint p( i, DesignSpace[i], ObjectiveSpace[i], Lambdas[i] );
	this->Corners.push_back( p );
      }
      
      this->isValid = false;
    }
    
    virtual std::vector< point_t* > GetNeighborLocOf(ind_t i){
      std::vector< point_t* > r( this->Points[i].NeighborP.size(), NULL );
      
      for(unsigned j=0; j<r.size(); j++) r[j] = &( this->Points[i].NeighborP[j]->ObjectiveCoords);
      return r;
    }
    
#ifdef SYNCH
    virtual void Generate() {
      std::vector< MeshPoint >::iterator i;
      for(i=Points.begin(); i!=Points.end(); i++){
	this->ShadowPoints.push_back( *i );
      }
      
      this->isValid = true;
    }
    virtual void UpdatePoint( ind_t i, point_t& d, point_t& o, point_t& l ) { 
      this->ShadowPoints[i].DesignCoords.assign(      d.begin(), d.end() );
      this->ShadowPoints[i].ObjectiveCoords.assign(   o.begin(), o.end() );
      this->ShadowPoints[i].LambdaCoords.assign(      l.begin(), l.end() );
    }
    virtual void Refresh() {
      assert( this->ShadowPoints.size() == this->Points.size() );
	for(unsigned i=0; i<this->ShadowPoints.size(); i++){
	  this->Points[i].DesignCoords.assign(    this->ShadowPoints[i].DesignCoords.begin(), 
						  this->ShadowPoints[i].DesignCoords.end() );
	  this->Points[i].ObjectiveCoords.assign( this->ShadowPoints[i].ObjectiveCoords.begin(), 
						  this->ShadowPoints[i].ObjectiveCoords.end() );
	  this->Points[i].LambdaCoords.assign(    this->ShadowPoints[i].LambdaCoords.begin(), 
						  this->ShadowPoints[i].LambdaCoords.end() );
	}
    }
#else
    virtual void Generate() {
      this->isValid = true;
    }
    virtual void UpdatePoint( ind_t i, point_t& d, point_t& o, point_t& l ) { 
      this->Points[i].DesignCoords.assign(      d.begin(), d.end() );
      this->Points[i].ObjectiveCoords.assign(   o.begin(), o.end() );
      this->Points[i].LambdaCoords.assign(      l.begin(), l.end() );
    }
#endif
    virtual bool Valid() { return isValid; }

    virtual void Print() {
      this->Refresh();

      printf("Corners: \n");
      for(unsigned i=0;i<Corners.size();i++) Corners[i].Print();
      printf("Points: \n");
      for(unsigned i=0;i<Points.size();i++) Points[i].Print();
    }
    
    void WriteOut(const char* filename ) {
      
      std::ofstream os( filename );
      os << "{" << std::endl;

      /*
	std::vector<std::string> cols;
	unsigned i;
	for(i=0;i<numcols;i++) {
	std::stringstream s;
	s << "\"" << i << "\"" ;
	cols.push_back( std::string( s.str() ) );
	}
	os << "\"columns\" : " << cols << "," << std::endl;
      */
      
      os << "\"name\" : \"opt-pilot v0.1\"," << std::endl;
      
      os << "\"solutions\" : ";
      os << this->Points;
      
      os << "}" << std::endl;
      os.close();
    }
  };
  
  class Simplex : public MeshBase {
  private:
  public:
    static int eta(int dim, int n) {
      return eta_num(dim, n) / eta_den(dim);
    }
    static int eta_num(int dim, int n) {
      if (dim == 1) return n;
      else return (n+dim -1) * eta_num(dim-1, n);
    }
    static int eta_den(int dim) {
      if (dim == 1) return 1;
      else return dim * eta_den(dim - 1);
    }
    static int coord_to_ind(std::vector< int > &coords, int n){
      return coord_to_ind_aux(coords, n, 0);
    }
    static int coord_to_ind_aux(std::vector< int > &coords, int n, int ind) {
      if ( ind == (int) coords.size()-1 ) return coords[ind];
      else {
	int this_coord = coords[ind];
	int dim = coords.size() - ind;
	return ( eta(dim, n) - eta(dim, n - this_coord) ) + coord_to_ind_aux(coords, n-this_coord, ind+1);
      }
    }

    typedef boost::function<int ( int )> func;
    static int bisection_search(func & f, int a, int b){
      assert ( (f)(a) * (f)(b) <= 0 );
      if( (f)(a) == 0 ) return a;
      if( (f)(b) == 0 ) return b;
      while ( fabs(b - a) > 1.0 ) {
	int sum = a+b;
	int c = ((sum)/2) + ( (sum % 2) ? 1 : 0 ) ;
	if( (f)(c) == 0 ) return c;
	if( (f)(a) * (f)(c) <= 0 ) b = c;
	else a = c;
      }
      return a;
    }
    static int eta_helper(int dim, int ind, int x){ return eta(dim, x) - ind; }
    static std::vector< int > ind_to_coord(int ind, int dim, int n){
      //flip the coordinate system
      ind = eta(dim, n) - 1 - ind;

      std::vector< int > coords(dim);
      int i;
      for(i=0; i<dim; i++){
	func f ( boost::bind( &eta_helper, dim-i, ind, _1 ) );
	int c = bisection_search(f, 0, n);
	coords[i] = c;
	ind -= eta(dim - i, c);
      }

      std::vector< int >::iterator it;
      for(it = coords.begin(); it != coords.end(); it++){
	*it = n - std::accumulate(coords.begin(), it, 0) - *it -1;
      }
      return coords;
    }

    typedef boost::function<double ( int )> baryTransFunc;
    static double barycentric_trans_helper(double n, int x){ return x*n; }
    static double barycentric_trans_helper_2(double n, double x, double y){ 
      return y + n*x; 
    }
    static bool isBoundaryPoint( std::vector< int > &coords, int n) {
      return ( std::accumulate( coords.begin(), coords.end(), 0 ) == (n - 1) );
    }
    static std::vector< std::vector< int > > getNeighborsAux( std::vector< int > coords, int n){
      std::vector< std::vector< int > > result;

      if( ! isBoundaryPoint( coords, n ) ) {
	unsigned i;
	for(i=0;i<coords.size();i++) {
	  if( coords[i] == 0 ) continue;

	  coords[i] += 1;
	  result.push_back( coords );					
	  coords[i] -= 2;
	  result.push_back( coords );
	  coords[i] += 1;
	}
      }
      else{
	coords.pop_back();

	if( !coords.empty() ) {
	  result = getNeighborsAux( coords, n );
	  unsigned i;
	  for(i=0;i<result.size();i++){
	    int sum = std::accumulate( result[i].begin(), result[i].end(), 0);
	    result[i].push_back( n - 1 - sum );
	  }
	}
      }
      return result;
    }
    static std::vector< int > getNeighbors( std::vector< int > coords, int n ){
      std::vector< int > neighbors;
      std::vector< std::vector< int > > NCoords(getNeighborsAux(coords, n) );
      unsigned i;
      for(i=0;i<NCoords.size();i++){ neighbors.push_back( coord_to_ind(NCoords[i], n) ); }
			
      return neighbors;
    }

    static std::vector< std::vector< int > > getNeighborhoodAux( std::vector< int > coords, int n){
      std::vector< std::vector< int > > result;

      if( ! isBoundaryPoint( coords, n ) ) {
	//Not a boundary point
	unsigned i;
	for(i=0;i<coords.size();i++) {
	  if( coords[i] == 0 ) {
	    //Only has neighborhood in positive direction
	    coords[i] += 1;
	    result.push_back( coords );	
	    coords[i] -= 1;
	  }
	  else {
	    //Has neighborhood in both directions
	    coords[i] += 1;
	    result.push_back( coords );					
	    coords[i] -= 2;
	    result.push_back( coords );
	    coords[i] += 1;
	  }
	}
      }
      else{
	//Is a boundary point
                
	unsigned c;
	for(c=0;c<coords.size();c++) {
	  if( coords[c] == 0 ) continue;

	  //Only has neighborhood in negative direction
	  coords[c] -= 1;
	  result.push_back( coords );	
	  coords[c] += 1;
	}
                
	coords.pop_back();

	//Find the diagonal "surface" neighbors
	if( !coords.empty() ) {
	  std::vector< std::vector< int > > surface;
	  surface = getNeighborhoodAux( coords, n );
	  unsigned i;
	  for(i=0;i<surface.size();i++){
	    int sum = std::accumulate( surface[i].begin(), surface[i].end(), 0);
	    surface[i].push_back( n - 1 - sum );
	  }
	  result.insert( result.end(), surface.begin(), surface.end() );
	}
      }
      return result;
    }
    static std::vector< int > getNeighborhood( std::vector< int > coords, int n ){
      std::vector< int > neighbors;
      std::vector< std::vector< int > > NCoords(getNeighborhoodAux(coords, n) );
      unsigned i;
      for(i=0;i<NCoords.size();i++){ neighbors.push_back( coord_to_ind(NCoords[i], n) ); }
      return neighbors;
    }

    //public:
    static std::set< int > GetNeighborIDs( int ID, int Dimension, int PPS ){
      //Get a vector of neighbors
      std::vector< int > vec (getNeighborhood( ind_to_coord( ID, Dimension, PPS ), PPS ) );
      //Copy to a set
      std::set< int > retval(vec.begin(), vec.end()); 
      return retval;
    }
    std::set< int > GetNeighborIDs( int ID ){
      //Get a vector of neighbors
      std::vector< int > vec (getNeighborhood( ind_to_coord( ID, this->MeshDim, this->PointsPerSide ), this->PointsPerSide ) );
      //Copy to a set
      std::set< int > retval(vec.begin(), vec.end()); 
      return retval;
    }

    int PointsPerSide;

    Simplex(	std::vector< std::vector< double > > DesignSpace,
		std::vector< std::vector< double > > ObjectiveSpace,
		std::vector< std::vector< double > > Lambdas, 
		int NumberOfPoints ) : MeshBase( DesignSpace, ObjectiveSpace, Lambdas),
				       PointsPerSide( NumberOfPoints ) {}

    virtual void Generate(){
      //Generate points
      baryTransFunc baryTrans( boost::bind( &barycentric_trans_helper, 1.0/(this->PointsPerSide-1), _1 ) );
      int i;
      for(i=0; i<eta(this->MeshDim, this->PointsPerSide); i++){

	//Get the point n-index
	std::vector< int > v( Mesh::Simplex::ind_to_coord( i, this->MeshDim, this->PointsPerSide) ) ;

	//Get the neighbor indicies
	std::vector< int > neighbors ( Mesh::Simplex::getNeighbors( v, this->PointsPerSide ) );

	//Generate the barycentric coordinates of this point
	v.push_back( PointsPerSide - 1 - std::accumulate(v.begin(), v.end(), 0) );
	std::vector< double > barycentric( v.size() );
	std::transform( v.begin(), v.end(), barycentric.begin(), baryTrans );

	//Generate point location in various spaces from
	//barycentric coordinates
	std::vector< double > design( this->DesignDim, 0.0 );
	std::vector< double > objective( this->ObjectiveDim, 0.0 );
	std::vector< double > lambda( this->ObjectiveDim, 0.0 );
	unsigned j;
	for(j=0; j<barycentric.size(); j++){
	  //printf("%lf\t", barycentric[j]);
	  boost::function<double ( double, double )> comb ( 
							   boost::bind( &barycentric_trans_helper_2, 
									barycentric[j], _1, _2) );
	  std::transform( this->Corners[j].DesignCoords.begin(), this->Corners[j].DesignCoords.end(), design.begin(), design.begin(), comb);
	  std::transform( this->Corners[j].ObjectiveCoords.begin(), this->Corners[j].ObjectiveCoords.end(), objective.begin(), objective.begin(), comb);
	  std::transform( this->Corners[j].LambdaCoords.begin(), this->Corners[j].LambdaCoords.end(), lambda.begin(), lambda.begin(), comb);
	}
	//printf("\n");

	MeshPoint m( i, design, objective, lambda );
	m.Neighbors = neighbors;

	this->Points.push_back( m );
      }

      MeshBase::Generate();
    }

    virtual void Refresh() {
      MeshBase::Refresh();
    }
  };

  class SimplexNodeSet {
  public:
    friend class Simplex;
    //Return the number of subsets
    static ind_t NumberOfSubsets( ind_t dim, ind_t SubsetsPerSide ) { return Simplex::eta( dim, SubsetsPerSide ); }
    //Map global coordinates to subset index
    static ind_t WhichSubset( coord_t coords, ind_t PointsPerSide, ind_t SubsetsPerSide ){
      ind_t PointsPerSubset = PointsPerSide / SubsetsPerSide;
      ind_t b = PointsPerSide  - PointsPerSubset * SubsetsPerSide;
      //Rescale coordinates to Node dimensions
      unsigned i;
      for(i=0; i<coords.size(); i++) {
	if( coords[i] > b*(PointsPerSubset+1) ) coords[i] = b + (coords[i] - b*(PointsPerSubset+1))/PointsPerSubset;
	else coords[i] = coords[i]/(PointsPerSubset+1);
      }
      //Get node index from node coordinates
      return Simplex::coord_to_ind( coords, SubsetsPerSide );
    }
  };

  class SimplexSubset : public Simplex {
  public:
    ind_t ID;
    ind_t SubsetsPerSide;
    coord_t origin;

    //TODO: rework this for better encapsulation
    std::vector< MeshPoint* > GhostNodes;
    std::map< ind_t, MeshPoint* > GlobalToGhost;
    std::map< ind_t, std::vector< MeshPoint* > > NeighborSubsetsToLocals;

    SimplexSubset(      std::vector< point_t > DesignSpace,
			std::vector< point_t > ObjectiveSpace,
			std::vector< point_t > Lambdas, 
			int NumberOfPoints, ind_t ID, ind_t SubsetsPerSide) : 
      Simplex( DesignSpace, ObjectiveSpace, Lambdas, NumberOfPoints ),
      ID( ID ), SubsetsPerSide( SubsetsPerSide ) {
      this->origin = Simplex::ind_to_coord( ID, this->MeshDim, this->SubsetsPerSide );
    }

    virtual void Generate(){
      std::map< ind_t, ind_t > GlobalToLocal;
      std::map< ind_t, ind_t > LocalToGlobal;

      //Generate points
      baryTransFunc baryTrans( boost::bind( &barycentric_trans_helper, 1.0/(this->PointsPerSide-1), _1 ) );
      int i;
      for(i=0; i<eta(this->MeshDim, this->PointsPerSide); i++){

	//Get the point n-index
	std::vector< int > v( Mesh::Simplex::ind_to_coord( i, this->MeshDim, this->PointsPerSide) ) ;

	//Only generate the point if it is in this subset
	if( this->ID == SimplexNodeSet::WhichSubset( v, this->PointsPerSide, this->SubsetsPerSide ) ){

	  //Get the neighbor indicies
	  std::vector< int > neighbors ( Mesh::Simplex::getNeighbors( v, this->PointsPerSide ) );

	  //Generate the barycentric coordinates of this point
	  v.push_back( PointsPerSide - 1 - std::accumulate(v.begin(), v.end(), 0) );
	  std::vector< double > barycentric( v.size() );
	  std::transform( v.begin(), v.end(), barycentric.begin(), baryTrans );

	  //Generate point location in various spaces from
	  //barycentric coordinates
	  std::vector< double > design( this->DesignDim, 0.0 );
	  std::vector< double > objective( this->ObjectiveDim, 0.0 );
	  std::vector< double > lambda( this->ObjectiveDim, 0.0 );
	  unsigned j;
	  for(j=0; j<barycentric.size(); j++){
	    boost::function<double ( double, double )> comb ( 
							     boost::bind( &barycentric_trans_helper_2, 
									  barycentric[j], _1, _2) );
	    std::transform( this->Corners[j].DesignCoords.begin(), this->Corners[j].DesignCoords.end(), design.begin(), design.begin(), comb);
	    std::transform( this->Corners[j].ObjectiveCoords.begin(), this->Corners[j].ObjectiveCoords.end(), objective.begin(), objective.begin(), comb);
	    std::transform( this->Corners[j].LambdaCoords.begin(), this->Corners[j].LambdaCoords.end(), lambda.begin(), lambda.begin(), comb);
	  }

	  MeshPoint m( i, design, objective, lambda );
	  m.Neighbors = neighbors;

	  //Temporary mapping for later optimizing member indexing
	  GlobalToLocal[ i ] = this->Points.size();
	  LocalToGlobal[ this->Points.size() ] = i;

	  this->Points.push_back( m );
	}
      }

      //Once the local points are generated, we need to 
      //re-index the neighbor references to the local scheme
      i = 0;
      std::vector< MeshPoint >::iterator p;
      for(p=this->Points.begin(); p!=this->Points.end(); p++, i++){

	//For each of the points in the neighborhood (adjacency graph)
	std::vector< ind_t > adjacent( Simplex::getNeighborhood( Simplex::ind_to_coord( p->ID, this->MeshDim, this->PointsPerSide ), this->PointsPerSide ) );
	std::vector< ind_t >::iterator a;
	for(a=adjacent.begin(); a!=adjacent.end(); a++){
	  //Is it non-local?
	  ind_t sub = SimplexNodeSet::WhichSubset( Mesh::Simplex::ind_to_coord( *a, this->MeshDim, this->PointsPerSide ), this->PointsPerSide, this->SubsetsPerSide );
	  if( sub != this->ID ){
	    //Is this current point a neighbor of the non-local adjacent point
	    std::vector< ind_t > a_n ( Simplex::getNeighbors( Simplex::ind_to_coord( *a, this->MeshDim, this->PointsPerSide ), this->PointsPerSide ) );
	    if( std::find( a_n.begin(), a_n.end(), p->ID ) != a_n.end() ){
	      //(in which case we need to send this location)
	      //So register with the ghost manager
	      if( NeighborSubsetsToLocals.find( sub ) !=  NeighborSubsetsToLocals.end() ){
		NeighborSubsetsToLocals[ sub ].push_back( &(*p) );
	      }
	      else{
		std::vector< MeshPoint* > v( 1, &(*p) );
		NeighborSubsetsToLocals[ sub ] = v;
	      }
	    }
	  }

	}

	unsigned j;
	for(j=0; j<p->Neighbors.size(); j++){
	  //Is the neighbor a local point?
	  std::map< ind_t, ind_t >::iterator it = GlobalToLocal.find( p->Neighbors[j] );
	  if( it != GlobalToLocal.end() ){
	    p->NeighborP.push_back( &(this->Points[it->second]) );
	  }
	  //If not, create and manage a ghost node
	  else{
	    MeshPoint* ghost;
	    point_t empty;

	    try{ 
	      ghost = new MeshPoint( p->Neighbors[j], empty, empty, empty );
	    }
	    catch ( std::bad_alloc& e)
	      {
		ghost = NULL;
		printf("Error allocating memory for ghost node %d\n", p->Neighbors[j] );
	      }
                         
	    //Register the ghost with the container (for destruction)
	    GhostNodes.push_back( ghost );
	    //Register the ghost with the point
	    p->NeighborP.push_back( ghost );
	    //Register the ghost with the ghost manager
	    GlobalToGhost[ p->Neighbors[j] ] = ghost;
	    //Register this local point with the ghost manager
	    ind_t ghostid = SimplexNodeSet::WhichSubset( Mesh::Simplex::ind_to_coord( p->Neighbors[j], this->MeshDim, this->PointsPerSide ), this->PointsPerSide, this->SubsetsPerSide );
	    //Have we seen this subset yet?
	    if( NeighborSubsetsToLocals.find( ghostid ) !=  NeighborSubsetsToLocals.end() ){
	      NeighborSubsetsToLocals[ ghostid ].push_back( &(*p) );
	    }
	    else{
	      std::vector< MeshPoint* > v( 1, &(*p) );
	      NeighborSubsetsToLocals[ ghostid ] = v;
	    }
	  }
	}
      }

      MeshBase::Generate();
    }

    virtual bool Valid() {
      if( MeshBase::Valid() ){
	bool valid = true;

	std::vector< MeshPoint* >::iterator i=GhostNodes.begin();
	while( valid && i!=GhostNodes.end() ){
	  valid &= ! (*i)->ObjectiveCoords.empty();

	  i++;
	}

	return valid;
      }
      else return false;
    }

    virtual void Refresh() {
      MeshBase::Refresh();
    }

    virtual std::vector< MeshPoint* > GetNeighborsOf(ind_t i){
      return std::vector< MeshPoint* > ( this->Points[i].NeighborP.begin(), this->Points[i].NeighborP.end() );
    }

    virtual void Print(){
      Simplex::Print();
      printf("Incoming Ghosts: \n");
      std::map< ind_t, MeshPoint* >::iterator i;
      for(i=GlobalToGhost.begin(); i!=GlobalToGhost.end(); i++){
	i->second->Print();
      }
      printf("Outgoing Ghosts: \n");
      std::map< ind_t, std::vector< MeshPoint* > >::iterator it;
      for(it=NeighborSubsetsToLocals.begin(); it!=NeighborSubsetsToLocals.end(); it++){
	printf("To %d: ", it->first);
	std::vector< MeshPoint* >::iterator ind;
	for(ind=it->second.begin(); ind!=it->second.end(); ind++) printf("%d ", (*ind)->ID);
	printf("\n");
      }
    }
    void WriteOut() {
      std::ostringstream fout;
      fout << "front." << this->ID;
      MeshBase::WriteOut( fout.str().c_str() );
      std::ostringstream mout;
      mout << "manifold." << this->ID;
      MeshBase::WriteOut( mout.str().c_str() );
    }
  };
}

#endif
