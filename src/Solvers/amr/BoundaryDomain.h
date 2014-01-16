#ifndef BOUNDARY_DOMAIN
#define BOUNDARY_DOMAIN

#include <vector>
#include <cmath>
#include <iostream> 
#include <map>
#include <string>
#include <assert.h> 
//#include "Structure/BoundaryGeometry.h"

class BoundaryDomain{

public:

   BoundaryDomain(std::vector<int>& nr, const std::vector<double>& hr) 
   {
      SemiMajor = 0.02; //bgeom->getA(); 
      SemiMinor = 0.02; //bgeom->getB(); 
     
      // setNr(nr);
      // setHr(hr);

      FindIntersection(nr,hr);
   };

    ~BoundaryDomain();

    /// BoundaryPointList maps from an (x,z) resp. (y,z) pair to double values (=intersections with boundary)
    typedef std::multimap< std::pair<int, int>, double >  BoundaryPointList;

    /// all intersection points with grid lines in X direction
    BoundaryPointList IntersectLoX, IntersectHiX;

    /// all intersection points with grid lines in Y direction
    BoundaryPointList IntersectLoY, IntersectHiY;    

    /// all intersection points with grid lines in Z direction
    BoundaryPointList IntersectLoZ, IntersectHiZ;    

    BoundaryPointList& GetIntersectLoX()
    {
       return IntersectLoX;
    }
    BoundaryPointList& GetIntersectHiX()
    {
       return IntersectHiX;
    }
    BoundaryPointList& GetIntersectLoY()
    {
       return IntersectLoX;
    }
    BoundaryPointList& GetIntersectHiY()
    {
       return IntersectHiX;
    }

    void FindIntersection(const std::vector<int>& nr_c,const std::vector<double>& hr_c){

  	   std::multimap< std::pair<int, int>, double >::iterator it;
  	   std::pair< std::multimap< std::pair<int, int>, double>::iterator, std::multimap< std::pair<int, int>, double>::iterator > ret;
 	   // bool hit;

	    IntersectLoX.clear();
	    IntersectHiX.clear();
	    IntersectLoY.clear();
	    IntersectHiY.clear();

	    double smajsq = SemiMajor*SemiMajor;
	    double sminsq = SemiMinor*SemiMinor;

 	   //calculate intersection with the ellipse
	    for(int z = 0; z < nr_c[2]; z++) 
            {
        	    for(int x = 0; x < nr_c[0]; x++) 
                    {

                	std::pair<int, int> pos(x, z);

			double posx = x*hr_c[0] -(nr_c[0]-1)*hr_c[0]/2.0;

                        if (posx <= -SemiMajor || posx >= SemiMajor)
                        {
                	    double ylo = 0.0;
                    	    double yhi = 0.0;
       	                    IntersectHiY.insert(std::pair< std::pair<int, int>, double >(pos, yhi));
                	    IntersectLoY.insert(std::pair< std::pair<int, int>, double >(pos, ylo));		

                        } else {

			    double yhi = std::abs(sqrt(sminsq - sminsq * posx * posx / smajsq));
   			    double ylo = -yhi;
       	                    IntersectHiY.insert(std::pair< std::pair<int, int>, double >(pos, yhi));
                	    IntersectLoY.insert(std::pair< std::pair<int, int>, double >(pos, ylo));		
			}
            	    }

            	   for(int y = 0; y < nr_c[1]; y++) 
                   {
                	std::pair<int, int> pos(y, z);
			double posy = y*hr_c[1] -(nr_c[1]-1)*hr_c[1]/2.0;

                        if (posy <= -SemiMinor || posy >= SemiMinor)
                        {
                	    double xlo = 0.0;
                    	    double xhi = 0.0;
       	                    IntersectHiX.insert(std::pair< std::pair<int, int>, double >(pos, xhi));
                	    IntersectLoX.insert(std::pair< std::pair<int, int>, double >(pos, xlo));		

                        } else {

                	    double xhi = std::abs(sqrt(smajsq - smajsq * posy * posy / sminsq));
                   	    double xlo = -xhi;
       	                    IntersectHiX.insert(std::pair< std::pair<int, int>, double >(pos, xhi));
                	    IntersectLoX.insert(std::pair< std::pair<int, int>, double >(pos, xlo));		
			}
             	   }
      	    }
    };

#if 0
   inline bool isInside(int x, int y, int z) {
    	double cx = (x - (nr[0]-1)/2)*hr[0];
   	double cy = (y - (nr[1]-1)/2)*hr[1];
	double val1 = 0; double val2 = 0;

	std::multimap < std::pair<int, int>, double >::iterator itr;
    	bool ret = false;

    	//check if x is inside with y,z coords
    	std::pair<int, int> coordyz(y, z);
    	itr = IntersectLoX.find(coordyz);
    	int count = IntersectLoX.count(coordyz);
    	if(count == 1)
         ret = (cx == itr->second);
    	else {
       	 val1 = itr->second;
       	 val2 = (++itr)->second;
       	 ret = (cx <= std::max(val1, val2)) && (cx >= std::min(val1, val2));
    	}

    	//check if y is inside with x,z coords
    	std::pair<int, int> coordxz(x, z);
   	itr = IntersectLoY.find(coordxz);
    	count = IntersectLoY.count(coordxz);
    	if(count == 1)
       	 ret = ret && (cy == itr->second);
   	else {
       	 val1 = itr->second;
       	 val2 = (++itr)->second;
       	 ret = ret && (cy <= std::max(val1, val2)) && (cy >= std::min(val1, val2));
	}

        return ret;
     };

   void LinearInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor){

        double dx_lo, dx_hi, dy_lo, dy_hi;

    	double cx = x*hr[0] - (nr[0]-1)*hr[0]/2.0;
    	double cy = y*hr[1] - (nr[1]-1)*hr[1]/2.0;

    	std::multimap< std::pair<int, int>, double >::iterator it;
    	std::pair< std::multimap< std::pair<int, int>, double>::iterator, std::multimap< std::pair<int, int>, double>::iterator > ret;

    	std::pair<int, int> coordyz(y, z);
   	ret = IntersectLoX.equal_range(coordyz);
    	for(it = ret.first; it != ret.second; ++it) {
        	if(fabs(it->second - cx) < hr[0]) {
           	 	dx_lo = it->second;
	  	 	dx_hi = -dx_lo;
            		break;
        	}
    	}

    	std::pair<int, int> coordxz(x, z);
    	ret = IntersectLoY.equal_range(coordxz);
    	for(it = ret.first; it != ret.second; ++it) {
        	if(fabs(it->second - cy) < hr[1]) {
            	dy_lo = it->second;
            	break;
        	}
    	}


        scaleFactor = 1.0;

    double dw=hr[0];
    double de=hr[0];
    double dn=hr[1];
    double ds=hr[1];
    C = 0.0;

    // we are a right boundary point
    if(!isInside(x+1,y,z)) {
        C += 1/((dx-cx)*de);
        E = 0.0;
    } else {
        C += 1/(de*de);
        E = -1/(de*de);
    }

    // we are a left boundary point
    if(!isInside(x-1,y,z)) {
        C += 1/((abs(dx)-abs(cx))*dw);
        W = 0.0;
    } else {
        C += 1/(dw*dw);
        W = -1/(dw*dw);
    }

    // we are a upper boundary point
    if(!isInside(x,y+1,z)) {
        C += 1/((dy-cy)*dn);
        N = 0.0;
    } else {
        C += 1/(dn*dn);
        N = -1/(dn*dn);
    }

    // we are a lower boundary point
    if(!isInside(x,y-1,z)) {
        C += 1/((abs(dy)-abs(cy))*ds);
        S = 0.0;
    } else {
        C += 1/(ds*ds);
        S = -1/(ds*ds);
    }

    F = -1/(hr[2]*hr[2]); 
    B = -1/(hr[2]*hr[2]);

    //XXX: In stand-alone only Dirichlet for validation purposes
    if(z == 0 || z == nr[2]-1) {

        // Dirichlet
        C += 2/hr[2]*1/hr[2];

        //C += 1/hr[2]*1/hr[2];

        // case where we are on the Neumann BC in Z-direction
        // where we distinguish two cases  
        if(z == 0)
            F = 0.0;
        else 
            B = 0.0;

        //for test no neumann 
        //C += 2/((hr[2]*nr[2]/2.0) * hr[2]);
        /*
           double d = hr[2]*(nr[2])/2;
           C += 2/(d * hr[2]);


        //neumann stuff
        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;

        scaleFactor *= 0.5;
        */

    } else 
        C += 2*1/hr[2]*1/hr[2];
};
#endif



#if 0
    void getBoundaryStencil(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor){

    // determine which interpolation method we use for points near the boundary
    if(interpolationMethod == "constant")
        ConstantInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor); 
    else if(interpolationMethod == "linear")
        LinearInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor); 
    else if(interpolationMethod == "quadratic")
        QuadraticInterpolation(x,y,z,W,E,S,N,F,B,C,scaleFactor); 
    // stencil center value has to be positive!
    assert(C > 0);
    };

    void getBoundaryStencil(int idx, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor){

    int x = 0, y = 0, z = 0;
    getCoord(idx,x,y,z);
    getBoundaryStencil(x,y,z,W,E,S,N,F,B,C,scaleFactor);
	};
    void getNeighbours(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B){
    if(x > 0)
        W = getIdx(x-1,y,z);
    else
        W = -1;
    if(x < nr[0]-1)
        E = getIdx(x+1,y,z);
    else
        E = -1;

    if(y < nr[1]-1)
        N = getIdx(x,y+1,z);
    else
        N = -1;
    if(y > 0)
        S = getIdx(x,y-1,z);
    else
        S = -1;

    if(z > 0)
        F = getIdx(x,y,z-1);
    else
        F = -1;
    if(z < nr[2]-1)
        B = getIdx(x,y,z+1);
    else
        B = -1;

};
    void getNeighbours(int idx, double& W, double& E, double& S, double& N, double& F, double& B) {

    int x = 0, y = 0, z = 0;
    getCoord(idx,x,y,z);
    getNeighbours(x,y,z,W,E,S,N,F,B);
};

    std::string getType() { return "Elliptic"; }
    int getNumXY(int z) { return nxy_m; }
    void setSemiMinor(double sm) {SemiMinor = sm;}
    void setSemiMajor(double sm) {SemiMajor = sm;}

    std::vector<int> getNr()    { return nr; }
    std::vector<double> getHr() { return hr; }
    void setNr(const std::vector<int>& nr_m)    { nr = nr_m; }
    void setHr(const std::vector<double>& hr_m) { hr = hr_m; }

    void getCoord(int idx, int& x, int& y, int& z) {
        int ixy = idx % nxy_m;
        int xy = CoordMap[ixy];
        x = xy % nr[0];
        y = (xy-x)/nr[0];
        z = (idx-ixy)/nxy_m + 1; 
    }

    inline int getIdx(int x, int y, int z) {
        if(isInside(x,y,z) && x >= 0 && y >= 0 && z >= 0)
            return IdxMap[toCoordIdx(x,y)]+(z-1)*nxy_m; 
        else 
            return -1;
    }
#endif

private:

    /// mapping (x,y,z) -> idx
    std::map<int, int> IdxMap;

    /// mapping idx -> (x,y,z)
    std::map<int, int> CoordMap;

    double SemiMajor;
    double SemiMinor;
    int nxy_m;
    std::string interpolationMethod ="linear";


    inline int toCoordIdx(int x, int y) { return y*nr[0] + x; }

#if 0
    /// different interpolation methods for boundary points 
    void ConstantInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor){

    scaleFactor = 1.0;

    W = -1/(hr[0]*hr[0]);
    E = -1/(hr[0]*hr[0]);
    N = -1/(hr[1]*hr[1]);
    S = -1/(hr[1]*hr[1]);
    F = -1/(hr[2]*hr[2]);
    B = -1/(hr[2]*hr[2]);
    C = 2/(hr[0]*hr[0]) + 2/(hr[1]*hr[1]) + 2/(hr[2]*hr[2]);

    // we are a right boundary point
    if(!isInside(x+1,y,z)) {
        E = 0.0;
    } 

    // we are a left boundary point
    if(!isInside(x-1,y,z)) {
        W = 0.0;
    } 

    // we are a upper boundary point
    if(!isInside(x,y+1,z)) {
        N = 0.0;
    } 

    // we are a lower boundary point
    if(!isInside(x,y-1,z)) {
        S = 0.0;
    } 

    // case where we are on Neumann BC in Z-direction
    // where we distinguish two cases  
    if(z == 0) 
        F = 0.0;
    if(z == nr[2]-1) 
        B = 0.0;

    //XXX: In stand-alone only Dirichlet for validation purposes
    // Neumann stuff
    //W /= 2.0;
    //E /= 2.0;
    //N /= 2.0;
    //S /= 2.0;
    //C /= 2.0;
};


    void QuadraticInterpolation(int x, int y, int z, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double& scaleFactor){

    double cx = x*hr[0] - nr[0]*hr[0]/2.0;
    double cy = y*hr[1] - nr[1]*hr[1]/2.0;

    //since every vector for elliptic domains has ALWAYS size 2 we 
    //can catch all cases manually

    std::multimap<int, double>::iterator it = IntersectLoX.find(y);
    if(cx < 0)
        it++; 

    it = IntersectLoY.find(x);
    if(cy < 0)
        it++;

    //double dw=hr[0];
    //double de=hr[0];
    //double dn=hr[1];
    //double ds=hr[1];
    W = 1.0;
    E = 1.0;
    N = 1.0;
    S = 1.0;
    F = 1.0;
    B = 1.0;
    C = 0.0;

    //TODO: finish cleanup
};
#endif

protected:
    /// number of mesh points in each direction
    std::vector<int> nr;
    /// mesh-spacings in each direction
    std::vector<double> hr;

};

#endif //#ifdef BOUNDARY_DOMAIN
