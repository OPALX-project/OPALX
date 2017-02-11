#include "Distribution.h"

#include "H5Reader.h"
#include <fstream>

// ----------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// ----------------------------------------------------------------------------

Distribution::Distribution():
    x_m(0),
    y_m(0),
    z_m(0),
    px_m(0),
    py_m(0),
    pz_m(0),
    q_m(0),
    mass_m(0),
    nloc_m(0),
    ntot_m(0)
{ }


void Distribution::uniform(double lower, double upper, size_t nloc, int seed) {
    
    nloc_m = nloc;
    
    std::mt19937_64 mt(0/*seed*/ /*0*/);
    
    std::uniform_real_distribution<> dist(lower, upper);
    
//     // assume that seed == rank of node
//     mt.discard(6 * (nloc + 1) * seed);
    
    // assume that seed == rank of node    
    // inefficient but only way to make sure that parallel distribution is equal to sequential
    for (size_t i = 0; i < 6 * nloc_m * seed; ++i)
        dist(mt);
    
    x_m.resize(nloc);
    y_m.resize(nloc);
    z_m.resize(nloc);
    
    px_m.resize(nloc);
    py_m.resize(nloc);
    pz_m.resize(nloc);
    
    q_m.resize(nloc);
    mass_m.resize(nloc);
    
    for (size_t i = 0; i < nloc; ++i) {
        x_m[i] = dist(mt);
        y_m[i] = dist(mt);
        z_m[i] = dist(mt);
        
        px_m[i] = dist(mt);
        py_m[i] = dist(mt);
        pz_m[i] = dist(mt);
        
        q_m[i] = 1.0;
        mass_m[i] = 1.0;
    }
}

void Distribution::gaussian(double mean, double stddev, size_t nloc, int seed) {
    
    nloc_m = nloc;
    
    std::mt19937_64 mt(0/*seed*/ /*0*/);
    
    std::normal_distribution<double> dist(mean, stddev);
    
//     // assume that seed == rank of node
//     mt.discard(6 * (nloc + 1) * seed);

    // assume that seed == rank of node    
    // inefficient but only way to make sure that parallel distribution is equal to sequential
    for (size_t i = 0; i < 6 * nloc_m * seed; ++i)
        dist(mt);
    
    x_m.resize(nloc);
    y_m.resize(nloc);
    z_m.resize(nloc);
    
    px_m.resize(nloc);
    py_m.resize(nloc);
    pz_m.resize(nloc);
    
    q_m.resize(nloc);
    mass_m.resize(nloc);
    
    for (size_t i = 0; i < nloc; ++i) {
        x_m[i] = dist(mt);
        y_m[i] = dist(mt);
        z_m[i] = dist(mt);
        
        px_m[i] = dist(mt);
        py_m[i] = dist(mt);
        pz_m[i] = dist(mt);
        
        
        q_m[i] = 1.0;
        mass_m[i] = 1.0;
    }
}

void Distribution::twostream(const Vector_t& lower, const Vector_t& upper,
                             const Vektor<std::size_t, 3>& nx, const Vektor<std::size_t, 3>& nv,
                             const Vektor<double, 3>& vmax, double alpha)
{
    Vektor<double, 3> hx = (upper - lower) / Vector_t(nx);
    Vektor<double, 3> hv = 2.0 * vmax / Vector_t(nv);
    
//     std::cout << hx << std::endl << hv << std::endl;
    
    double thr = 1.0e-12;
    
    double factor = 1.0 / ( M_PI * 30.0 );
    
    nloc_m = 0;

    /* we parallelize only the longitudinal direction
     * since there we use the most grid points
     */
    int nMaxInitializer = std::min(nx[2], nv[2]);
    if ( Ippl::myNode() < nMaxInitializer ) {

        // number of processes really used for initialization of particles
        int nInitializer = ( Ippl::getNodes() < nMaxInitializer ) ? Ippl::getNodes() : nMaxInitializer;

        for (std::size_t i = 0; i < nx[0]; ++i) {
            for (std::size_t j = 0; j < nx[1]; ++j) {
                for ( std::size_t k = 0; k < nx[2]; ++k) {

                    Vektor<double, 3> pos = Vektor<double,3>(
                                                (0.5 + i) * hx[0] + lower[0],
                                                (0.5 + j) * hx[1] + lower[1],
                                                (0.5 + k) * hx[2] + lower[2]
                                            );
                    
                    for (std::size_t iv = 0; iv < nv[0]; ++iv) {
                        for (std::size_t jv = 0; jv < nv[1]; ++jv) {
                            for (std::size_t kv = 0; kv < nv[2]; ++kv) {

                                if ( Ippl::myNode() != int(kv) % nInitializer )
                                    continue;

                                Vektor<double, 3> vel = -vmax + hv *
                                                        Vektor<double, 3>(iv + 0.5,
                                                                          jv + 0.5,
                                                                          kv + 0.5);
                                
                                double v2 = vel[0] * vel[0] +
                                            vel[1] * vel[1] +
                                            vel[2] * vel[2];
                                
                                double f = factor * std::exp(-0.5 * v2) *
                                           (1.0 + alpha * std::cos(k * pos[2])) *
                                           (1.0 + /*0.*/5.0 * vel[2] * vel[2]);
                                
                                double m = hx[0] * hv[0] *
                                           hx[1] * hv[1] *
                                           hx[2] * hv[2] * f;
                                           
                                
                                if ( m > thr ) {
                                    ++nloc_m;
                                    x_m.push_back( pos[0] );
                                    y_m.push_back( pos[1] );
                                    z_m.push_back( pos[2] );
                                    
                                    px_m.push_back( vel[0] );
                                    py_m.push_back( vel[1] );
                                    pz_m.push_back( vel[2] );
                                    
//                                     total_charge += m;
                                    
                                    q_m.push_back( -m );
                                    mass_m.push_back( m );
                                }
                            }
                        }
                    }
                }
            }
        }
//         std::cout << "Total charge: " << total_charge << std::endl;
    }
}

void Distribution::uniformPerCell(const Array<Geometry>& geom,
                                  const Array<BoxArray>& ba,
                                  const Vektor<std::size_t, 3>& nr,
                                  std::size_t nParticles, int seed) {
    
    nloc_m = 0;
    
    std::mt19937_64 mt(0);
    std::uniform_real_distribution<> dist(0.0, 1.0);
    std::vector< Vektor<double, 3> > rn(nParticles);
        
    for (std::size_t i = 0; i < nParticles; ++i) {
	for (int d = 0; d < 3; ++d) {
	    rn[i](d) = dist(mt);
	}
    }
    
    double dx[3] = {0.0, 0.0, 0.0};     // cell size (a bit smaller such that particles not a cell boundary)
    double cidx[3] = {0, 0, 0}; // max. cell index for a level
        
    for (int d = 0; d < 3; ++d) 
	dx[d] = geom[0].CellSize(d);
    
        // map [0, 1] --> to cell dimension [0 + 0.25 * dx, 0.75 * dx]
    std::vector< Vektor<double, 3> > mapped2cell(nParticles);
    
    for (std::size_t pi = 0; pi < nParticles; ++pi) {
	for (int d = 0; d < 3; ++d)
	    mapped2cell[pi](d) = dx[d] * (0.5 * rn[pi](d) + 0.25);
    }
    
    for (int j = 0; j < ba[0].size(); ++j) {
	
	// make sure every processs gets a different chunk to initialize
	if ( Ippl::myNode() != j % Ippl::getNodes() )
	    continue;

	Box bx = ba[0].get(j);
	
	for (int k = bx.loVect()[0]; k <= bx.hiVect()[0]; ++k) {
	    for (int l = bx.loVect()[1]; l <= bx.hiVect()[1]; ++l) {
		for (int m = bx.loVect()[2]; m <= bx.hiVect()[2]; ++m) {
                        
		    // assign particle position
		    for (std::size_t pi = 0; pi < nParticles; ++pi) {
			// [index space] --> [physical domain]
			double kk = geom[0].ProbLength(0) / nr[0] * k + geom[0].ProbLo(0);
			double ll = geom[0].ProbLength(1) / nr[1] * l + geom[0].ProbLo(1);
			double mm = geom[0].ProbLength(2) / nr[2] * m + geom[0].ProbLo(2);
                        
			x_m.push_back( kk + mapped2cell[pi](0) );
			y_m.push_back( ll + mapped2cell[pi](1) );
			z_m.push_back( mm  + mapped2cell[pi](2) );
                        
			px_m.push_back( 1.0 );
			py_m.push_back( 1.0 );
			pz_m.push_back( 1.0 );
			q_m.push_back( 1.0 );
			mass_m.push_back( 1.0 );
                        
			++nloc_m;
		    }
		}
	    }
	}
    }
}


void Distribution::readH5(const std::string& filename, int step) {
    H5Reader h5(filename);
    
    h5.open(step);
    
    
    
    nloc_m = h5.getNumParticles();
    ntot_m = nloc_m;
    
    size_t numParticlesPerNode = nloc_m / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = nloc_m - 1;

    nloc_m = lastParticle - firstParticle + 1;
    
    x_m.resize(nloc_m);
    y_m.resize(nloc_m);
    z_m.resize(nloc_m);
    
    px_m.resize(nloc_m);
    py_m.resize(nloc_m);
    pz_m.resize(nloc_m);
    
    q_m.resize(nloc_m);
    mass_m.resize(nloc_m);
    
    h5.read(x_m, px_m,
            y_m, py_m,
            z_m, pz_m,
            q_m,
            mass_m,
            firstParticle,
            lastParticle);
    
    h5.close();
}


void Distribution::injectBeam(
#ifdef IPPL_AMR
    PartBunchAmr< ParticleAmrLayout<double, BL_SPACEDIM> >& bunch,
#else
    PartBunchBase& bunch,
#endif
    bool doDelete, std::array<double, 3> shift)
{
    
    // destroy all partcles
#ifndef IPPL_AMR
    if ( doDelete && bunch.getLocalNum() )
        bunch.destroyAll();
#endif
    
    // previous number of particles
    int prevnum = bunch.getLocalNum();

    // create memory space
    bunch.create(nloc_m);

    for (int i = nloc_m - 1; i >= 0; --i) {
#ifdef IPPL_AMR
      bunch.R[i + prevnum] = Vector_t(x_m[i] + shift[0], y_m[i] + shift[1], z_m[i] + shift[2]);
      bunch.P[i + prevnum] = Vector_t(px_m[i], py_m[i], pz_m[i]);
      bunch.qm[i + prevnum] = q_m[i];
      bunch.mass[i + prevnum] = mass_m[i];
#else
        bunch.setR(Vector_t(x_m[i] + shift[0], y_m[i] + shift[1], z_m[i] + shift[2]), i + prevnum);
        bunch.setP(Vector_t(px_m[i], py_m[i], pz_m[i]), i + prevnum);
        bunch.setQM(q_m[i], i + prevnum);
        bunch.setMass(mass_m[i], i + prevnum);
#endif
        
        x_m.pop_back();
        y_m.pop_back();
        z_m.pop_back();
        px_m.pop_back();
        py_m.pop_back();
        pz_m.pop_back();
        
        q_m.pop_back();
        mass_m.pop_back();
    }
}


void Distribution::setDistribution(
#ifdef IPPL_AMR
    PartBunchAmr< ParticleAmrLayout<double, BL_SPACEDIM> >& bunch,
#else
    PartBunchBase& bunch,
#endif
    const std::string& filename, int step)
{
    
    readH5(filename, step);
    
    for (unsigned int i = 0; i < bunch.getLocalNum(); ++i)
#ifdef IPPL_AMR
        bunch.R[i] = Vector_t(x_m[i], y_m[i], z_m[i]);
#else
        bunch.setR(Vector_t(x_m[i], y_m[i], z_m[i]), i);
#endif
}

void Distribution::print2file(std::string pathname) {
    
    for (int n = 0; n < Ippl::getNodes(); ++n) {
        
        if ( n == Ippl::myNode() ) {
            
            std::ofstream out;
            switch (n) {
                case 0:
                    out.open(pathname);
                    out << ntot_m << std::endl;
                    break;
                default:
                    out.open(pathname, std::ios::app);
                    break;
            }
            
            for (std::size_t i = 0; i < x_m.size(); ++i)
                out << x_m[i] << " " << px_m[i] << " "
                    << y_m[i] << " " << py_m[i] << " "
                    << z_m[i] << " " << pz_m[i] << std::endl;
            
            out.close();
        }
        
        Ippl::Comm->barrier();
    }
}
