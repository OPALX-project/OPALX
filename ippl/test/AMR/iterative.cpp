/*!
 * @file iterative.cpp
 * @author Matthias Frey
 * @date 19. Oct. 2016, LBNL
 * @details Solve \f$\Delta\phi = -\rho\f$ in 3D iteratively [nsteps] with
 *          Dirichlet boundary conditions (zero) \n
 *          Domain [-0.05, 0.05]^3\n
 *          There are several right-hand sides available.
 * 
 * Compiling:
 *      g++ -std=c++11 iterative.cpp -o iterative
 * 
 * Running:
 *      ./iterative [#gridpoints] [nsteps]
 * 
 * Visualization:
 *      The program writes 4 files (phi.dat, ex.dat, ey.dat, ez.dat) that
 *      can be visualized using the Python script vis_iter.py
 *      (run: python vis_iter.py)
 * @brief Solve \f$\Delta\phi = -\rho\f$ in 3D iteratively
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>


typedef std::vector<double> vector_t;
typedef std::vector<vector_t> matrix_t;
typedef std::vector<matrix_t> tensor_t;

/*! Solves Poisson equation \f$ \Delta\phi = -\rho \f$ in 3D
 *  using the Jacobi method. It uses Dirichlet boundary
 *  conditions (zero), where the boundary is the edge.
 * @param phi is the potential (unknowns)
 * @param h is the mesh spacing
 * @param n is the number of discretization points
 * @param rho is the right-hand side
 * @param nSteps is the number of iteration steps
 */
void jacobi(tensor_t& phi, const double& h, const int& n, const tensor_t& rho, const int& nSteps) {
    
    tensor_t newphi(n + 2,
                 matrix_t(n + 2,
                          vector_t(n + 2)
                         )
                );
    
    
    double fac = - h * h / 6.0;
    double inv = 1.0 / 6.0;
    
    // solve (boundary is zero)
    /*
     * boundary is the edge and not the nodal point
     * we can get a zero edge when we imply
     * that
     * \phi[i-1] = -\phi[i] for i = 1 and i = n
     * 
     * (also for j and k direction)
     */
    for (int t = 0; t < nSteps; ++t) {
        for (int i = 1; i < n + 1; ++i)
            for (int j = 1; j < n + 1; ++j)
                for (int k = 1; k < n + 1; ++k) {
                    newphi[i][j][k] = fac * rho[i-1][j-1][k-1] + inv * (phi[i+1][j][k] + phi[i-1][j][k] +
                                                                  phi[i][j+1][k] + phi[i][j-1][k] +
                                                                  phi[i][j][k+1] + phi[i][j][k-1]
                                                                 );
                }
        std::swap(phi, newphi);
        
        
        // update boundary
        for (int i = 0; i < n + 2; ++i)
            for (int j = 0; j < n + 2; ++j) {
                phi[i][j][0] = -phi[i][j][1];
                phi[i][j][n+1] = -phi[i][j][n];
            }
        
        for (int i = 0; i < n + 2; ++i)
            for (int k = 0; k < n + 2; ++k) {
                phi[i][0][k] = -phi[i][1][k];
                phi[i][n+1][k] = -phi[i][n][k];
            }
        
        for (int j = 0; j < n + 2; ++j)
            for (int k = 0; k < n + 2; ++k) {
                phi[0][j][k] = -phi[1][j][k];
                phi[n+1][j][k] = -phi[n][j][k];
            }
    }
}

// 6.12.2016, http://www.physics.buffalo.edu/phy410-505/2011/topic3/app1/index.html
/*!
 * successive over relaxation, not working
 * @param phi is the potential (unknowns)
 * @param h is the mesh spacing
 * @param n is the number of discretization points
 * @param rho is the right-hand side
 * @param nSteps is the number of iteration steps
 */
void sor(tensor_t& phi, const double& h, const int& n, const double& rho, const int& nSteps) {
    
    double w = 2.0 / ( 1.0 + M_PI / double(n) );
    
    tensor_t newphi(n + 2,
                 matrix_t(n + 2,
                          vector_t(n + 2)
                         )
                );
    
    
    double fac = - h * h / 6.0;
    double inv = 1.0 / 6.0;
    
    // solve (boundary is zero)
    /*
     * boundary is the edge and not the nodal point
     * we can get a zero edge when we imply
     * that
     * \phi[i-1] = -\phi[i] for i = 1 and i = n
     * 
     * (also for j and k direction)
     */
    for (int t = 0; t < nSteps; ++t) {
        for (int i = 1; i < n + 1; ++i)
            for (int j = 1; j < n + 1; ++j)
                for (int k = 1; k < n + 1; ++k) {
                    newphi[i][j][k] = (1.0 - w) * phi[i][j][k]
                                      + w * fac * rho + w * inv * (phi[i+1][j][k] + phi[i-1][j][k] +
                                                                   phi[i][j+1][k] + phi[i][j-1][k] +
                                                                   phi[i][j][k+1] + phi[i][j][k-1]
                                                                   );
                }
        std::swap(phi, newphi);
        
        
        // update boundary
        for (int i = 0; i < n + 2; ++i)
            for (int j = 0; j < n + 2; ++j) {
                phi[i][j][0] = -phi[i][j][1];
                phi[i][j][n+1] = -phi[i][j][n];
            }
        
        for (int i = 0; i < n + 2; ++i)
            for (int k = 0; k < n + 2; ++k) {
                phi[i][0][k] = -phi[i][1][k];
                phi[i][n+1][k] = -phi[i][n][k];
            }
        
        for (int j = 0; j < n + 2; ++j)
            for (int k = 0; k < n + 2; ++k) {
                phi[0][j][k] = -phi[1][j][k];
                phi[n+1][j][k] = -phi[n][j][k];
            }
    }
}

/*! Write the electric field and the potential on the grid point
 *  to files.
 * @param phi is the potential on the grid
 * @param h is the mesh spacing
 * @param n is the number of discretization points
 */
void write(tensor_t& phi, const double& h, const int& n) {
    // just write interior points
    std::ofstream pout("phi.dat");
    for (int i = 1; i < n + 1; ++i)
        for (int j = 1; j < n + 1; ++j)
            for (int k = 1; k < n + 1; ++k)
                pout << i << " " << j << " " << k << " "
                     << phi[i][j][k] << std::endl;
    
    pout.close();
    
    /*
     * compute electric field in center
     * (longitudinal direction)
     * using central difference
     */
    int half = 0.5 * n;
    
    // electric field in x
    std::ofstream exout("ex.dat");
    for (int i = 1; i < n + 1; ++i)
        for (int j = 1; j < n + 1; ++j)
            exout << i - 1 << " " << j - 1 << " "
                  << 0.5 * (phi[i + 1][j][half] - phi[i - 1][j][half]) / h
                  << std::endl;
    
    exout.close();
    
    // electric field in y
    std::ofstream eyout("ey.dat");
    for (int i = 1; i < n + 1; ++i)
        for (int j = 1; j < n + 1; ++j)
            eyout << i - 1 << " " << j - 1 << " "
                  << 0.5 * (phi[i][j + 1][half] - phi[i][j - 1][half]) / h
                  << std::endl;
    
    eyout.close();
    
    
    // electric field in z (in x half)
    half = 0.5 * n;
    std::ofstream ezout("ez.dat");
    for (int j = 1; j < n + 1; ++j)
        for (int k = 1; k < n + 1; ++k)
            ezout << j - 1 << " " << k - 1 << " "
                  << 0.5 * (phi[half][j][k + 1] - phi[half][j][k - 1]) / h
                  << std::endl;
    
    ezout.close();
}

/*!
 * Initialize the charge distribution on the grid with
 * charge 1.0 for each grid point.
 * @param rho is the right-hand side to be filled
 * @param a is the domain size (cubic) [m]
 * @param R is the radius of the sphere [m]
 * @param n is the number of discretization points
 */
void initSphereOnGrid(tensor_t& rho, double a, double R, int n) {
    
    double eps = 8.854187817e-12;
    double q = 2.78163e-15; // 1.112650056e-14;
    int num = 0;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
//                 rho[i][j][k] = -1.0;
                double x = 2.0 * a / double(n - 1) * i - a;
                double y = 2.0 * a / double(n - 1) * j - a;
                double z = 2.0 * a / double(n - 1) * k - a;
                
                if ( x * x + y * y + z * z <= R * R ) {
                    rho[i][j][k] = -1.0; //- q / (4.0 * M_PI * eps);
                    ++num;
                }
            }
        }
    }
    
    
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                
                double x = 2.0 * a / double(n - 1) * i - a;
                double y = 2.0 * a / double(n - 1) * j - a;
                double z = 2.0 * a / double(n - 1) * k - a;
                
                if ( x * x + y * y + z * z <= R * R ) {
//                     rho[i][j][k] /= double(num);
                    sum += rho[i][j][k];// * (4.0 * M_PI * eps);
                }
            }
        }
    }
    
    std::cout << "Total Charge: " << sum << " [C]" << std::endl;
    std::cout << "#Grid non-zero points: " << num << std::endl;
}

/*!
 * Initialize the charge distribution on the grid with
 * charge 1.0 on the whole domain.
 * @param rho is the right-hand side to be filled
 * @param n is the number of discretization points
 */
void initMinusOneEverywhere(tensor_t& rho, int n) {
    
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                rho[i][j][k] = -1.0;
                sum += 1.0;
            }
        }
    }
    std::cout << "Total Charge: " << sum << " [C]" << std::endl;
}

/*!
 * Initialize the charge distribution on the grid within
 * a sphere of radius R. The total charge is
 * \f$ q = 2.78163\cdot 10^{-15} \f$.
 * @param rho is the right-hand side to be filled
 * @param a is the domain size (cubic) [m]
 * @param R is the radius of the sphere [m]
 * @param n is the number of discretization points
 */
void initSphere(tensor_t& rho, double a, double R, int n) {
    
    // vacuum permittivity [C / (V m)]
    long double eps = 8.854187817e-12;
    
    // total charge
    long double q = 2.78163e-15; // 1.112650056e-14;
    int num = 0;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                double x = 2.0 * a / double(n - 1) * i - a;
                double y = 2.0 * a / double(n - 1) * j - a;
                double z = 2.0 * a / double(n - 1) * k - a;
                
                if ( x * x + y * y + z * z <= R * R ) {
                    rho[i][j][k] = - q / (4.0 * M_PI * eps);
                    ++num;
                }
            }
        }
    }
    
    
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                
                double x = 2.0 * a / double(n - 1) * i - a;
                double y = 2.0 * a / double(n - 1) * j - a;
                double z = 2.0 * a / double(n - 1) * k - a;
                
                if ( x * x + y * y + z * z <= R * R ) {
                    rho[i][j][k] /= double(num);
                    sum += rho[i][j][k] * (4.0 * M_PI * eps);
                }
            }
        }
    }
    
    std::cout << "Total Charge: " << sum << " [C]" << std::endl;
}

/*!
 * Interpolate the particle \f$ (x, y, z) \f$ to the nearest
 * grid point \f$ (i, j, k) \f$. Called by initSphereNGP().
 * @param i is the grid point in x-direction (computed)
 * @param j is the grid point in y-direction (computed)
 * @param k is the grid point in z-direction (computed)
 * @param x is the horizontal position of the particle [m]
 * @param y is the vertical position of the particle [m]
 * @param z is the longitudinal position of the particle [m]
 * @param h is the mesh spacing
 * @param a is side length of the cubic domain [m]
 * @param n is the number of discretization points
 */
void nearestGridPoint(int& i, int& j, int& k,
                      double x, double y, double z,
                      double h, double a, int n)
{
    // map from [-a, a] --> [0, n - 1]
    double itmp = (n - 1) / ( 2.0 * a ) * x + 0.5 * (n - 1);
    double jtmp = (n - 1) / ( 2.0 * a ) * y + 0.5 * (n - 1);
    double ktmp = (n - 1) / ( 2.0 * a ) * z + 0.5 * (n - 1);
    
    i = itmp;
    j = jtmp;
    k = ktmp;
    
    if ( std::abs( h - itmp +  i  ) >= 0.5 * h )
        ++i;
    
    if ( std::abs( h - jtmp +  j  ) >= 0.5 * h )
        ++j;
    
    if ( std::abs( h - ktmp +  k  ) >= 0.5 * h )
        ++k;
}

/*!
 * Initialize a particle distribution of \f$ 10^6 \f$ particles
 * within a sphere of radius R and assign the charges to the grid
 * using nearest grid point interpolation. The total charge
 * is defined to be \f$ q = 2.78163\cdot 10^{-15} \f$.
 * @param rho is the right-hand side to be filled
 * @param a is the side length of the cubic domain [m]
 * @param R is the radius of the sphere [m]
 * @param n is the number of discretization points
 */
void initSphereNGP(tensor_t& rho, double a, double R, int n) {
    
    double h = 2 * a / double(n);
    int nParticles = 1e6;
    
    std::mt19937_64 eng;
    std::uniform_real_distribution<> ph(-1.0, 1.0);
    std::uniform_real_distribution<> th(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<> u(0.0, 1.0);
    
    long double eps = 8.854187817e-12;
    long double qi = 4.0 * M_PI * eps * R * R / double(nParticles);
    
    for (int pi = 0; pi < nParticles; ++pi) {
        // 17. Dec. 2016,
        // http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
        // http://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
        double phi = std::acos( ph(eng) );
        double theta = th(eng);
        double radius = R * std::cbrt( u(eng) );
        
        double x = radius * std::cos( theta ) * std::sin( phi );
        double y = radius * std::sin( theta ) * std::sin( phi );
        double z = radius * std::cos( phi );
        
        int i = 0, j = 0, k = 0;
        nearestGridPoint(i, j, k, x, y, z, h, a, n);
        rho[i][j][k] -= qi / (4.0 * M_PI * eps);
    }
}
    
int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        std::cerr << "./iterative [nx] [nsteps]" << std::endl;
        return -1;
    }
    
    int n = std::atoi(argv[1]); // grid points
    double R = 0.005; // m
    double a = 0.05;  // m
    
    double h = 2.0 * a / double(n); // mesh size
    int nSteps = std::atoi(argv[2]);
    
    // charge density
    tensor_t rho(n,
                 matrix_t(n,
                          vector_t(n)
                 )
             );
    
//     initMinusOneEverywhere(rho, n);
//     initSphereOnGrid(rho, a, R, n);
    initSphere(rho, a, R, n);
//     initSphereNGP(rho, a, R, n);
    
    // unknowns (initialized to zero by default)
    /*
     * the dimension in each direction has to be
     * + 2 since we have Dirichlet boundary
     * conditions and we just consider interior points
     */
    tensor_t phi(n + 2,
                 matrix_t(n + 2,
                          vector_t(n + 2)
                         )
                );
    
    jacobi(phi, h, n, rho, nSteps);
//     sor(phi, h, n, rho, nSteps);
    
    write(phi, h, n);
    
    return 0;
}