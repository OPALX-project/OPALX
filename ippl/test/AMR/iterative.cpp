/*!
 * @file iterative.cpp
 * @author Matthias Frey
 * @date 19. Oct. 2016, LBNL
 * @details Solve \f$\Delta\phi = -1\f$ in 3D iteratively [nsteps] with
 * Dirichlet boundary conditions (zero)
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
 * @brief Solve \f$\Delta\phi = -1\f$ in 3D iteratively
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>


typedef std::vector<double> vector_t;
typedef std::vector<vector_t> matrix_t;
typedef std::vector<matrix_t> tensor_t;

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
// successive over relaxation, not working
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
    initSphereOnGrid(rho, a, R, n);
    
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