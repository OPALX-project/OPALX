/* Author: Matthias Frey
 * Date: 19. Oct. 2016, LBNL
 * 
 * Solve \Lap(\phi) = -1 in 3D iteratively [nsteps] with
 * Dirichlet boundary conditions (zero)
 * 
 * Compiling:
 *      g++ -std=c++11 iterative.cpp -o iterative
 * 
 * Running:
 *      ./iterative [#gridpoints x] [#gridpoints y] [#gridpoints z] [nsteps]
 * 
 * Visualization:
 *      The program writes 4 files (phi.dat, ex.dat, ey.dat, ez.dat) that
 *      can be visualized using the Python script vis_iter.py
 *      (run: python vis_iter.py)
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>


typedef std::vector<double> vector_t;
typedef std::vector<vector_t> matrix_t;
typedef std::vector<matrix_t> tensor_t;


int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        std::cerr << "./iterative [nx] [nsteps]" << std::endl;
        return -1;
    }
    
    int n = std::atoi(argv[1]); // grid points
    double h = 1.0 / double(n); // mesh size
    
    // charge density initialized with -1
    double rho = -1.0;
    
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
    for (int t = 0; t < std::atoi(argv[2]); ++t) {
        for (int i = 1; i < n + 1; ++i)
            for (int j = 1; j < n + 1; ++j)
                for (int k = 1; k < n + 1; ++k) {
                    newphi[i][j][k] = fac * rho + inv * (phi[i+1][j][k] + phi[i-1][j][k] +
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
    
    return 0;
}