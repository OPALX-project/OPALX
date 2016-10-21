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
    
    if (argc != 5) {
        std::cerr << "./iterative [nx] [ny] [nz] [nsteps]" << std::endl;
        return -1;
    }
    
    int nr[3] = {
        std::atoi(argv[1]),
        std::atoi(argv[2]),
        std::atoi(argv[3]),
    };
    
    std::cout << "nx = " << nr[0] << std::endl
              << "ny = " << nr[1] << std::endl
              << "nz = " << nr[2] << std::endl;
    
    // charge density initialized with -1
    double rho = -1.0;
    
    // unknowns (initialized to zero by default)
    /*
     * the dimension in each direction has to be
     * + 2 since we have Dirichlet boundary
     * conditions and we just consider interior points
     */
    tensor_t phi(nr[0] + 2,
                 matrix_t(nr[1] + 2,
                          vector_t(nr[2] + 2)
                         )
                );
    
    
    tensor_t newphi(nr[0] + 2,
                 matrix_t(nr[1] + 2,
                          vector_t(nr[2] + 2)
                         )
                );
    
    
    // mesh sizes
    double h[3] = { 1.0 / double(nr[0]),
                    1.0 / double(nr[1]),
                    1.0 / double(nr[2])
    };
    
    double inx = 1.0 / (h[0] * h[0]);
    double iny = 1.0 / (h[1] * h[1]);
    double inz = 1.0 / (h[2] * h[2]);
    
    double fac = -1.0 / (2.0 * ( inx + iny + inz ));
    
    // solve (boundary is zero)
    for (int t = 0; t < std::atoi(argv[4]); ++t) {
        for (int i = 1; i < nr[0] + 1; ++i)
            for (int j = 1; j < nr[1] + 1; ++j)
                for (int k = 1; k < nr[2] + 1; ++k) {
                    newphi[i][j][k] = fac * ( rho
                                                    - inx * (phi[i+1][j][k] + phi[i-1][j][k])
                                                    - iny * (phi[i][j+1][k] + phi[i][j-1][k])
                                                    - inz * (phi[i][j][k+1] + phi[i][j][k-1])
                                                    );
                }
        std::swap(phi, newphi);
    }
    
    
    
    // just write interior points
    std::ofstream pout("phi.dat");
    for (int i = 1; i < nr[0] + 1; ++i)
        for (int j = 1; j < nr[1] + 1; ++j)
            for (int k = 1; k < nr[2] + 1; ++k)
                pout << i << " " << j << " " << k << " "
                     << phi[i][j][k] << std::endl;
    
    pout.close();
    
    /*
     * compute electric field in center
     * (longitudinal direction)
     * using central difference
     */
    int half = 0.5 * nr[2];
    
    // electric field in x
    std::ofstream exout("ex.dat");
    for (int i = 1; i < nr[0] + 1; ++i)
        for (int j = 1; j < nr[1] + 1; ++j)
            exout << i - 1 << " " << j - 1 << " "
                  << 0.5 * (phi[i + 1][j][half] - phi[i - 1][j][half]) / h[0]
                  << std::endl;
    
    exout.close();
    
    // electric field in y
    std::ofstream eyout("ey.dat");
    for (int i = 1; i < nr[0] + 1; ++i)
        for (int j = 1; j < nr[1] + 1; ++j)
            eyout << i - 1 << " " << j - 1 << " "
                  << 0.5 * (phi[i][j + 1][half] - phi[i][j - 1][half]) / h[1]
                  << std::endl;
    
    eyout.close();
    
    
    // electric field in z (in x half)
    half = 0.5 * nr[0];
    std::ofstream ezout("ez.dat");
    for (int j = 1; j < nr[1] + 1; ++j)
        for (int k = 1; k < nr[2] + 1; ++k)
            ezout << j - 1 << " " << k - 1 << " "
                  << 0.5 * (phi[half][j][k + 1] - phi[half][j][k - 1]) / h[2]
                  << std::endl;
    
    ezout.close();
    
    return 0;
}