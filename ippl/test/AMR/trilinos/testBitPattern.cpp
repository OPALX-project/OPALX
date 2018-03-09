/*!
 * @file testBitPattern.cpp
 * @author Matthias Frey
 * 
 * @details Quadratic Lagrange interpolation for
 * coarse-to-fine interface in 3D. The interpolation
 * requires 9 coefficients, thus, 9 coarse cells
 * shouldn't be refined. With std::bitset we can
 * find the best pattern.
 * 
 *  ___________________________________
 * |      |      |      |      |      |
 * |  20  |  21  |  22  |  23  |  24  |
 * |______|______|______|______|______|
 * |      |      |      |      |      |
 * |  15  |  16  |  17  |  18  |  19  |
 * |______|______|______|______|______|
 * |      |      |      |      |      |
 * |  10  |  11  |  12  |  13  |  14  |
 * |______|______|______|______|______|
 * |      |      |      |      |      |
 * |  05  |  06  |  07  |  08  |  09  |
 * |______|______|______|______|______|
 * |      |      |      |      |      |
 * |  00  |  01  |  02  |  03  |  04  |
 * |______|______|______|______|______|
 * 
 * 
 * The interface cell is at 12.
 */

#include <bitset>
#include <iostream>
#include <iterator>
#include <vector>

#include <random>

void random(std::bitset<25>& test) {
//     std::mt19937_64 mt;
    std::random_device rd;
    std::uniform_int_distribution<> dist(0, 3);
    
    for (uint i = 0; i < test.size(); ++i) {
        test[i] = dist(rd);
    }
}
        

void fill(std::vector< std::bitset<25> >& pattern) {
    // cross pattern --> 6, 7, 8, 11, 12, 13, 16, 17, 18
    std::bitset<25> b0("0000001110011100111000000");
    pattern.push_back(b0);
    
    unsigned long int ib0 = b0.to_ulong();
    std::cout << b0 << " " << b0.to_ulong() << " " << ib0 << std::endl;
    
    // T pattern --> 1, 2, 3, 6, 7, 8, 11, 12, 13
    std::bitset<25> b1("0000000000011100111001110");
    pattern.push_back(b1);
    unsigned long int ib1 = b1.to_ulong();
    std::cout << b1 << " " << b1.to_ulong() << " " << ib1 << std::endl;
    
    // T on head pattern --> 11, 12, 13, 16, 17, 18, 21, 22, 23
    std::bitset<25> b2("0111001110011100000000000");
    pattern.push_back(b2);
    unsigned long int ib2 = b2.to_ulong();
    std::cout << b2 << " " << b2.to_ulong() << " " << ib2  << std::endl;
    
    // upper right corner pattern --> 2, 3, 4, 7, 8, 9, 17, 18, 19
    std::bitset<25> b3("0000011100000001110011100");
    pattern.push_back(b3);
    unsigned long int ib3 = b3.to_ulong();
    std::cout << b3 << " " << b3.to_ulong() << " " << ib3 << std::endl;
    
    // upper left corner pattern --> 0, 1, 2, 5, 6, 7, 10, 11, 12
    std::bitset<25> b4("0000000000001110011100111");
    pattern.push_back(b4);
    unsigned long int ib4 = b4.to_ulong();
    std::cout << b4 << " " << b4.to_ulong() << " " << ib4  << std::endl;
    
    // mirrored L pattern --> 10, 11, 12, 15, 16, 17, 20, 21, 22
    std::bitset<25> b5("0011100111001110000000000");
    pattern.push_back(b5);
    unsigned long int ib5 = b5.to_ulong();
    std::cout << b5 << " " << b5.to_ulong() << " " << ib5  << std::endl;
    
    // L pattern --> 12, 13, 14, 17, 18, 19, 22, 23, 24
    std::bitset<25> b6("1110011100111000000000000");
    pattern.push_back(b6);
    unsigned long int ib6 = b6.to_ulong();
    std::cout << b6 << " " << b6.to_ulong() << " " << ib6 << std::endl;
    
    // left hammer pattern --> 7, 8, 9, 12, 13, 14, 17, 18, 19
    std::bitset<25> b7("0000011100111001110000000");
    pattern.push_back(b7);
    unsigned long int ib7 = b7.to_ulong();
    std::cout << b7 << " " << b7.to_ulong() << " " << ib7  << std::endl;
    
    // right hammer pattern --> 5, 6, 7, 10, 11, 12, 15, 16, 17
    std::bitset<25> b8("0000000111001110011100000");
    pattern.push_back(b8);
    unsigned long int ib8 = b8.to_ulong();
    std::cout << b8 << " " << b8.to_ulong() << " " << ib8  << std::endl;
}
 
int main() 
{
    std::vector< std::bitset<25> > pattern;
    
    fill(pattern);
    
    std::bitset<25> test;
    random(test);
    
    std::cout << "Comparison:" << std::endl;
    std::cout << test << " " << test.to_ulong() << std::endl;
    
    
    bool found = false;
    std::vector< std::bitset<25> >::iterator it = std::begin(pattern);
    
    while ( !found && it != std::end(pattern) ) {
        if ( it->to_ulong() == (test & *it).to_ulong() )
            break;
        ++it;
    }
    
    switch ( it->to_ulong() ) {
        case 473536:
            std::cout << "Cross pattern." << std::endl;
            found = true;
            break;
        case 14798:
            std::cout << "T pattern." << std::endl;
            break;
        case 15153152:
            std::cout << "T on head pattern." << std::endl;
            break;
        case 918428:
            std::cout << "upper right corner pattern." << std::endl;
            break;
        case 7399:
            std::cout << "upper left corner pattern." << std::endl;
            break;
        case 7576576:
            std::cout << "mirrored L pattern." << std::endl;
            break;
        case 30306304:
            std::cout << "L pattern." << std::endl;
            break;
        case 947072:
            std::cout << "left hammer pattern." << std::endl;
            break;
        case 236768:
            std::cout << "right hammer pattern." << std::endl;
            break;
        default:
            std::cout << "Unknown pattern." << std::endl;
            break;
    }
}