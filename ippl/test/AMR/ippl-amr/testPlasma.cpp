#include "Ippl.h"
#include <iostream>

#include "PlasmaPIC.h"


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    amrex::Initialize(argc,argv, false);
    
    Inform msg(argv[0]);
    
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    
    IpplTimings::startTimer(mainTimer);
    
    try {
        
        PlasmaPIC pic;
        
        pic.execute();
        
    } catch(const std::exception& ex) {
        std::cout << ex.what() << std::endl;
    }
    
    IpplTimings::stopTimer(mainTimer);
    
    IpplTimings::print();
    
    amrex::Finalize(true);
    
    return 0;
}
