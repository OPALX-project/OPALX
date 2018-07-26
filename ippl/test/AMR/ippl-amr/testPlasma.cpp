#include "Ippl.h"

#include "PlasmaPIC.h"


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    char* in_file_name = argv[1];
    amrex::Initialize(argc, argv /*false*/);
    
    Inform msg(argv[0]);
    
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    
    IpplTimings::startTimer(mainTimer);
    
    
    try {
        
        PlasmaPIC pic;
        
        pic.execute(msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    IpplTimings::stopTimer(mainTimer);
    
    IpplTimings::print();
    
    amrex::Finalize(true);
    
    return 0;
}
