// g++ Tests/generatesim.cpp -I .

#include "Util/GenerateSimulation.h"

#include <string>
#include <iostream>
#include <map>

int main(int argc, char** argv) {

    std::map<std::string, std::string> userVars;
    userVars.insert(std::pair<string, string>("NSLICE", "200"));

    GenerateSimulation *gs = new GenerateSimulation("Tests/gentest.tmpl", "Tests/gentest.data", userVars);
    gs->writeInputFile("Tests/gentest.in");
    delete gs;

    return 0;
}
