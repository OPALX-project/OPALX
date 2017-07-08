#include "Util/CmdArguments.h"

void CmdArguments::addArguments(int argc, char **argv) {

    for(int i=1; i<argc; i++) {
        std::string arg = argv[i];
        std::string name, value;
        this->split(name, value, arg);
        arguments_.insert(std::pair<std::string, std::string>(name, value));
    }
}

void CmdArguments::split(std::string &name,
                         std::string &value, std::string arg) {

    size_t pos = arg.find("=");
    //strip leading '--' and '='
    name = arg.substr(2, pos - 2);
    value = arg.substr(pos + 1);

    boost::trim(name);
    boost::trim(value);
}
