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

char** CmdArguments::getArguments() const {
    const unsigned int size = arguments_.size();
    char** args = new char*[2 * size];

    unsigned int i = 0;
    auto it = arguments_.cbegin();
    const auto end = arguments_.cend();
    for (; it != end; ++ it) {
        const std::string &key = it->first;
        char* argname = new char[key.length() + 1];
        strcpy(argname, key.c_str());
        args[i ++] = argname;

        const std::string &value = it->second;
        char* argvalue = new char[value.length() + 1];
        strcpy(argvalue, value.c_str());
        args[i ++] = argvalue;
    }

    return args;
}