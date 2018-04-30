#include <vector>
#include <string>
#include <sstream>

#include "boost/uuid/sha1.hpp"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

/**
 *  \brief Generates a hash name.
 *
 *  Concatenates and hashes a vector of strings.
 */
class UUIDNameGenerator {

public:

    static std::string generate(std::vector<std::string> arguments,
                                size_t world_pid = 0) {

        std::string hash_string = "";

        foreach( std::string arg, arguments ) {
            hash_string += arg;
        }

        // append PID
        std::ostringstream pidStr;
        pidStr << world_pid;
        hash_string += "_" + pidStr.str();

        boost::uuids::detail::sha1 hash;
        hash.process_bytes(hash_string.c_str(), hash_string.size());
        unsigned int digest[5];
        hash.get_digest(digest);
        char hash_final[20];
        for(int i=0; i < 5; i++) {
            const char* tmp     = reinterpret_cast<char*>(digest);
            hash_final[i*4]     = tmp[i*4+3];
            hash_final[i*4 + 1] = tmp[i*4+2];
            hash_final[i*4 + 2] = tmp[i*4+1];
            hash_final[i*4 + 3] = tmp[i*4];
        }

        std::ostringstream hash_out;
        hash_out << std::hex;

        for(int i=0; i < 20; i++) {
            hash_out << ((hash_final[i] & 0x000000F0) >> 4);
            hash_out <<  (hash_final[i] & 0x0000000F);
        }

        return hash_out.str();
    }

private:

    UUIDNameGenerator() {}
    ~UUIDNameGenerator() {}

};
