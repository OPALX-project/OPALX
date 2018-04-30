#include <vector>
#include <string>
#include <sstream>

#include "boost/uuid/uuid.hpp"
#include "boost/uuid/uuid_io.hpp"
#include "boost/uuid/name_generator.hpp"
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

        for (const std::string &arg: arguments ) {
            hash_string += arg;
        }

        // append PID
        std::ostringstream pidStr;
        pidStr << world_pid;
        hash_string += "_" + pidStr.str();

        boost::uuids::name_generator_sha1 hash(boost::uuids::ns::dns());
        return boost::uuids::to_string(hash(hash_string));
    }

private:

    UUIDNameGenerator() {}
    ~UUIDNameGenerator() {}

};