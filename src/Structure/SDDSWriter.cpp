#include "SDDSWriter.h"

#include <queue>

SDDSWriter(const std::string& fname)
    : DataSink()
    , os_m(fname.c_str(), std::ios::out)
    , indent_m("        ")
{ }


void SDDSWriter::rewindLines(const std::string &fname, size_t numberOfLines) const {
    if (numberOfLines == 0) return;

    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;

    fs.open (fname.c_str(), std::fstream::in);

    if (!fs.is_open()) return;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fname.c_str(), std::fstream::out);

    if (!fs.is_open()) return;

    while (allLines.size() > numberOfLines) {
        fs << allLines.front() << "\n";
        allLines.pop();
    }
    fs.close();
}


void SDDSWriter::replaceVersionString(const std::string &fname) const {

    if (Ippl::myNode() == 0) {
        std::string versionFile;
        SDDS::SDDSParser parser(fname);
        parser.run();
        parser.getParameterValue("revision", versionFile);

        std::string line;
        std::queue<std::string> allLines;
        std::fstream fs;

        fs.open (fname.c_str(), std::fstream::in);

        if (!fs.is_open()) return;

        while (getline(fs, line)) {
            allLines.push(line);
        }
        fs.close();


        fs.open (fname.c_str(), std::fstream::out);

        if (!fs.is_open()) return;

        while (allLines.size() > 0) {
            line = allLines.front();

            if (line != versionFile) {
                fs << line << "\n";
            } else {
                fs << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << "\n";
            }

            allLines.pop();
        }

        fs.close();
    }
}


void SDDSWriter::addDescription(const std::string& text,
                                const std::string& content)
{
    os_m << "&description\n"
         << indent_m << "text=\"" << text << "\",\n"
         << indent_m << "contents=\"" << content << "\"\n"
         << "&end\n";
}


void SDDSWriter::addParameter(const std::string& name,
                              const std::string& type,
                              const std::string& desc)
{
    os_m << "&parameter\n"
         << indent_m << "name=" << name << ",\n"
         << indent_m << "type=" << type << ",\n"
         << indent_m << "description=\"" << desc << "\"\n"
         << "&end\n";
}


void SDDSWriter::addColumn(const std::string& name,
                           const std::string& type,
                           const std::string& unit,
                           const std::string& desc)
{
    static short column = 1;
    
    os_m << "&column\n"
         << indent_m << "name=" << name << ",\n"
         << indent_m << "type=" << type << ",\n"
         << indent_m << "units=" << unit << ",\n"
         << indent_m << "description=\"" << column++ << " " << desc << "\"\n"
         << "&end\n";
}


void SDDSWriter::addData(const std::string& mode,
                         const short& no_row_counts)
{
    os_m << "&data\n"
         << indent_m << "mode=" << mode << ",\n"
         << indent_m << "no_row_counts=" << no_row_counts << "\n"
         << "&end\n";
}


template<typename T>
void SDDSWriter::writeValue(const T& value) {
    os_m << value << std::setw(10) << "\t";
}
