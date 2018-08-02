#include <fstream>
#include <iostream>
#include <sstream>

#include "boost/algorithm/string.hpp"
#include "boost/tuple/tuple.hpp"

#include "Expression/SumErrSq.h"
#include "Expression/FromFile.h"

#include "CommentAnnotatedInputFileParser.h"
#include "Util/OptPilotException.h"

typedef std::pair<std::string, DVar_t>  namedDVar_t;

CommentAnnotatedInputFileParser::CommentAnnotatedInputFileParser(
                                    std::string filename,
                                    std::string comment_symbol,
                                    functionDictionary_t known_expr_funcs)
    : InputFileParser(filename, known_expr_funcs)
    , comment_symbol_(comment_symbol)
{}

CommentAnnotatedInputFileParser::~CommentAnnotatedInputFileParser()
{}

void CommentAnnotatedInputFileParser::getProblem(Expressions::Named_t &objectives,
                                     Expressions::Named_t &constraints,
                                     DVarContainer_t &dvars) {

    for(Expressions::SingleNamed_t obj : nobjectives_) {
        if(objectives_.count(obj.first) > 0)
            objectives.insert(Expressions::SingleNamed_t(obj.first, obj.second));
    }

    for(Expressions::SingleNamed_t constr : nconstraints_) {
        if(constraints_.count(constr.first) > 0)
            constraints.insert(Expressions::SingleNamed_t(constr.first, constr.second));
    }

    for(namedDVar_t dvar : ndvars_) {
        if(dvars_.count(dvar.first) > 0)
            dvars.insert(namedDVar_t(dvar.first, dvar.second));
    }
}


void CommentAnnotatedInputFileParser::doParse() {

    std::ifstream file;
    file.open(filename_.c_str(), std::ios::in);
    if(!file)
        throw OptPilotException("CommentAnnotatedInputFileParser::doParse()",
                                "Unable to open file " + filename_);

    std::string stat_filename = filename_;
    boost::trim_right_if(stat_filename, boost::is_any_of(".in"));
    stat_filename = stat_filename.append(".stat");

    while(file) {

        std::string line;
        std::getline(file, line, '\n');

        size_t fpos = line.find("OBJECTIVE,");
        if(fpos != std::string::npos) {
            std::string name = getAttribute(line, "name");
            std::string expr = getAttribute(line, "EXPR");
            nobjectives_.insert(
                    Expressions::SingleNamed_t(
                        name, new Expressions::Expr_t(expr, known_expr_funcs_)));
        }

        fpos = line.find("CONSTRAINT,");
        if(fpos != std::string::npos) {
            std::string name = getAttribute(line, "name");
            std::string expr = getAttribute(line, "EXPR");
            nconstraints_.insert(Expressions::SingleNamed_t(
                        name, new Expressions::Expr_t(expr, known_expr_funcs_)));
        }

        fpos = line.find("DVAR,");
        if(fpos != std::string::npos) {
            std::string name = getAttribute(line, "name");
            std::string var = getAttribute(line, "VARIABLE");
            std::string lowerbound = getAttribute(line, "LOWERBOUND");
            std::string upperbound = getAttribute(line, "UPPERBOUND");

            double lb,ub;
            std::istringstream str2lb(lowerbound);
            std::istringstream str2ub(upperbound);
            str2lb >> lb;
            str2ub >> ub;

            DVar_t tmp = boost::make_tuple( var, lb, ub );
            ndvars_.insert(namedDVar_t(name, tmp));
        }


        fpos = line.find("OPTIMIZE,");
        if(fpos != std::string::npos) {

            //FIXME: ignored for now
            //string optname_ = getAttribute(line, "name");
            //string listObjectives_ = getAttribute(line, "OBJECTIVES");
            //string listConstraints_ = getAttribute(line, "CONSTRAINTS");
            //string listDVars_ = getAttribute(line, "DVARS");
        } else {

            fpos = line.find("OBJECTIVES");
            if(fpos != std::string::npos)
                objectives_ = getListAttribute(line);

            fpos = line.find("CONSTRAINTS");
            if(fpos != std::string::npos)
                constraints_ = getListAttribute(line);

            fpos = line.find("DVARS");
            if(fpos != std::string::npos)
                dvars_ = getListAttribute(line);
        }
    }
}

std::set<std::string> CommentAnnotatedInputFileParser::getListAttribute(
        std::string str) {

    size_t sstart = str.find("(") + 1;
    size_t send = str.find(")");
    std::string args = str.substr(sstart, (send-sstart));

    StringList_t attributes;
    boost::split(attributes, args, boost::is_any_of(","),
                 boost::token_compress_on);

    std::set<std::string> ret;
    for(std::string list_attribute_name : attributes) {
        boost::trim(list_attribute_name);
        ret.insert(list_attribute_name);
    }
    return ret;
}

std::string CommentAnnotatedInputFileParser::getAttribute(
        std::string str, std::string attribute) {

    // handle name attribute differently since it does not follow the default
    // "name = value" attribute style
    if(attribute.compare("name") == 0) {

        // get string on left side of colon
        StringList_t res;
        boost::split(res, str, boost::is_any_of(":"),
                     boost::token_compress_on);
        std::string ret = res[0];

        // remove comment symbols at the beginning of the line
        if(comment_symbol_ == "") {
            // trim left space
            boost::trim_left(ret);
            return ret;
        }

        // get string on right side of comment
        boost::split(res, ret, boost::is_any_of(comment_symbol_),
                     boost::token_compress_on);
        return res[1];

    } else {

        size_t fpos = str.find(attribute + "=");
        if(fpos == std::string::npos)
            throw OptPilotException("CommentAnnotatedInputFileParser::getAttribute",
                                    "Attribute " + attribute + " not found!");

        std::string ret = str.substr(fpos, std::string::npos);

        //XXX: since function arguments use the ',' delemitter we have to
        // split by searching in reversed direction the last ',' from the next
        // attribute.
        size_t end = std::string::npos;
        size_t next_attribute = ret.find('=', fpos+1);
        if(next_attribute != std::string::npos) {
            ret.erase(next_attribute, std::string::npos);
            end = ret.rfind(",");
        } else {
            end = ret.rfind(";");
            if(end == std::string::npos)
                throw OptPilotException("CommentAnnotatedInputFileParser::getAttribute",
                                        "Attribute not properly terminated!");
        }

        StringList_t res;
        ret.erase(end, std::string::npos);
        boost::split(res, ret, boost::is_any_of("="),
                     boost::token_compress_on);
        ret = res[1];

        // strip quotes
        boost::trim_left_if(ret, boost::is_any_of("\""));
        boost::trim_right_if(ret, boost::is_any_of("\""));

        return ret;
    }
}

