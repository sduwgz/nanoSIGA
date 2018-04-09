#include "mole.h"
#include "runner.h"
#include "correct_builder.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/assign.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.Correct"));

class Correcting : public Runner {
public:
    virtual int run(const Properties options, const Arguments& arg);
private:
    Correcting() : Runner("a:f:o:t:p:", boost::assign::map_list_of('a', "algorithm")('f', "overlapfile")('o', "prefix")('t', "threads")('p', "parameter")) {
        RunnerManager::instance()->install("correct", this);
    }

    int checkOptions(const Properties options, const Arguments& args) const {
        if(options.find("h") != options.not_found() || args.size() == 0) {
            printHelps();
            return 1;
        }
        return 0;
    }
    int printHelps() const {
        std::cout << boost::format(
                "nanoARCS correct [OPTION] ... READSFILE\n"
                "Correct and vote for a center for each clusters.\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -a, --algorithm=A                use A correct function (default: vote)\n"
                "      -f, --overlapfile=file           input file\n"
                "      -t, --threads=NUM                use NUM threads to construct the index (default: 1)\n"
                "      -o, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
                "\n"
                ) << std::endl;
        return 256;
    }

    static Correcting _runner;
};

Correcting Correcting::_runner;
int Correcting::run(const Properties options, const Arguments& args) {
    int r = 0;
    if((r = checkOptions(options, args)) == 1) {
        return r;
    }
    std::string input = args[0];
    std::string moleFile = args[1];
    LOG4CXX_INFO(logger, boost::format("input file is: %s") % input);
    std::string output = boost::filesystem::path(input).stem().string();
    if(options.find("prefix") != options.not_found()) {
        output = options.get< std::string >("prefix");
    }
    LOG4CXX_INFO(logger, boost::format("output file is: %s") % output);
    std::string parameter_file = "parameters.ini";
    if(options.find("parameter") != options.not_found()) {
        parameter_file = options.get< std::string > ("parameter");
    }
    
    CorrectBuilder builder(parameter_file, output);
    if(!builder.build(moleFile, input, options.get< double >("minscore", 15.0), output, options.get< size_t >("threads", 1))) {
        LOG4CXX_ERROR(logger, boost::format("Failed to correct form %s alignments") % input);
        r = -1;
    }
    return r;
}
