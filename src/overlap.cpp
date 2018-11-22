#include "runner.h"
#include "mole.h"
#include "map.h"
#include "overlap_builder.h"

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.Overlap"));

class Overlapping : public Runner {
public:
    virtual int run(const Properties options, const Arguments& arg);
private:
    Overlapping() : Runner("m:p:o:t:SR", boost::assign::map_list_of('m',"minscore")('p', "parameter")('o', "prefix")('t', "threads")('S', "hash")('R', "reverse")) {
        RunnerManager::instance()->install("overlap", this);
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
                "nanoARCS overlap [OPTION] ... READSFILE\n"
                "Compute pairwise overlap between all the molecule in READS\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -t, --threads=NUM                use NUM threads to construct the index (default: 1)\n"
                "      -m, --min-score=f                minimum aligment score required between two reads (default: 5)\n"
                "      -S, --hash                       use hash strategy to speed up\n"
                "      -R, --reverse                    reverse mole as a new mole\n"
                "      -o, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
                "\n"
                ) << std::endl;
        return 256;
    }

    static Overlapping _runner;
};
Overlapping Overlapping::_runner;

int Overlapping::run(const Properties options, const Arguments& arg) {
    int r = 0;
    if((r = checkOptions(options, arg)) == 1) {
        return r;
    }
    std::string parameter_file = "parameters.ini";
    if(options.find("parameter") != options.not_found()) {
        parameter_file = options.get< std::string > ("parameter");
    }

    std::string input = arg[0];
    LOG4CXX_INFO(logger, boost::format("input file is: %s") % input);
    std::string output = boost::filesystem::path(input).stem().string();
    if(options.find("prefix") != options.not_found()) {
        output = options.get< std::string >("prefix");
    }
    LOG4CXX_INFO(logger, boost::format("output file is: %s") % output);

    bool useHash = false, useReverse = false;
    if(options.find("hash") != options.not_found()) {
        useHash = true;
        LOG4CXX_INFO(logger, boost::format("hash strategy has been used to speed up."));
    }
    if(options.find("reverse") != options.not_found()) useReverse = true;
    OverlapBuilder builder(parameter_file, output);
    if(!builder.build(input, options.get< double >("minscore", 15.0), output, options.get< size_t >("threads", 1), 1, useReverse, useHash)) {
        LOG4CXX_ERROR(logger, boost::format("Failed to build overlap from %s MOLECULEs") % input);
        r = -1;
    }
    return r;
}
