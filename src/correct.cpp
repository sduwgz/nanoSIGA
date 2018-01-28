#include "mole.h"
#include "runner.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::getLogger("nanoARCS.Correct"));

class Correcting : public Runner {
public:
    virtual int run(const Properties options, const Arguments& arg);
private:
    Correcting() : Runner("a:f:o:t:", boost::assign::map_list_of('a', "algorithm")('f', "overlapfile")('o', "prefix")('t', "threads")) {
        RunnerManager::instance()->install("corrct", this);
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
    if((r = checkOptions(options, arg)) == 1) {
        return r;
    }
    std::string input = args[0];
    LOG4CXX_INFO(logger, boost::format("input file is: %s") % input);
    std::string output = boost::filesystem::path(input).stem().string();
    if(options.find("prefix") != options.not_found()) {
        output = options.get< std::string >("prefix");
    }
    LOG4CXX_INFO(logger, boost::format("output file is: %s") % output);
    CorrectBuilder builder(output);
    if(!builder.build(input, options.get< double >("minscore", 15.0), output, options.get< size_t >("threads", 1))) {
        LOG4CXX_ERROR(logger, boost::format("Failed to build overlap from %s MOLECULEs") % input);
        r = -1;
    }
    return r;
}
