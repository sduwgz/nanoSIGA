#include "runner.h"
#include "mole.h"
#include "map.h"
#include "cluster_builder.h"

#include <vector>
#include <map>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/assign.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.cluster"));

class Clustering : public Runner {
public:
    virtual int run(const Properties options, const Arguments& arg);
private:
    Clustering() : Runner("m:s:o:", boost::assign::map_list_of('m', "mincluster")('s', "minscore")('o', "prefix")) {
        RunnerManager::instance()->install("cluster", this);
    }
    int checkOptions(const Properties options, const Arguments& arg) {
        if(options.find("h") != options.not_found() || arg.size() == 0) {
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
                "      -m, --mincluster=NUM             size of min clusters have NUM molecule(default: 2)\n"
                "      -s, --minscore=S                 size of min clusters have NUM molecule(default: 2)\n"
                "      -o, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
                "\n"
                ) << std::endl;
        return 256;
    }

    static Clustering  _runner;
};
Clustering Clustering::_runner;

int Clustering::run(const Properties options, const Arguments& arg) {
    int r = 0;
    if((r = checkOptions(options, arg)) == 1) {
        return r;
    }
    std::string input = arg[0];
    LOG4CXX_INFO(logger, boost::format("input file is: %s") % input);
    std::string output = boost::filesystem::path(input).stem().string();
    if(options.find("prefix") != options.not_found()) {
        output = options.get< std::string >("prefix");
    }
    LOG4CXX_INFO(logger, boost::format("output file is: %s") % output);
    ClusterBuilder builder(output);
    if(!builder.build(input, options.get< double >("minscore", 23.0), options.get< double >("mincluster", 2), output)) {
        LOG4CXX_ERROR(logger, boost::format("Failed to build overlap from %s MOLECULEs") % input);
        r = -1;
    }
    return r;
}
