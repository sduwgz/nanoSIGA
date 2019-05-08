#include "mole.h"
#include "runner.h"
#include "refine_builder.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/assign.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.Refine"));

class Refining : public Runner {
public:
    virtual int run(const Properties options, const Arguments& arg);
private:
    Refining() : Runner("i:o:p:", boost::assign::map_list_of('i', "contigFile")('o', "prefix")('p', "parameter")) {
        RunnerManager::instance()->install("refine", this);
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
                "nanoARCS refine [OPTION] ... READSFILE\n"
                "Refinement contig by alignment.\n"
                "\n"
                "      -h, --help                       display this help and exit\n"
                "\n"
                "      -i, --contigfile=file            contig file[BNX]\n"
                "      -p, --parameterfile=file         parameter file[INI]\n"
                "      -o, --prefix=PREFIX              write refined result to file using PREFIX instead of prefix of READSFILE\n"
                "\n"
                ) << std::endl;
        return 256;
    }

    static Refining _runner;
};

Refining Refining::_runner;
int Refining::run(const Properties options, const Arguments& args) {
    int r = 0;
    if((r = checkOptions(options, args)) == 1) {
        return r;
    }
    std::string contigFile = options.get< std::string >("contigFile");
    std::string moleFile = args[0];
    LOG4CXX_INFO(logger, boost::format("contig file is: %s") % contigFile);
    std::string output = boost::filesystem::path(contigFile).stem().string() + "refined";
    if(options.find("prefix") != options.not_found()) {
        output = options.get< std::string >("prefix");
    }
    LOG4CXX_INFO(logger, boost::format("output file is: %s") % output);
    std::string parameter_file = "parameters.ini";
    if(options.find("parameter") != options.not_found()) {
        parameter_file = options.get< std::string > ("parameter");
    }
    
    RefineBuilder builder(parameter_file);
    if(!builder.build(contigFile, moleFile, output)) {
        LOG4CXX_ERROR(logger, boost::format("Failed to refine contigs.") % contigFile);
        r = -1;
    }
    return r;
}
