#include "overlap_builder.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.overlap"));

bool OverlapBuilder::build(const std::string& input, double minScore, const std::string& output, size_t threads, int trim, int reverseLabel) const {
    std::vector<Mole> moleSet;
    if (boost::filesystem::exists(input)) {
        std::ifstream moleInstream(input.c_str());
        MoleReader mReader(moleInstream);
        Mole m;
        while(mReader.read(m)) {
            moleSet.push_back(m);
            if (reverseLabel == 1) {
                Mole reMole = m.reverseMole();
                moleSet.push_back(reMole);
            }
        }
        int moleNumber = moleSet.size();
        if (moleNumber == 0) {
            LOG4CXX_WARN(logger, "no mole is in moleSet");
            return false;
        } else {
            LOG4CXX_INFO(logger, boost::format("%s moles have been inited.") % moleNumber);
        } 
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % input);
        return false;
    }
    std::ofstream overlapOutstream(output.c_str());
    for(int i = 0; i < moleSet.size(); ++ i) {
        for(int j = i + 1; j < moleSet.size(); ++ j) {
            Alignment ret = _maptool.localDPscore(moleSet[i], moleSet[j]);
            if(trim == 1) {
                ret.trim();
            }
            if(ret.score > minScore) {
                ret.print(overlapOutstream, false);
            }
        }
    }
    return true;
}
