#include "overlap_builder.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.overlap"));

void start(const Map* maptool, const std::vector<Mole>* moleSetPtr, std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId) {
    const std::vector<Mole>& moleSet = *moleSetPtr;
    std::vector<Alignment>& alignments = *alignmentsPtr;
    for(int i = threadId; i < moleSet.size(); i += threads) {
        for(int j = i + 1; j < moleSet.size(); ++ j) {
            Alignment align = maptool->localDPscore(moleSet[i], moleSet[j]);
            //a hard threshold
            //if(align.score < minScore || align.score / align.alignedMole1.size() < (minScore / 10)) {
            if(align.score < minScore) {
                continue;
            }
            alignments.push_back(align);
        }
    }
}
void alignment(const std::string parameterFile, const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, int threads, double minScore) {
    const std::vector<Mole>* moleSetPtr = &moleSet;
    std::vector<std::vector<Alignment>> threadAlignments(threads, std::vector<Alignment>());
    Map* maptool = Map::instance(parameterFile);
    boost::thread_group group;
    LOG4CXX_INFO(logger, boost::format("Alignment using %d threads.") % threads);
    for(int i = 0; i < threads; ++ i) {
        std::vector<Alignment>* alignmentsPtr = &threadAlignments[i];
        group.create_thread(boost::bind(start, maptool, moleSetPtr, alignmentsPtr, minScore, threads, i));
    }
    group.join_all();
    for(int i = 0; i < threads; ++ i) {
        for(auto al : threadAlignments[i]) {
            alignments.push_back(al);
        }
    }
}

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
    std::vector<Alignment> alignments;
    alignment(_parameterFile, moleSet, alignments, threads, minScore);
    //TODO
    for(int i = 0; i < alignments.size(); ++ i) {
        Alignment ret = alignments[i];
        if(trim == 1) {
            ret.trimHead();
            //ret.trimTail();
        }
        if(ret.score > minScore) {
            ret.print(overlapOutstream, false);
        }
    }
    return true;
}
