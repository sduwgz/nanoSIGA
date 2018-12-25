#include "overlap_builder.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.overlap"));

void start(const Map* maptool, const Index* index, const std::vector<Mole>& moleSet, std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId, int batch) {
    //const std::vector<Mole>& moleSet = *moleSetPtr;
    std::vector<Alignment>& alignments = *alignmentsPtr;
    int dpCount = 0;
    int batchSize = std::min((batch + 1) * BATCH, (int)moleSet.size());
    for(int i = batch * BATCH + threadId; i < batchSize; i += threads) {
        if(i % threads == 1 && ((i / threads) % 100 == 0))
            std::cout << i / threads << std::endl;
        if(moleSet[i].size() < 6) continue;
        auto hitMole = index->query(moleSet[i]);
        if(hitMole.size() > 30) {
            for(int j = i + 1; j < moleSet.size(); ++ j) {
                if(hitMole.find(moleSet[j].getID()) != hitMole.end()) {
                    Alignment align = maptool->localDPscore(moleSet[i], moleSet[j]);
                    //a hard threshold
                    if(align.getScore() < minScore || align.getScore() / align.size() < (minScore / 20)) {
                        continue;
                    }
                    alignments.push_back(align);
                }
            }
        } else {
            for(int j = 0; j < moleSet.size(); ++ j) {
                Alignment align = maptool->localDPscore(moleSet[i], moleSet[j]);
                //a hard threshold
                if(align.getScore() < minScore || align.getScore() / align.size() < (minScore / 20)) {
                    continue;
                }
                alignments.push_back(align);
            }
        }
    }
}
void alignment(const std::string parameterFile, const std::vector<Mole>& moleSet, std::ofstream& overlapOutstream, int threads, double minScore, int trim=1, bool useHash=false) {
    //build index
    Index *index = new Index(3);
    if(useHash) {
        index->build(moleSet);
        LOG4CXX_INFO(logger, boost::format("Hash table has been built."));
    }
    Map* maptool = Map::instance(parameterFile);
    //const std::vector<Mole>* moleSetPtr = &moleSet;
    int n = moleSet.size() / BATCH + 1;
    for(int batch = 0; batch < n; ++ batch) {
        std::vector<Alignment> alignments;
        std::vector<std::vector<Alignment>> threadAlignments(threads, std::vector<Alignment>());
        boost::thread_group group;
        LOG4CXX_INFO(logger, boost::format("Alignment using %d threads.") % threads);
        for(int i = 0; i < threads; ++ i) {
            std::vector<Alignment>* alignmentsPtr = &threadAlignments[i];
            group.create_thread(boost::bind(start, maptool, index, moleSet, alignmentsPtr, minScore, threads, i, batch));
        }
        group.join_all();
        for(int i = 0; i < threads; ++ i) {
            for(auto al : threadAlignments[i]) {
                alignments.push_back(al);
            }
        }
        for(int i = 0; i < alignments.size(); ++ i) {
            Alignment ret = alignments[i];
            if(trim == 1) {
                ret.trimHead();
                //ret.trimTail();
            }
            if(ret.score > minScore) {
                ret.print(overlapOutstream, true);
            }
        }
    }
}

bool OverlapBuilder::build(const std::string& input, double minScore, const std::string& output, size_t threads, int trim, bool reverseLabel, bool useHash) const {
    std::vector<Mole> moleSet;
    if (boost::filesystem::exists(input)) {
        std::ifstream moleInstream(input.c_str());
        MoleReader mReader(moleInstream);
        Mole m;
        while(mReader.read(m)) {
            if(m.size() < 10 || m.size() > 40) continue;
            moleSet.push_back(m);
            if (reverseLabel) {
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
    //std::vector<Alignment> alignments;
    alignment(_parameterFile, moleSet, overlapOutstream, threads, minScore, 1, useHash);
    //TODO
    return true;
}
