#include "overlap_builder.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.overlap"));

void OverlapBuilder::start(std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId, int batch) {
    //std::vector<Alignment>& alignments = *alignmentsPtr;
    int batchEnd = std::min((batch + 1) * BATCH, (int)_moleSet.size());
    for(int i = batch * BATCH + threadId; i < batchEnd; i += threads) {
        auto hitMole = _index.query(_moleSet[i]);
        if(hitMole.size() > 30) {
            for(int j = i + 1; j < _moleSet.size(); ++ j) {
                if(hitMole.find(_moleSet[j].getID()) != hitMole.end()) {
                    Alignment align = _maptool.localDPscore(_moleSet[i], _moleSet[j]);
                    //a hard threshold
                    if(align.getScore() < minScore || align.getScore() / align.size() < (minScore / 20)) {
                        continue;
                    }
                    alignmentsPtr->push_back(align);
                }
            }
        } else {
            for(int j = 0; j < _moleSet.size(); ++ j) {
                Alignment align = _maptool.localDPscore(_moleSet[i], _moleSet[j]);
                //a hard threshold
                if(align.getScore() < minScore || align.getScore() / align.size() < (minScore / 20)) {
                    continue;
                }
                alignmentsPtr->push_back(align);
            }
        }
    }
}
void OverlapBuilder::alignment(std::ofstream& overlapOutstream, int threads, double minScore, int trim=1) {
    //divide moles into batch to prevent memoery overflow
    for(int i = 0; i < _moleSet.size() / BATCH + 1; ++ i) {
        std::vector<std::vector<Alignment>> threadAlignments(threads, std::vector<Alignment>());
        boost::thread_group group;
        for(int j = 0; j < threads; ++ j) {
            std::vector<Alignment>* alignmentsPtr = &threadAlignments[j];
            group.create_thread(boost::bind(&OverlapBuilder::start, this, alignmentsPtr, minScore, threads, j, i));
        }
        group.join_all();
        //collect all alignments into alignments, and output to overlap file
        for(int j = 0; j < threads; ++ j) {
            for(auto al : threadAlignments[j]) {
                if(trim == 1) {
                    al.trimHead();
                }
                if(al.score > minScore) {
                    al.print(overlapOutstream, true);
                }
            }
        }
    }
}

bool OverlapBuilder::build(const std::string& input, double minScore, const std::string& output, size_t threads, int trim, bool reverse, bool useHash) {
    if (boost::filesystem::exists(input)) {
        LOG4CXX_WARN(logger, boost::format("Get moleculars from %s.") % input);
        std::ifstream moleInstream(input.c_str());
        MoleReader mReader(moleInstream);
        Mole m;
        while(mReader.read(m)) {
            //40 is a threshold for chimeric moleculars
            if(m.size() < 10 || m.size() > 40) continue;
            _moleSet.push_back(m);
            if (reverse) {
                Mole reMole = m.reverseMole();
                _moleSet.push_back(reMole);
            }
        }
        if (_moleSet.size() == 0) {
            LOG4CXX_WARN(logger, "No mole is in moleSet.");
            return false;
        } else {
            LOG4CXX_INFO(logger, boost::format("%s moles have been inited.") % _moleSet.size());
        } 
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % input);
        return false;
    }
    std::ofstream overlapOutstream(output.c_str());
    LOG4CXX_INFO(logger, boost::format("Alignment using %d threads.") % threads);
    //build index
    if(useHash) {
        buildIndex();
        LOG4CXX_INFO(logger, boost::format("Hash table has been built."));
    }
    //instantiate a DP solver
    alignment(overlapOutstream, threads, minScore, trim);
    //TODO
    return true;
}
