#include "refine_builder.h"

#include <iostream>
#include <fstream>
#include <queue>
#include <utility>
#include <map>
#include <unordered_map>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.refine"));
void mergeExpand(const Alignment& align, const Mole& m1, const Mole& m2, std::vector<int>& newDistance, bool flag) {
    //debug
    std::cout << "merge: " << m1.getID() << " " << m2.getID() << std::endl;
    for(int d : m1.getData()) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
    for(int d : m2.getData()) {
        std::cout << d << " ";
    }
    std::cout << std::endl;

    //Case 1: 
    if(align.mole1Start != 0 && align.mole1End != m1.size()) {
        if(flag) {
            for(int k = 0; k < align.mole1Start; ++ k) {
                newDistance.push_back(m1.getInterval(k));
            }
            for(int k = 0; k < align.alignedMole1.size(); ++ k) {
                if(align.alignedMole1[k].size() == 1 && align.alignedMole2[k].size() == 1)
                    newDistance.push_back((align.alignedMole1[k][0] + align.alignedMole2[k][0]) / 2);
                else if(align.alignedMole1[k].size() != 1)
                    newDistance.insert(newDistance.end(), align.alignedMole1[k].begin(), align.alignedMole1[k].end());
                else if(align.alignedMole2[k].size() != 1)
                    newDistance.insert(newDistance.end(), align.alignedMole2[k].begin(), align.alignedMole2[k].end());
            }
            for(int k = align.mole1End; k < m1.size(); ++ k)
                newDistance.push_back(m1.getInterval(k));
        } else {
            for(int k = 0; k < m1.size(); ++ k)
                newDistance.push_back(m1.getInterval(k));
        }
    }
    //Case 2:
    if(align.mole1Start == 0 && align.mole1End != m1.size()) {

        for(int k = 0; k < align.mole2Start; ++ k) {
            newDistance.push_back(m2.getInterval(k));
        }
        for(int k = 0; k < align.alignedMole1.size(); ++ k) {
            if(align.alignedMole1[k].size() == 1 && align.alignedMole2[k].size() == 1)
                newDistance.push_back((align.alignedMole1[k][0] + align.alignedMole2[k][0]) / 2);
            else if(align.alignedMole1[k].size() != 1)
                newDistance.insert(newDistance.end(), align.alignedMole1[k].begin(), align.alignedMole1[k].end());
            else if(align.alignedMole2[k].size() != 1)
                newDistance.insert(newDistance.end(), align.alignedMole2[k].begin(), align.alignedMole2[k].end());
        }
        for(int k = align.mole1End; k < m1.size(); ++ k)
            newDistance.push_back(m1.getInterval(k));
    }
    //Case 3:
    if(align.mole1Start != 0 && align.mole1End == m1.size()) {
        for(int k = 0; k < align.mole1Start; ++ k) {
            newDistance.push_back(m1.getInterval(k));
        }
        for(int k = 0; k < align.alignedMole1.size(); ++ k) {
            if(align.alignedMole1[k].size() == 1 && align.alignedMole2[k].size() == 1)
                newDistance.push_back((align.alignedMole1[k][0] + align.alignedMole2[k][0]) / 2);
            else if(align.alignedMole1[k].size() != 1)
                newDistance.insert(newDistance.end(), align.alignedMole1[k].begin(), align.alignedMole1[k].end());
            else if(align.alignedMole2[k].size() != 1)
                newDistance.insert(newDistance.end(), align.alignedMole2[k].begin(), align.alignedMole2[k].end());
        }
        for(int k = align.mole2End; k < m2.size(); ++ k)
            newDistance.push_back(m2.getInterval(k));
    }
    std::cout << "merge" << std::endl;
    for(int d : newDistance) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
}

bool RefineBuilder::build(const std::string& input, const std::string& moleFile, const std::string& output) const {
    // read mole
    std::vector<Mole> moleSet;
    if(boost::filesystem::exists(moleFile)) {
        std::ifstream moleInstream(moleFile.c_str());
        int count = -1;
        MoleReader mReader(moleInstream);
        Mole m;
        while(mReader.read(m)) {
            if(m.size() < 10 || m.size() > 40) continue;
            moleSet.push_back(m);
        }
        int moleNumber = moleSet.size();
        if (moleNumber == 0) {
            LOG4CXX_WARN(logger, "no mole is in moleSet");
            return false;
        } else {
            LOG4CXX_WARN(logger, boost::format("%s moles have been inited.") % moleNumber);
        }
        moleInstream.close();
    }
    // read contig
    std::vector<Mole> contigSet;
    if(boost::filesystem::exists(input)) {
        std::ifstream contigInstream(input.c_str());
        int count = -1;
        MoleReader mReader(contigInstream);
        Mole contig;
        while(mReader.read(contig)) {
            if(contig.size() < 10) continue;
            contigSet.push_back(contig);
        }
        int contigNumber = contigSet.size();
        if (contigNumber == 0) {
            LOG4CXX_WARN(logger, "no contig is in contig file.");
            return false;
        } else {
            LOG4CXX_WARN(logger, boost::format("%s contigs have been inited.") % contigNumber);
        }
        contigInstream.close();

    }
    // merge contigs
    // sort by length
    std::sort(contigSet.begin(), contigSet.end(), [](const Mole& m1, const Mole& m2) {return m1.size() > m2.size();});
    std::set<int> usedContig;
    Map* maptool = Map::instance(_parameter_file);
    std::string contigAlign = output + ".align";
    std::ofstream os(contigAlign.c_str());
    for(int i = 0; i < contigSet.size(); ++ i) {
        std::cout << i << std::endl;
        if(usedContig.find(i) != usedContig.end()) continue;
        bool flag = true;
        while(flag) {
            flag = false;
            for(int j = i + 1; j < contigSet.size(); ++ j) {
                if(usedContig.find(j) != usedContig.end()) continue;
                Mole recontig = contigSet[j].reverseMole();
                Alignment align = maptool->localDPscore(contigSet[i], contigSet[j]);
                Alignment realign = maptool->localDPscore(contigSet[i], recontig);
                if(align.score >= 10) {
                    usedContig.insert(j);
                    std::vector<int> newDistance;
                    align.print(os, true);
                    mergeExpand(align, contigSet[i], contigSet[j], newDistance, true);
                    contigSet[i].setData(newDistance);
                    flag = true;
                } else if(realign.score >= 10) {
                    usedContig.insert(j);
                    std::vector<int> newDistance;
                    realign.print(os, true);
                    mergeExpand(realign, contigSet[i], recontig, newDistance, true);
                    contigSet[i].setData(newDistance);
                    flag = true;
                }
            }
        }
    }
    std::cout << "refined" << std::endl;
    std::vector<Mole> refinedContigSet;
    for(int i = 0; i < contigSet.size(); ++ i) {
        if(usedContig.find(i) == usedContig.end()) 
            refinedContigSet.push_back(contigSet[i]);
    }

    // mole stitch
    std::map<std::string, std::set<std::string>> moleContigs;
    for(int i = 0; i < moleSet.size(); ++ i) {
        std::cout << "mole stitch: " << i << std::endl;
        for(int j = 0; j < refinedContigSet.size(); ++ j) {
            Mole recontig = refinedContigSet[j].reverseMole();
            Alignment align = maptool->localDPscore(moleSet[i], refinedContigSet[j]);
            Alignment realign = maptool->localDPscore(moleSet[i], recontig);
            if(align.score >= 15 && (align.mole2Start == 0 || align.mole2End == refinedContigSet[j].size())) {
                moleContigs[align.mole1Id].insert(align.mole2Id);
                std::cout << align.mole1Id << " " << j + 1 << std::endl;
                align.print(os, true);
            } else if(realign.score >= 15 && (realign.mole2Start == 0 || realign.mole2End == refinedContigSet[j].size())) {
                moleContigs[realign.mole1Id].insert(realign.mole2Id);
                std::cout << align.mole1Id << " " << j + 1 << std::endl;
                realign.print(os, true);
            }
        }
        if(moleContigs[moleSet[i].getID()].size() >= 2) {
            std::cout << "connect" << std::endl;
            
        }
    }
    
    // mole extension
    for(int i = 0; i < refinedContigSet.size(); ++ i) {
        std::cout << "mole refine: " << i << std::endl;
        for(int j = 0; j < moleSet.size(); ++ j) {
            Mole remole = moleSet[j].reverseMole();
            Alignment align = maptool->localDPscore(refinedContigSet[i], moleSet[j]);
            Alignment realign = maptool->localDPscore(refinedContigSet[i], remole);
            if(align.score >= 1000) {
                std::vector<int> newDistance;
                align.print(os, true);
                mergeExpand(align, refinedContigSet[i], moleSet[j], newDistance, false);
                refinedContigSet[i].setData(newDistance);
            } else if(realign.score >= 1000) {
                std::vector<int> newDistance;
                realign.print(os, true);
                mergeExpand(realign, refinedContigSet[i], remole, newDistance, false);
                refinedContigSet[i].setData(newDistance);
            }
        }
    }
    std::string outputbnx = output + ".bnx";
    std::string outputcmap = output + ".cmap";
    std::ofstream bnxos(outputbnx.c_str());
    std::ofstream cmapos(outputcmap.c_str());
    MoleWriter mwriter1(bnxos);
    MoleWriter mwriter2(cmapos);
    for(auto contig : refinedContigSet) {
        for(auto it : contig.getData()) {
            std::cout << it << " ";
        }
        std::cout << std::endl;
        mwriter1.writebnx(contig);
        mwriter2.writecmap(contig);
    }
    bnxos.close();
    cmapos.close();
    os.close();
    return true;
}
