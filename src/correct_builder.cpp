#include "correct_builder.h"

#include <iostream>
#include <fstream>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#include <log4cxx/logger.h>

typedef std::vector<int> Center;
typedef std::map<std::string, Mole> MoleSet;
typedef std::vector<std::vector<Index>> IndexesSet;
typedef std::map<int, std::vector<Center>> CenterSet;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.correct"));

Index center(const std::vector<Index>& indexes) {
    Index centerIndex("-1", 0, 0);
    for(auto index : indexes) {
        if(index.endSite - index.startSite > centerIndex.endSite - centerIndex.startSite)
            centerIndex = index;
    }
    return centerIndex;
}
int clacDistance(std::vector<int> distance) {
    int tempDis = 0;
    if(distance.size() == 2)
        tempDis = std::accumulate(distance.begin(), distance.end(), 0.0) / distance.size();
    else if(distance.size() == 3) {
        int ta = abs(distance[0] - distance[1]);
        int tb = abs(distance[2] - distance[1]);
        int tc = abs(distance[2] - distance[0]);
        if(ta < tb && ta < tc) {
            tempDis = (distance[0] + distance[1]) / 2;
        } else if(tb < ta && tb < tc)
            tempDis = (distance[2] + distance[1]) / 2;
        else
            tempDis = (distance[2] + distance[0]) / 2;
    } else {
        std::sort(distance.begin(), distance.end());
        int maxCluster = 0, maxIndex = 0;
        for(int i = 0; i < distance.size(); ++ i) {
            int count = 0;
            for(int j = i; j < distance.size(); ++ j) {
                if(distance[j] - distance[i] < 500) {
                    count ++;
                }
            }
            if(maxCluster < count) {
                maxCluster = count;
                maxIndex = i;
            }
        }
        if(maxCluster == 1) {
            maxIndex = 0, maxCluster = distance.size();
        }
        tempDis = std::accumulate(distance.begin() + maxIndex, distance.begin() + maxIndex + maxCluster, 0.0) / (maxCluster);
    }
    return tempDis;
}
Center vote(std::vector<Alignment>& alignments) {
    for(Alignment al : alignments) {
        al.print(std::cout, true);
    }
    std::vector<std::vector<int>> clusterMole;
    std::vector<int> mole;
    for(Alignment &al : alignments) {
        al.trimTail();
        if(al.alignedMole1.size() > mole.size()) {
            mole.clear();
            for(Fragment frag : al.alignedMole1) {
                for(int dis : frag) {
                    mole.push_back(dis);
                }
            }
        }
    }
    if(mole.size() == 0) {
        LOG4CXX_WARN(logger, boost::format("template mole is not found"));
    }
    clusterMole.push_back(mole);
    for(Alignment al : alignments) {
        mole.clear();
        for(Fragment frag : al.alignedMole2) {
            for(int dis : frag) {
                mole.push_back(dis);
            }
        }
        clusterMole.push_back(mole);
    }
    LOG4CXX_DEBUG(logger, boost::format("cluster size %d, center length %d") % clusterMole.size() % clusterMole[0].size());

    int clusterSize = alignments.size() + 1;
    std::vector<std::map<int, int>> correspondAlignment;
    for(int i = 0; i < clusterSize - 1; ++ i) {
        int p = 0, q = 0;
        std::map<int, int> correspondSite;
        correspondSite[p] = q;
        for(int j = 0; j < alignments[i].alignedMole2.size(); ++ j) {
            p += alignments[i].alignedMole1[j].size();
            q += alignments[i].alignedMole2[j].size();
            correspondSite[p] = q;
        }
        for(int j = 0; j < clusterMole[0].size() + 1; ++ j) {
            if(j < p && correspondSite.find(j) == correspondSite.end()) {
                correspondSite[j] = -1;
            }
        }
        correspondAlignment.push_back(correspondSite);
    }
    LOG4CXX_DEBUG(logger, boost::format("correspondAlignment size %d") % correspondAlignment.size());
    std::vector<int> preSite(clusterSize, 0);
    std::vector<int> nextSite(clusterSize, 0);
    Center voteCenter;
    int count = 0;
    while(count < clusterMole[0].size()) {
        count ++;
        nextSite[0] = count;
        for(int i = 1; i < clusterSize; ++ i) {
            nextSite[i] = correspondAlignment[i - 1][count];
        }
        int missMatch = std::count(nextSite.begin(), nextSite.end(), -1);
        int cutOff = std::count(nextSite.begin(), nextSite.end(), 0);
        if(cutOff > clusterSize / 2 && clusterSize - cutOff < 2) {
            break;
        }
        if(missMatch >= (clusterSize - cutOff) / 2) {
            continue;
        }
        
        LOG4CXX_DEBUG(logger, boost::format("#site %d, missMatch %d, cutOff %d") % count % missMatch % cutOff);
        std::vector<int> betweenSiteMount;
        for(int i = 0; i < preSite.size(); ++ i) {
            betweenSiteMount.push_back(nextSite[i] - preSite[i]);
        }
        int voteSiteMount = -1;
        int tempMount = -1;
        for(int i = 1; i < 4; ++ i) {
            if(tempMount <= std::count(betweenSiteMount.begin(), betweenSiteMount.end(), i)) {
                tempMount = std::count(betweenSiteMount.begin(), betweenSiteMount.end(), i); 
                voteSiteMount = i;
             }
        }
        std::vector<int> distance;
        std::vector<int> midSite1;
        std::vector<int> midSite2;
        LOG4CXX_DEBUG(logger, boost::format("vote site is %d") % voteSiteMount);
        if(voteSiteMount == 1) {
            for(int i = 0; i < clusterSize; ++ i) {
                if(nextSite[i] != 0 && nextSite[i] != -1) {
                    int tempSum = accumulate(clusterMole[i].begin() + preSite[i], clusterMole[i].begin() + nextSite[i], 0);
                    if(tempSum > 0) {
                        distance.push_back(tempSum);
                    }
                }
            }
        }
        if(voteSiteMount == 2) {
            for(int i = 0; i < clusterSize; ++ i) {
                if(nextSite[i] != 0 && nextSite[i] != -1) {
                    int tempSum = accumulate(clusterMole[i].begin() + preSite[i], clusterMole[i].begin() + nextSite[i], 0);
                    if(tempSum > 0) {
                        distance.push_back(tempSum);
                    }
                }
                if(nextSite[i] - preSite[i] == 2) {
                    if(clusterMole[i][preSite[i]] > 0) {
                        midSite1.push_back(clusterMole[i][preSite[i]]);
                    }
                }
            }
        }
        if(voteSiteMount == 3) {
            for(int i = 0; i < clusterSize; ++ i) {
                if(nextSite[i] != 0 && nextSite[i] != -1) {
                    int tempSum = accumulate(clusterMole[i].begin() + preSite[i], clusterMole[i].begin() + nextSite[i], 0);
                    if(tempSum > 0) {
                        distance.push_back(tempSum);
                    }
                }
                if(nextSite[i] - preSite[i] == 3) {
                    if(clusterMole[i][preSite[i]] > 0) {
                        midSite1.push_back(clusterMole[i][preSite[i]]);
                    }
                    if(clusterMole[i][preSite[i] + 1] > 0) {
                        midSite2.push_back(clusterMole[i][preSite[i] + 1]);
                    }
                }
            }
        }
        for(int d : distance) {
            std::cout << d << " ";
        }
        std::cout << std::endl;
        for(int d : midSite1) {
            std::cout << d << " ";
        }
        std::cout << std::endl;
        for(int d : midSite2) {
            std::cout << d << " ";
        }
        std::cout << std::endl;
        int s1 = 0, s2 = 0;
        if(midSite1.size() != 0) {
            s1 = std::accumulate(midSite1.begin(), midSite1.end(), 0.0) / midSite1.size();
            voteCenter.push_back(s1);
        }
        if(midSite2.size() != 0) {
            s2 = std::accumulate(midSite2.begin(), midSite2.end(), 0.0) / midSite2.size();
            voteCenter.push_back(s2);
        }
        int tempDis = clacDistance(distance) - s1 - s2;
        voteCenter.push_back(std::max(tempDis , 1));
        for(int i = 0; i < preSite.size(); ++ i) {
            if(nextSite[i] != -1) {
                preSite[i] = nextSite[i];
            } else {
                clusterMole[i][preSite[i]] -= tempDis;
            }
            nextSite[i] = 0;
        }
    }
    for(int d : voteCenter) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
    return voteCenter;
}
void divide(const std::vector<Alignment>& alignments, std::map<int, std::vector<Alignment>>& dividedAlignments) {
    int clusterSize = alignments.size();
    for(int i = 0; i < clusterSize; ++ i) {
        int templateStart = alignments[i].mole1Start;
        int othersStart = alignments[i].mole2Start;
        dividedAlignments[templateStart].push_back(alignments[i]);
    }
}
void alignment(const Map* maptool, const std::map<std::string, Mole>& moleSet, const std::vector<Index>& indexes, const Index& centerIndex, std::vector<Alignment>& alignments) {
    if(moleSet.find(centerIndex.moleId) != moleSet.end()) {
        auto it = moleSet.find(centerIndex.moleId);
        Mole centerMole = it->second;
        std::string newId = centerIndex.moleId + "|";
        newId += std::to_string(centerIndex.startSite);
        centerMole.setID(newId);
        centerMole.shear(centerIndex.startSite);
        for(Index index : indexes) {
            if(index.moleId == centerIndex.moleId) continue;
            LOG4CXX_DEBUG(logger, boost::format("Alignment between index %s and center %s.") % index.toString() % centerIndex.toString());
            if(moleSet.find(index.moleId) != moleSet.end()) {
                it = moleSet.find(index.moleId);
                Mole m2 = it->second;
                m2.shear(index.startSite);
                newId = index.moleId + "|";
                newId += std::to_string(index.startSite);
                m2.setID(newId);
                Alignment align = maptool->localDPscore(centerMole, m2);
                if(align.score <= 15) {
                    continue;
                }
                alignments.push_back(align);
            } else {
                LOG4CXX_WARN(logger, boost::format("Index %s is not in cluster.") % index.toString());
            }
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("center %s is not in cluster.") % centerIndex.toString());
    }
}
//void start(const Map* maptool, const MoleSet* moleSetPtr, const IndexesSet* indexesSetPtr, CenterSet* centerSetPtr, int threadNum, int threadId) {
void start(const Map* maptool, const MoleSet& moleSet, const IndexesSet& indexesSet, CenterSet* centerSetPtr, int threadNum, int threadId) {
    //const MoleSet& moleSet = *moleSetPtr;
    //const IndexesSet& indexesSet = *indexesSetPtr;
    CenterSet& centerSet = *centerSetPtr;
    for(int i = threadId; i < indexesSet.size(); i += threadNum) {
        auto indexes = indexesSet[i];
        Index centerIndex = center(indexes);
        std::vector<Alignment> alignments;
        alignment(maptool, moleSet, indexes, centerIndex, alignments);
        std::map<int, std::vector<Alignment>> dividedAlignments;
        divide(alignments, dividedAlignments);
        for(auto it : dividedAlignments) {
            Center voteCenter = vote(it.second);
            centerSet[i].push_back(voteCenter);
        }
    }
}
bool CorrectBuilder::build(const std::string& moleFile, const std::string& clusterFile, double minScore, const std::string& output, int thread) const {
    MoleSet moleSet;
    if (boost::filesystem::exists(moleFile)) {
        std::ifstream moleInstream(moleFile.c_str());
        MoleReader mReader(moleInstream);
        Mole m;
        int reverseLabel = 0;
        while(mReader.read(m)) {
            moleSet[m.getID()] = m;
            if (reverseLabel == 1) {
                Mole reMole = m.reverseMole();
                moleSet[reMole.getID()] = reMole;
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
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % moleFile);
        return false;
    }
    if (boost::filesystem::exists(clusterFile)) {
        std::vector<Index> indexes;
        IndexesSet indexesSet;
        std::ifstream clusterInstream(clusterFile.c_str());
        std::ofstream resOutstream(output.c_str());
        std::string line;
        while(std::getline(clusterInstream, line)) {
            indexes.clear();
            boost::trim(line);
            std::vector<std::string> indexVec;
            boost::split(indexVec, line, boost::is_any_of(" "), boost::token_compress_on);
            for(std::string ind : indexVec) {
                std::vector<std::string> indexString;
                LOG4CXX_DEBUG(logger, boost::format("an index %s.") % ind);
                boost::split(indexString, ind, boost::is_any_of("|"), boost::token_compress_on);
                indexes.emplace_back(indexString[0], boost::lexical_cast<int>(indexString[1]), moleSet[indexString[0]].size());
            } 
            indexesSet.push_back(indexes);
        }
        CenterSet centerSet;
        boost::thread_group group;
        std::vector<CenterSet> threadCenterSet(thread);
        Map* maptool = Map::instance(_parameterFile);
        const MoleSet* moleSetPtr = &moleSet;
        const IndexesSet* indexesSetPtr = &indexesSet;
        for(int i = 0; i < thread; ++ i) {
            CenterSet* centerSetPtr = &threadCenterSet[i];
            //group.create_thread(boost::bind(start, maptool, moleSetPtr, indexesSetPtr, centerSetPtr, thread, i));
            group.create_thread(boost::bind(start, maptool, moleSet, indexesSet, centerSetPtr, thread, i));
        }
        group.join_all();
        for(int i = 0; i < thread; ++ i) {
            std::cout << i << " " << threadCenterSet[i].size() << std::endl;

            for(auto it = threadCenterSet[i].begin(); it != threadCenterSet[i].end(); ++ it) {
                centerSet.insert(*it);
            }
        }
        for(int i = 0; i < centerSet.size(); ++ i) {
            if(centerSet[i].size() != 1) {
                LOG4CXX_DEBUG(logger, boost::format("cluster %s is a mixCluster, divided into %d clusters.") % i % centerSet[i].size());
            }
            for(auto it = centerSet[i].begin(); it != centerSet[i].end(); ++ it) {
                resOutstream << i << '-' << it - centerSet[i].begin() <<  '\n';
                Center voteCenter = *it;
                for(int i = 0; i < voteCenter.size(); ++ i) {
                    resOutstream << voteCenter[i] << " ";
                }
                resOutstream << '\n';
            }
        }
        return true;
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % clusterFile);
        return false;
    }
}
