#include "contig_builder.h"

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

typedef std::vector<std::set<Vertex>> Components;
typedef std::map<int, std::set<int> > SiteGraph;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.contig"));
void ContigBuilder::start(const std::vector<Mole>* moleSetPtr, std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId) const {
    const std::vector<Mole>& moleSet = *moleSetPtr;
    std::vector<Alignment>& alignments = *alignmentsPtr;
    for(int i = threadId; i < moleSet.size(); i += threads) {
        for(int j = i + 1; j < moleSet.size(); ++ j) {
            Alignment align = _maptool.localDPscore(moleSet[i], moleSet[j]);
            //a hard threshold
            if(align.getScore() < minScore || align.getScore() / align.size() < (minScore / 20)) {
                continue;
            }
            alignments.push_back(align);
        }
    }
}
void ContigBuilder::alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, int threads, double minScore) const {
    const std::vector<Mole>* moleSetPtr = &moleSet;
    std::vector<std::vector<Alignment>> threadAlignments(threads, std::vector<Alignment>());
    boost::thread_group group;
    LOG4CXX_INFO(logger, boost::format("Alignment using %d threads.") % threads);
    for(int i = 0; i < threads; ++ i) {
        std::vector<Alignment>* alignmentsPtr = &threadAlignments[i];
        group.create_thread(boost::bind(&ContigBuilder::start, this, moleSetPtr, alignmentsPtr, minScore, threads, i));
    }
    group.join_all();
    for(int i = 0; i < threads; ++ i) {
        for(auto &al : threadAlignments[i]) {
            al.trimHead();
            al.trimTail();
            alignments.push_back(al);
            al.print(std::cout, true);
        }
    }
}


void ContigBuilder::constructGraph(const std::vector<Alignment>& alignments, Graph& graph1, Graph& graph2, edgeSet& edges) const {
    //filter alignment using LSAT
    std::map<std::string, std::vector<int>> moleAlignment;
    

    //connect the sites
    for(Alignment al : alignments) {
        std::string m1 = al.mole1Id;
        std::string m2 = al.mole2Id;
        int s1 = al.mole1Start;
        int s2 = al.mole2Start;
        Vertex vertex1(m1, s1);
        Vertex vertex2(m2, s2);
        std::string edge1 = vertex1.to_string() + "-" + vertex2.to_string();
        std::string edge2 = vertex2.to_string() + "-" + vertex1.to_string();

        std::cout << edge1 << "-" << 3.6 << std::endl;
        graph1[vertex1].insert(vertex2);
        graph1[vertex2].insert(vertex1);
        if(s1 == 0 || s2 == 0) {
            graph2[vertex1].insert(vertex2);
            graph2[vertex2].insert(vertex1);
        }
        for(int i = 0; i < al.alignedMole1.size(); ++ i) {
            Fragment f1 = al.alignedMole1[i], f2 = al.alignedMole2[i];
            double leftScore = _maptool.validScore(f1, f2);
            double rightScore = 0;
            if(i !=  al.alignedMole1.size() - 1) {
                Fragment f3 = al.alignedMole1[i + 1], f4 = al.alignedMole2[i + 1];
                rightScore = _maptool.validScore(f3, f4);
            }
            double edgeScore = std::min(leftScore, rightScore);
            s1 += al.alignedMole1[i].size();
            s2 += al.alignedMole2[i].size();
            if(edgeScore < 1.5) continue;
            vertex1.startSite = s1;
            vertex2.startSite = s2;
            graph1[vertex1].insert(vertex2);
            graph1[vertex2].insert(vertex1);
            //the score of edge is used as a threshold in construct graph.
            std::string edge1 = vertex1.to_string() + "-" + vertex2.to_string();
            std::string edge2 = vertex2.to_string() + "-" + vertex1.to_string();

            std::cout << edge1 << "-" << edgeScore << std::endl;

            edges[edge1] = al.score;
            edges[edge2] = al.score;
        }
    }
}
void bfsSearch(Graph& graph, edgeSet& edges, int minCluster, Components& comps) {
    //local soft adjustable threshold
    Graph refinedGraph;
    for(auto vertex : graph) {
        Vertex v = vertex.first;
        for(Vertex u : vertex.second) {
            std::vector<Vertex> set_intersection, set_union;
            std::set_intersection(graph[v].begin(), graph[v].end(), graph[u].begin(), graph[u].end(), std::back_inserter(set_intersection));
            std::set_union(graph[v].begin(), graph[v].end(), graph[u].begin(), graph[u].end(), std::back_inserter(set_union));
            //remove edges
            if(set_intersection.size() == 2 && (graph[v].size() > 4 || graph[u].size() > 4)) {
                continue;
            }

            if(set_intersection.size() && set_union.size() / set_intersection.size() < 4) {
                refinedGraph[v].insert(u);
            }
        }

    }
    std::cout << "refine finished" << std::endl;
    //bfs visit to find all connected components.
    std::map<Vertex, int> colors;
    for(auto it : refinedGraph) {
        colors[it.first] = 0;
    }
    std::queue<Vertex> q;
    std::set<Vertex> comp;
    int count = 0;
    for(auto it : refinedGraph) {
        if(colors[it.first] != 0) continue;
        comp.clear();
        q.push(it.first);
        while(!q.empty()) {
            Vertex tmp = q.front();
            colors[tmp] = 2;
            comp.insert(tmp);
            auto tmpAdj = refinedGraph.find(tmp);
            q.pop();
            for(auto iv : tmpAdj->second) {
                if(colors[iv] == 0) {
                    q.push(iv);
                    colors[iv] = 1;
                }
            }
        }
        /*
         * edge is not used
         *
        if(comp.size() == 2) {
            std::string edge = comp[0].moleId + std::to_string(comp[0].startSite) + comp[1].moleId + std::to_string(comp[1].startSite);
            if(edges[edge] < 30)
                continue;
        }
        */
        count ++;
        LOG4CXX_DEBUG(logger, boost::format("#%d cluster, cluster size is %d.") % count % comp.size());
        //filter small comp
        if(comp.size() <= 2) continue;
        comps.push_back(comp);
    }
}
void divide2(std::set<Vertex>& component, const Graph& graph, Components& dividedComponents) {
    std::map<std::string, int> vertexInMole;
    bool isMixed = false;
    for(Vertex vertex : component) {
        if(vertexInMole[vertex.moleId] != 0) isMixed = true;
        auto it = graph.find(vertex);
        vertexInMole[vertex.moleId] += it->second.size();
    }
    // not mix
    if(! isMixed) {
        dividedComponents.push_back(component);
        return;
    }
    // find a template
    std::string templateMole;
    int maxSupp = -1;
    for(auto it : vertexInMole) {
        if(it.second > maxSupp) {
            maxSupp = it.second;
            templateMole = it.first;
        }
    }
    // divide the component
    for(auto v : component) {
        if(v.moleId == templateMole) {
            dividedComponents.push_back(std::set<Vertex>({v}));
            for(auto u : component) {
                auto it = graph.find(v);
                if(it->second.find(u) != it->second.end()) {
                    dividedComponents.back().insert(u);
                }
            }
        }
    }
    // debug information
    if(dividedComponents.size() > 1) {
        std::cout << "divide: " << dividedComponents.size() << "\n";
        for(Vertex v : component) {
            std::cout << v.to_string() << ",";
        }
        std::cout << "\n";
        for(int i = 0; i < dividedComponents.size(); ++ i) {
            for(Vertex v : dividedComponents[i]) {
                std::cout << v.to_string() << ",";
            }
            std::cout << "\n";
        }
    }
}
void divide1(std::set<Vertex>& component, const Graph& graph, Components& dividedComponents) {
    std::map<std::string, std::vector<Vertex> > vertexInMole;
    for(Vertex vertex : component) {
        if(vertexInMole.find(vertex.moleId) == vertexInMole.end()) {
            vertexInMole[vertex.moleId] = std::vector<Vertex>();
        }
        vertexInMole[vertex.moleId].push_back(vertex);
    }
    std::map<int, int> vertexCount;
    for(auto it : vertexInMole) {
        if(vertexCount.find(it.second.size()) == vertexCount.end()) {
            
            vertexCount[it.second.size()] = 0;
        }
        vertexCount[it.second.size()] ++;
    }
    int maxCount = -1;
    int maxCountVertex = -1;

    for(auto it : vertexCount) {
        if(maxCount < it.second) {
            maxCount = it.second;
            maxCountVertex = it.first;
        }
    }
    for(auto it : vertexInMole) {
        if(it.second.size() == maxCountVertex) {
            for(int i = 0; i < maxCountVertex; ++ i) {
                if(dividedComponents.size() <= i) {
                    dividedComponents.push_back(std::set<Vertex>());
                }
                dividedComponents[i].insert(it.second[i]);
            }
        }
    }

    // debug information
    if(dividedComponents.size() > 1) {
        std::cout << "divide: " << dividedComponents.size() << "\n";
        for(Vertex v : component) {
            std::cout << v.to_string() << ",";
        }
        std::cout << "\n";
        for(int i = 0; i < dividedComponents.size(); ++ i) {
            for(Vertex v : dividedComponents[i]) {
                std::cout << v.to_string() << ",";
            }
            std::cout << "\n";
        }
    }
    
}
void EMCalling() {
    
}
void topoSort(SiteGraph& spreGraph, SiteGraph& spostGraph, std::vector<int>& spath) {
    // get start
    std::cout << "contig topo" << std::endl;
    for(auto site : spreGraph) {
        for(auto nsite : site.second) {
            std::cout << site.first << "->" << nsite << ";" << std::endl;
        }
        if(spostGraph.find(site.first) == spostGraph.end()) {
            spath.push_back(site.first);
            for(auto nextsite : spostGraph) {
                if(nextsite.second.find(site.first) != nextsite.second.end()) {
                    spostGraph[nextsite.first].erase(site.first);
                }
            }
        }
    }
    while(spostGraph.size() != 0) {
        int tempSite = -1;
        for(auto site : spostGraph) {
            if(site.second.size() == 0) {
                tempSite = site.first;
                break;
            }
        }
        if(tempSite == -1 && spostGraph.size() != 0) {
            std::cout << spath.size() << " " << spostGraph.size() << std::endl;
            LOG4CXX_INFO(logger, boost::format("there is a circle in siteGraph."));
            for(auto site : spostGraph) {
                if(site.second.size() != 0) {
                    std::cout << site.first << " [style=\"filled\", fillcolor=\"red\"];" << std::endl;
                }
            }
            std::cout << std::endl;
            break;
        }
        spath.push_back(tempSite);
        spostGraph.erase(tempSite);
        for(int site : spreGraph[tempSite]) {
            spostGraph[site].erase(tempSite);
        }
    }
}

std::vector<std::pair<int, int> > nextSites(const Components& comps, const std::set<Vertex>& curcomp, const std::set<int>& moleComps, std::vector<Mole>& moleSet) {
    //TODO: slove wrong connection
    int nextcomp = 1;
    std::vector<std::pair<int, int>> res;
    std::vector<int> distances;
    int minDistance = INT_MAX;
    std::map<std::string, int> siteOnmole;
    for(int i : moleComps) {
        auto comp = comps[i];
        distances.clear();
        bool label = false;
        for(Vertex ns : comp) {
            for(Vertex ps : curcomp) {
                if(ns.moleId == ps.moleId && ns.startSite > ps.startSite) {
                    distances.push_back(ns.startSite - ps.startSite);
                    if(siteOnmole[ps.moleId] == 0) siteOnmole[ps.moleId] = ns.startSite;
                    else siteOnmole[ps.moleId] = std::min(siteOnmole[ps.moleId], ns.startSite);
                }
                // detect circle
                if(ns.moleId == ps.moleId && ns.startSite <= ps.startSite) {
                    label = true;
                    break;
                }
            }
        }
        if(label || distances.size() == 0) continue;
        res.push_back(std::make_pair(i, distances.size()));
    }
    auto f = [&](std::pair<int, int>& a){
        for(Vertex v : comps[a.first]) {
            if(siteOnmole.find(v.moleId) != siteOnmole.end() && siteOnmole[v.moleId] == v.startSite) return false;
        }
        return true;
    };
    res.erase(std::remove_if(res.begin(), res.end(), f), res.end());
    //debug information
    std::unordered_map<std::string, size_t> idIndex;
    for(int i = 0; i < moleSet.size(); ++ i) {
        idIndex[moleSet[i].getID()] = i;
    }
    std::cout << "next comps:" << std::endl;
    for(Vertex ns : curcomp) {
        std::cout << ns.moleId << "|" << ns.startSite << ",";
    }
    std::cout << std::endl;
    if(nextcomp != -1) {
        for(auto it : res) {
            nextcomp = it.first;
            for(Vertex ns : comps[nextcomp]) {
                std::cout << ns.moleId << "|" << ns.startSite << ",";
            }
            std::cout << std::endl;
            std::vector<int> distance;
            std::cout << nextcomp << std::endl;
            for(Vertex ns : comps[nextcomp]) {
                for(Vertex ps : curcomp) {
                    if(ns.moleId == ps.moleId) {
                        for(int ss = ps.startSite; ss < ns.startSite; ++ ss) {
                            std::cout << ps.to_string() << " " << ns.to_string() << " " << moleSet[idIndex[ps.moleId]].getInterval(ss) << std::endl;
                        }
                        distance.push_back(moleSet[idIndex[ps.moleId]].getInterval(ps.startSite));
                    }
                }
            }
            for(int i = 0; i < distance.size(); ++ i) {
                std::cout << distance[i] << " ";
            }
            std::cout << accumulate(distance.begin(), distance.end(), 0) / distance.size() << std::endl;
        }
    }
    return res;
}

void component(const SiteGraph& sgraph, std::vector<std::set<int>>& parts) {
    // bfs search
    std::map<int, int> colors;
    for(auto it : sgraph) {
        colors[it.first] = 0;
    }
    std::queue<int> q;
    std::set<int> comp;
    int count = 0;
    for(auto it : sgraph) {
        if(colors[it.first] != 0) continue;
        comp.clear();
        q.push(it.first);
        while(!q.empty()) {
            int tmp = q.front();
            colors[tmp] = 2;
            comp.insert(tmp);
            auto tmpAdj = sgraph.find(tmp);
            q.pop();
            for(auto iv : tmpAdj->second) {
                if(colors[iv] == 0) {
                    q.push(iv);
                    colors[iv] = 1;
                }
            }
        }
        count ++;
        LOG4CXX_DEBUG(logger, boost::format("#%d cluster, cluster size is %d.") % count % comp.size());
        parts.push_back(comp);
    }
}

bool connect(const Components& comps, std::vector<Mole>& moleSet, std::vector<std::vector<int>>& contigs) {
    //construct contig graph
    SiteGraph sGraph, spreGraph, spostGraph;
    std::map<std::string, std::vector<int>> moleInComps;
    for(int i = 0; i < comps.size(); ++ i) {
        for(Vertex v : comps[i]) {
            moleInComps[v.moleId].push_back(i);
        }
    }
    std::vector<int> contig;
    for(int i = 0; i < comps.size(); ++ i) {
        auto curcomp = comps[i];
        std::set<int> moleComps;
        for(Vertex v : curcomp) {
            moleComps.insert(moleInComps[v.moleId].begin(), moleInComps[v.moleId].end());
        }
        auto nextComps = nextSites(comps, curcomp, moleComps, moleSet);
        // add a edge between comps
        if(nextComps.size() != 0) {
            // maxSupport is the support number of most confident nextSite
            int maxSupport = std::max_element(nextComps.begin(), nextComps.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b){return a.second < b.second;})->second;
            /*
            an edge filter
            
                   [s1]
                   /  \
                 5/    \1
                 /      \
              [n1]      [n2]
            
            [s1] -- [n2] will be removed as its maxSupport / support number = 5
            */
            for(auto nextcomp : nextComps) {
                if(comps[i].size() / comps[nextcomp.first].size() > 9 || comps[nextcomp.first].size() / comps[i].size() > 9 || nextcomp.first == i || (maxSupport >= 5 && nextcomp.second == 1) || maxSupport / comps[nextcomp.first].size() > 5) continue;
                // debug information for dot graphviz
                std::cout << i << "--" << nextcomp.first << " [label=" << nextcomp.second << "];" << std::endl;
                // constuct graph for longest path
                spreGraph[i].insert(nextcomp.first);
                spostGraph[nextcomp.first].insert(i);
                sGraph[i].insert(nextcomp.first);
                sGraph[nextcomp.first].insert(i);
            }
        }
    }
    /*
     * all components mix in one dp process is not a good idea
     * if one component has a circle, there will be no contigs 
     */
    std::vector<std::set<int>> parts;
    // bfs for contig
    component(sGraph, parts);
    LOG4CXX_DEBUG(logger, boost::format("site graph has %d components.") % parts.size());

    // build subgraph, get topo order
    for(auto part : parts) {
        contig.clear();
        std::vector<int> spath;
        SiteGraph partPreGraph, partPostGraph;
        for(auto node : spreGraph) {
            if(part.find(node.first) != part.end()) {
                partPreGraph.insert(node);
            }
        }
        for(auto node : spostGraph) {
            if(part.find(node.first) != part.end()) {
                partPostGraph.insert(node);
            }
        }
        topoSort(partPreGraph, partPostGraph, spath);

        // exists a circle, drop out this contig
        if(spath.size() != partPreGraph.size()) {
            LOG4CXX_INFO(logger, boost::format("there is a circle in siteGraph."));
            //continue;
        }
        // debug information
        std::cout << "spath: " << std::endl;
        for(int i = 0; i < spath.size(); ++ i) {
            std::cout << spath[i] << " ";
        }
        std::cout << std::endl;

        // get longest path
        std::map<int, int> dpLength;
        std::map<int, int> dpPath;
        for(int i = 0; i < spath.size(); ++ i) {
             dpPath[spath[i]] = -1;
        }
        for(int i = 0; i < spath.size(); ++ i) {
            for(int postSite : spreGraph[spath[i]]) {
                int weight = 0;
                for(Vertex ns : comps[postSite]) {
                    for(Vertex ps : comps[spath[i]]) {
                        if(ns.moleId == ps.moleId && ns.startSite - ps.startSite == 1) {
                            weight ++;
                        }
                    }
                }
                if(dpLength[spath[i]] + weight > dpLength[postSite]) {
                    dpLength[postSite] = dpLength[spath[i]] + weight;
                    dpPath[postSite] = spath[i];
                }
            }
        }
        // debug information
        std::cout << "dpLength: " << std::endl;
        for(auto it : dpLength) {
            std::cout << it.first << "-" << it.second << " ";
        }
        std::cout << std::endl;
        std::cout << "dpPath: " << std::endl;
        for(auto it : dpPath) {
            std::cout << it.first << "-" << it.second << " ";
        }
        std::cout << std::endl;

        // trace back
        int pathEnd = -1;
        int maxLength = -1;
        std::vector<int> contigPath;
        // end of longest path
        for(int i = 0; i < spath.size(); ++ i) {
            if(dpLength[spath[i]] > maxLength && part.find(spath[i]) != part.end()) {
                maxLength = dpLength[spath[i]];
                pathEnd = spath[i];
            }
        }
        // get path
        while(pathEnd != -1) {
            contigPath.push_back(pathEnd);
            pathEnd = dpPath[pathEnd];
        }
        // debug information
        std::reverse(contigPath.begin(), contigPath.end());
        for(int i = 0; i < contigPath.size(); ++ i) {
            std::cout << contigPath[i] << " ";
        }
        std::cout << std::endl;
        contigs.push_back(contigPath);
    }
}

int calcDistance(std::vector<int> distance) {
    /*
     * claculate the interval length
     */
    int tempDis = 0;
    if(distance.size() == 1)
        tempDis = distance[0];
    if(distance.size() == 2)
        // just return average of them
        tempDis = std::accumulate(distance.begin(), distance.end(), 0.0) / distance.size();
    else if(distance.size() == 3) {
        // return average of two of them which are closer
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
        // cluster these numbers by a threshold(default: 500)
        // return the average of the biggest cluster
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
        // max cluster size is 1
        if(maxCluster == 1) {
            maxIndex = 0, maxCluster = distance.size();
        }
        // cluster average
        tempDis = std::accumulate(distance.begin() + maxIndex, distance.begin() + maxIndex + maxCluster, 0.0) / (maxCluster);
    }
    std::cout << tempDis << std::endl;
    return tempDis;
}

std::vector<int> intervalCalling(const std::set<Vertex>& startComp, const std::set<Vertex>& endComp, const Graph& graph, const std::vector<Mole>& moles) {
    /*
     * call the interval between start and end comp
     */
    std::map<std::string, std::pair<int, int>> moleInterval;
    for(auto v : startComp) {
        moleInterval[v.moleId] = std::make_pair(v.startSite, -1);
    }
    for(auto v : endComp) {
        if(moleInterval.find(v.moleId) != moleInterval.end()) 
            moleInterval[v.moleId].second = v.startSite;
    }

    // get a center mole
    std::map<std::string, int> moleConnect;
    int maxConnect = -1;
    std::string centerMole;
    for(auto interval : moleInterval) {
        Vertex v1(interval.first, interval.second.first), v2(interval.first, interval.second.second);
        auto it1 = graph.find(v1), it2 = graph.find(v2);
        if(it1 != graph.end() && it2 != graph.end()) {
            moleConnect[interval.first] += it1->second.size();
            moleConnect[interval.first] += it2->second.size();
        }
        // record max connection
        if(moleConnect[interval.first] > maxConnect) {
            maxConnect = moleConnect[interval.first];
            centerMole = interval.first;
        }
    }
    LOG4CXX_DEBUG(logger, boost::format("Center mole is %s, and its connection is %d.") % centerMole % maxConnect);
    // moles ID to index
    std::unordered_map<std::string, size_t> idIndex;
    for(size_t i = 0; i < moles.size(); ++ i) {
        idIndex[moles[i].getID()] = i;
    }

    // call the interval site by site
    // TODO: handle deletion in center mole
    std::vector<int> contig;
    auto preComp = startComp;
    std::cout << "centerMole" << centerMole << std::endl;
    for(int i = moleInterval[centerMole].first; i < moleInterval[centerMole].second; ++ i) {
        auto nextSite = Vertex(centerMole, i + 1);
        std::cout << nextSite.to_string() << std::endl;
        auto it = graph.find(nextSite);
        if(it == graph.end()) continue;
        auto nextComp = it->second;
        nextComp.insert(nextSite);
        /*
        if(nextComp.size() <= startComp.size() / 5) {
            LOG4CXX_DEBUG(logger, boost::format("%s is a insertion.") % nextSite.to_string());
            continue;
        }
        */
        // handle insertion
        std::vector<int> moleSites;
        std::vector<int> distances;
        for(auto v : preComp) {
            for(auto u : nextComp) {
                if(v.moleId == u.moleId) {
                    moleSites.push_back(u.startSite - v.startSite);
                    int intervalDis = 0;
                    for(int i = v.startSite; i < u.startSite; ++ i) {
                        intervalDis += moles[idIndex[v.moleId]].getInterval(i);
                    }
                    distances.push_back(intervalDis);
                }
            }
        }
        for(auto d : distances) {
            std::cout << d << " ";
        }
        contig.push_back(calcDistance(distances));
        std::cout << contig.back() << std::endl;
        preComp = nextComp;
    }
    // debug information
    for(auto pv : startComp) {
        std::cout << pv.to_string() << " ";
    }
    std::cout << std::endl;
    for(auto nv : endComp) {
        std::cout << nv.to_string() << " ";
    }
    std::cout << std::endl;

    for(int i = 0; i < contig.size(); ++ i) {
        std::cout << contig[i] << " ";
    }
    std::cout << "contig" << std::endl;
    return contig;
}

void outputBNX(const std::vector<std::vector<int>>& contigs, const std::string& output) {
    std::ofstream os(output.c_str());
    for(size_t i = 0; i < contigs.size(); ++ i) {
        int startPos = 0;
        int totalLength = std::accumulate(contigs[i].begin(), contigs[i].end(), 0);
        os << "0" << "\t" << i << "\t" << totalLength << "\t" << "15" << "\t" << "15" << "\t" << contigs[i].size() + 1 << "\t" << i << "1" << "\t" << "-1" << "\t" << CHIP_ID << "2" << "\t" << "1" << "\t" << "1" << "\n";
        os << "1" << "\t" << startPos << "\t";
        for(size_t j = 0; j < contigs[i].size(); ++ j) {
            startPos += contigs[i][j];
            os << startPos << "\t";
        }
        os << startPos << "\n";
        os << "QX11" << "\t" << "15";
        for(size_t j = 0; j < contigs[i].size(); ++ j) {
            os << "\t" << "15";
        }
        os << "\n";
        os << "QX12" << "\t" << "0.03";
        for(size_t j = 0; j < contigs[i].size(); ++ j) {
            os << "\t" << "0.03";
        }
        os << "\n";
    }
    os.close();
}
void outputCMP(const std::vector<std::vector<int>>& contigs, const std::string& output) {
    std::ofstream os(output.c_str());
    os << CMP_HEADER;
    for(size_t i = 0; i < contigs.size(); ++ i) {
        int startPos = 0;
        int totalLength = std::accumulate(contigs[i].begin(), contigs[i].end(), 0);
        os << i + 1 << "\t" << totalLength << "\t" << contigs[i].size() + 1 << "\t" << 1 << "\t" << "1" << "\t" << startPos << "\t" << "1.0" << "\t" << "1" << "\t" << "1" << "\n";
        for(size_t j = 0; j < contigs[i].size(); ++ j) {
            startPos += contigs[i][j];
            os << i + 1 << "\t" << totalLength << "\t" << contigs[i].size() + 1 << "\t" << j + 2 << "\t" << "1" << "\t" << startPos << "\t" << "1.0" << "\t" << "1" << "\t" << "1" << "\n";
        }
        os << i + 1 << "\t" << totalLength << "\t" << contigs[i].size() + 1 << "\t" << contigs[i].size() + 2 << "\t" << "0" << "\t" << startPos << "\t" << "1.0" << "\t" << "1" << "\t" << "1" << "\n";
    }
    os.close();
}
void outputContig(const std::vector<std::vector<int>>& contigs, const std::string& output) {
    std::string bnxFile = output + ".bnx";
    outputBNX(contigs, bnxFile);
    std::string cmpFile = output + ".cmap";
    outputCMP(contigs, cmpFile);
}

bool ContigBuilder::build(const std::string& input, const std::string& output, const std::string& alignmentFile, double minScore, int threads) const {

    // read mole
    LOG4CXX_WARN(logger, boost::format("alignment score threshold is %s") % minScore);
    std::vector<Mole> moleSet;
    int reverseLabel = 1;
    if(boost::filesystem::exists(input)) {
        std::ifstream moleInstream(input.c_str());
        int count = -1;
        /* 
        std::string line;
        while(std::getline(moleInstream, line)) {
            count ++;
            std::string lineName = std::to_string(count);
            std::getline(moleInstream, line);
            boost::trim(line);
            std::vector<std::string> lineSplit;
            boost::split(lineSplit, line, boost::is_any_of(" \t"), boost::token_compress_on);
            std::vector<int> tempMole;
            BOOST_FOREACH(std::string interval, lineSplit) {
                tempMole.push_back(boost::lexical_cast< int > (interval));
            }
            Mole m(lineName, tempMole);
            moleSet.push_back(m);
        }
        */
        MoleReader mReader(moleInstream);
        Mole m;
        while(mReader.read(m)) {
            if(m.size() < 10) continue;
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
            LOG4CXX_WARN(logger, boost::format("%s moles have been inited.") % moleNumber);
        } 
    }

    // read alignment
    std::vector<Alignment> alignments;
    if(boost::filesystem::exists(alignmentFile)) {
        std::ifstream alignInstream(alignmentFile.c_str());
        AlignmentReader aReader(alignInstream);
        Alignment al;
        while(aReader.read(al)) {
            if(al.getScore() < minScore || al.getScore() / al.size() < minScore / 20) continue;
            alignments.push_back(al);
        }
        int alignmentNumber = alignments.size();
        if (alignmentNumber == 0) {
            LOG4CXX_WARN(logger, "no alignment is in alignmentFile");
            return false;
        } else {
            LOG4CXX_INFO(logger, boost::format("%s alignments have been inited.") % alignmentNumber);
        } 
    } else {
        // no alignment file
        alignment(moleSet, alignments, threads, minScore);
    }

    // construct graph, g1 for connect, g2 for cluster
    Graph g1, g2;
    edgeSet edges;
    std::cout << "construct graph" << std::endl;
    constructGraph(alignments, g1, g2, edges);
    std::cout << "construct graph" << std::endl;
    
    // bfs to cluster
    Components comps;
    bfsSearch(g1, edges, 2, comps);
    std::cout << "bfs result:" << std::endl;
    for(auto cp : comps) {
        for(Vertex ns : cp) {
            std::cout << ns.moleId << "|" << ns.startSite << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "bfs result" << std::endl;

    // divide
    Components dividedComponents;
    Components tempComp;
    for(auto comp : comps) {
        tempComp.clear();
        divide2(comp, g1, tempComp);
        LOG4CXX_DEBUG(logger, boost::format("divide into %d components.") % tempComp.size());
        for(auto cp : tempComp) {
            dividedComponents.push_back(cp);
        }
    }
    // debug information
    std::cout << "comps result:" << std::endl;
    for(auto cp : dividedComponents) {
        for(Vertex ns : cp) {
            std::cout << ns.moleId << "|" << ns.startSite << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "comps result" << std::endl;

    // contiging
    std::vector<std::vector<int>> siteContigs;
    connect(dividedComponents, moleSet, siteContigs);
    LOG4CXX_DEBUG(logger, boost::format("there are %d contigs.") % siteContigs.size());
    
    // interval call
    std::vector<std::vector<int>> contigs;
    for(int i = 0; i < siteContigs.size(); ++ i) {
        std::vector<int> contig;
        for(int j = 0; j < siteContigs[i].size() - 1; ++ j) {
            auto intervals = intervalCalling(dividedComponents[siteContigs[i][j]], dividedComponents[siteContigs[i][j + 1]], g1, moleSet);
            for(int k : intervals) contig.push_back(k);
        }
        if(contig.size() >= 20)
            contigs.push_back(contig);
    }
    outputContig(contigs, output);
    return true;
}
