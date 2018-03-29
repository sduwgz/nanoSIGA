#include "contig_builder.h"

#include <iostream>
#include <fstream>
#include <queue>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include <log4cxx/logger.h>

typedef std::map<Vertex, std::set<Vertex>> Graph;
typedef std::vector<std::set<Vertex>> Components;
typedef std::map<std::string, double> edgeSet;
typedef std::map<int, std::set<int> > SiteGraph;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.contig"));
void ContigBuilder::start(const std::vector<Mole>* moleSetPtr, std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId) const {
    const std::vector<Mole>& moleSet = *moleSetPtr;
    std::vector<Alignment>& alignments = *alignmentsPtr;
    for(int i = threadId; i < moleSet.size(); i += threads) {
        for(int j = i + 1; j < moleSet.size(); ++ j) {
            Alignment align = _maptool.localDPscore(moleSet[i], moleSet[j]);
            //a hard threshold
            if(align.score < minScore || align.score / align.alignedMole1.size() < (minScore / 10)) {
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
        }
    }
}


void constructGraph(const std::vector<Alignment>& alignments, Graph& graph, edgeSet& edges) {
    //connect the sites which are corresponded in alignments.
    for(Alignment al : alignments) {
        std::string m1 = al.mole1Id;
        std::string m2 = al.mole2Id;
        int s1 = al.mole1Start;
        int s2 = al.mole2Start;
        Vertex vertex1(m1, s1);
        Vertex vertex2(m2, s2);
        graph[vertex1].insert(vertex2);
        graph[vertex2].insert(vertex1);
        for(int i = 0; i < al.alignedMole1.size(); ++ i) {
            s1 += al.alignedMole1[i].size();
            s2 += al.alignedMole2[i].size();
            vertex1.startSite = s1;
            vertex2.startSite = s2;
            graph[vertex1].insert(vertex2);
            graph[vertex2].insert(vertex1);
            //the score of edge is used as a threshold in construct graph.
            //std::string edge = vertex1.moleId + std::to_string(vertex1.startSite) + vertex2.moleId + std::to_string(vertex2.startSite);
            //edges[edge] = al.score;
        }
    }
}
void bfsSearch(const Graph& graph, edgeSet& edges, int minCluster, Components& comps) {
    //bfs visit to find all connected components.
    std::map<Vertex, int> colors;
    for(auto it : graph) {
        colors[it.first] = 0;
    }
    std::queue<Vertex> q;
    std::set<Vertex> comp;
    int count = 0;
    for(auto it : graph) {
        if(colors[it.first] != 0) continue;
        comp.clear();
        q.push(it.first);
        while(!q.empty()) {
            Vertex tmp = q.front();
            colors[tmp] = 2;
            comp.insert(tmp);
            auto tmpAdj = graph.find(tmp);
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
        //small comp is filtered
        if(comp.size() < 4) continue;
        comps.push_back(comp);
    }
}
void divide(std::set<Vertex>& component, const Graph& graph, Components& dividedComponents) {
    //there are some mixed components, divide them into several components.
    std::map<std::string, std::vector<Vertex> > vertexInMole;
    std::cout << component.size() << std::endl;
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
    //debug information
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

void topoSort(SiteGraph& spreGraph, SiteGraph& spostGraph, std::vector<int>& spath) {
    //Toposort maybe drop in infinite loops
    for(auto site : spreGraph) {
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
            LOG4CXX_INFO(logger, boost::format("there is a circle in siteGraph."));
            break;
        }
        std::cout << tempSite << std::endl;
        spath.push_back(tempSite);
        spostGraph.erase(tempSite);
        for(int site : spreGraph[tempSite]) {
            spostGraph[site].erase(tempSite);
        }
    }
}

std::vector<int> nextSites(const Components& comps, const std::set<Vertex>& curcomp, std::vector<Mole>& moleSet) {
    //TODO: how to slove wrong connection
    int nextcomp = -1;
    std::vector<int> res;
    std::vector<int> distances;
    int maxOne = 0;
    int minDistance = INT_MAX;
    for(int i = 0; i < comps.size(); ++ i) {
        auto comp = comps[i];
        distances.clear();
        for(Vertex ns : comp) {
            for(Vertex ps : curcomp) {
                if(ns.moleId == ps.moleId && ns.startSite > ps.startSite) {
                    distances.push_back(ns.startSite - ps.startSite);
                }
                if(distances.size() == 0) continue;

                int d = accumulate(distances.begin(), distances.end(), 0) / distances.size();
                //
                if(d == minDistance) {
                    res.push_back(i);
                } else if(d < minDistance) {
                    minDistance = d;
                    nextcomp = i;
                    res.clear();
                    res.push_back(i);
                }
            }
        }
    }
    //debug information
    std::cout << "next comps:" << std::endl;
    for(Vertex ns : curcomp) {
        std::cout << ns.moleId << "|" << ns.startSite << ",";
    }
    std::cout << std::endl;
    if(nextcomp != -1) {
        for(Vertex ns : comps[nextcomp]) {
            std::cout << ns.moleId << "|" << ns.startSite << ",";
        }
        std::cout << std::endl;
        std::vector<int> distance;

        std::cout << comps.size() << std::endl;
        std::cout << nextcomp << std::endl;
        for(Vertex ns : comps[nextcomp]) {
            for(Vertex ps : curcomp) {
                if(ns.moleId == ps.moleId) {
                    std::cout << ps.to_string() << " " << ns.to_string() << std::endl;
                    distance.push_back(moleSet[boost::lexical_cast<int>(ps.moleId)].getInterval(ps.startSite));
                }
            }
        }
        for(int i = 0; i < distance.size(); ++ i) {
            std::cout << distance[i] << " ";
        }
        std::cout << accumulate(distance.begin(), distance.end(), 0) / distance.size() << std::endl;
    }
    return res;
}

void component(const SiteGraph& sgraph, std::vector<std::set<int>>& parts) {
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
    SiteGraph sGraph, spreGraph, spostGraph;
    std::vector<int> contig;
    for(int i = 0; i < comps.size(); ++ i) {
        auto curcomp = comps[i];
        std::vector<int> nextComps = nextSites(comps, curcomp, moleSet);
        if(nextComps.size() != 0) {
            for(int nextcomp : nextComps) {
                if(nextcomp == i) continue;
                spreGraph[i].insert(nextcomp);
                spostGraph[nextcomp].insert(i);
                sGraph[i].insert(nextcomp);
                sGraph[nextcomp].insert(i);
            }
        }
    }
    //all components mix in one dp process is not a good idea.
    //if one component has a circle, there will be no contig. 
    std::vector<std::set<int>> parts;
    component(sGraph, parts);
    LOG4CXX_DEBUG(logger, boost::format("site graph has %d components.") % parts.size());
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
        if(spath.size() != partPreGraph.size()) {
            LOG4CXX_INFO(logger, boost::format("there is a circle in siteGraph."));
            continue;
        }

        std::cout << "spath: " << std::endl;
        for(int i = 0; i < spath.size(); ++ i) {
            std::cout << spath[i] << " ";
        }
        std::cout << std::endl;

        //dynamic programming for longest path
        std::map<int, int> dpLength;
        std::map<int, int> dpPath;
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
                std::cout << spath[i] << "-" << postSite << "-" << weight << std::endl;
                if(dpLength[spath[i]] + weight > dpLength[postSite]) {
                    dpLength[postSite] = dpLength[spath[i]] + weight;
                    dpPath[postSite] = spath[i];
                }
            }
        }
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
        //trace back
        int pathEnd = -1;
        int maxLength = -1;
        std::vector<int> contigPath;
        for(int i = 0; i < spath.size(); ++ i) {
            if(dpLength[spath[i]] > maxLength && part.find(spath[i]) != part.end()) {
                maxLength = dpLength[spath[i]];
                pathEnd = spath[i];
            }
        }
        std::cout << "maxLength: " << maxLength << "pathEnd: " << pathEnd << std::endl;
        while(pathEnd != 0) {
            contigPath.push_back(pathEnd);
            pathEnd = dpPath[pathEnd];
        }
        for(int i = 0; i < contigPath.size(); ++ i) {
            std::cout << contigPath[i] << " ";
        }
        std::cout << std::endl;
        for(int i = contigPath.size() - 2; i >= 0; -- i) {
            std::cout << "site: " << contigPath[i] << std::endl;
            std::vector<int> distance;
            for(Vertex ns : comps[contigPath[i]]) {
                for(Vertex ps : comps[contigPath[i + 1]]) {
                    if(ns.moleId == ps.moleId && ns.startSite - ps.startSite == 1) {
                        std::cout << ps.to_string() << " " << ns.to_string() << std::endl;
                        distance.push_back(moleSet[boost::lexical_cast<int>(ns.moleId)].getInterval(ps.startSite));
                    }
                }
            }
            std::cout << "distance size: " << distance.size() << std::endl;
            if(distance.size() != 0) {
                //distance average is not a fair value
                contig.push_back(std::accumulate(distance.begin(), distance.end(), 0.0) / distance.size());
                for(int dis : distance) {
                    std::cout << dis << " ";
                }
                std::cout << std::endl;
            }

        }
        contigs.push_back(contig);
    }
    return true;
}

bool ContigBuilder::build(const std::string& input, const std::string& output, double minScore, int threads) const {
    std::vector<Mole> moleSet;
    if(boost::filesystem::exists(input)) {
        std::ifstream moleInstream(input.c_str());
        std::string line;
        int count = -1;
        while(std::getline(moleInstream, line)) {
            //using count as moleName is convince, but is hard to dubug.
            count ++;
            std::string lineName = std::to_string(count);
            std::getline(moleInstream, line);
            boost::trim(line);
            std::vector<std::string> lineSplit;
            boost::split(lineSplit, line, boost::is_any_of(" "), boost::token_compress_on);
            std::vector<int> tempMole;
            BOOST_FOREACH(std::string interval, lineSplit) {
                tempMole.push_back(boost::lexical_cast< int > (interval));
            }
            Mole m(lineName, tempMole);
            moleSet.push_back(m);
        }
    }
    std::vector<Alignment> alignments;
    alignment(moleSet, alignments, threads, minScore);
    Graph g;
    edgeSet edges;
    constructGraph(alignments, g, edges);
    Components comps;
    bfsSearch(g, edges, 2, comps);

    //divide
    Components dividedComponents;
    Components tempComp;
    
    for(auto comp : comps) {
        tempComp.clear();
        divide(comp, g, tempComp);
        LOG4CXX_DEBUG(logger, boost::format("divide into %d components.") % tempComp.size());
        for(auto cp : tempComp) {
            dividedComponents.push_back(cp);
        }
    }
    std::vector<std::vector<int>> contigs;
    connect(dividedComponents, moleSet, contigs);
    for(int i = 0; i < contigs.size(); ++ i) {
        std::cout << "contig_" << i << std::endl;
        for(int j = 0; j < contigs[i].size(); ++ j) {
            std::cout << contigs[i][j] << " ";
        }
        std::cout << std::endl;
    }
    LOG4CXX_DEBUG(logger, boost::format("there are %d components.") % contigs.size());
}

