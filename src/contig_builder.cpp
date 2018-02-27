#include "contig_builder.h"

#include <iostream>
#include <fstream>
#include <queue>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <log4cxx/logger.h>

typedef std::map<Vertex, std::set<Vertex>> Graph;
typedef std::vector<std::set<Vertex>> Components;
typedef std::map<std::string, double> edgeSet;
typedef std::map<int, std::set<int> > SiteGraph;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.contig"));

void ContigBuilder::alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, double minScore) const {
    for(int i = 0; i < moleSet.size(); ++ i) {
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
        comps.push_back(comp);
    }
}
void divide(std::set<Vertex>& component, const Graph& graph) {
    //there are some mixed components, divide them into several components.
    std::map<std::string, std::vector<Vertex> > vertexInMole;
    for(Vertex vertex : component) {
        vertexInMole[vertex.moleId].push_back(vertex);
    }
    std::map<int, int> vertexCount;
    for(auto it : vertexInMole) {
        vertexCount[it.second.size()] ++;
    }
    int maxCount = -1;
    int maxCountVertex = -1;

    for(auto it : vertexCount) {
        if(maxCount > it.second) {
            maxCount = it.second;
            maxCountVertex = it.first;
        }
    }
    std::vector<std::set<Vertex>> dividedComponents(maxCountVertex, std::set<Vertex>());
    for(auto it : vertexInMole) {
        if(it.second.size() == maxCountVertex) {
            for(int i = 0; i < maxCountVertex; ++ i) {
                dividedComponents[i].insert(it.second[i]);
            }
        }
    }
}

void topoSort(SiteGraph& spreGraph, SiteGraph& spostGraph, std::vector<int>& spath) {
    //Toposort maybe drop in infinite loops
    while(spostGraph.size() != 0) {
        int tempSite = -1;
        for(auto site : spostGraph) {
            if(spostGraph[site.first].size() == 0) {
                tempSite = site.first;
                break;
            }
        }
        if(tempSite == -1 && spostGraph.size() != 0) {
            LOG4CXX_INFO(logger, boost::format("there is a circle in siteGraph."));
            break;
        }
        spath.push_back(tempSite);
        spostGraph.erase(tempSite);
        for(int site : spreGraph[tempSite]) {
            spostGraph[site].erase(tempSite);
        }
    }
}

int nextSites(const Components& comps, const std::set<Vertex>& curcomp) {
    //TODO: how to slove wrong connection
    int nextcomp = -1;
    std::vector<int> distances;
    int maxOne = 0;
    for(int i = 0; i < comps.size(); ++ i) {
        auto comp = comps[i];
        distances.clear();
        for(Vertex ns : comp) {
            for(Vertex ps : curcomp) {
                if(ns.moleId == ps.moleId && ns.startSite > ps.startSite) {
                    distances.push_back(ns.startSite - ps.startSite);
                }
                if(distances.size() == 0) continue;
                if(std::count(distances.begin(), distances.end(), 1) > maxOne) {
                    maxOne = std::count(distances.begin(), distances.end(), 1);
                    nextcomp = i;
                }
            }
        }
    }
    return nextcomp;
}

bool connect(const Components& comps, std::vector<Mole>& moleSet, std::vector<int>& contig) {
    SiteGraph spreGraph, spostGraph;
    for(int i = 0; i < comps.size(); ++ i) {
        auto curcomp = comps[i];
        int nextcomp = nextSites(comps, curcomp);
        if(nextcomp != -1) {
            spreGraph[i].insert(nextcomp);
            spostGraph[nextcomp].insert(i);
        }
    }
    std::vector<int> spath;
    topoSort(spreGraph, spostGraph, spath);
    if(spath.size() != spreGraph.size()) {
        LOG4CXX_INFO(logger, boost::format("there is a circle in siteGraph."));
        return false;
    }
    //dynamic programming for longest path
    std::vector<int> dpLength(spath.size(), 0);
    std::vector<int> dpPath(spath.size(), -1);
    dpLength[spath[0]] = 0;
    dpPath[spath[0]] = -1;
    for(int i = 1; i < spath.size(); ++ i) {
        for(int postSite : spostGraph[spath[i]]) {
            int weight = 0;
            for(Vertex ps : comps[postSite]) {
                for(Vertex ns : comps[spath[i]]) {
                    if(ns.moleId == ps.moleId && ns.startSite - ps.startSite == 1) {
                        weight ++;
                    }
                }
            }
            if(dpLength[postSite] + weight > dpLength[spath[i]]) {
                dpLength[spath[i]] = dpLength[postSite] + weight;
                dpPath[spath[i]] = postSite;
            }
        }
    }
    //trace back
    int pathEnd = -1;
    int maxLength = -1;
    std::vector<int> contigPath;
    for(int i = 0; i < spath.size(); ++ i) {
        if(dpLength[i] > maxLength) {
            maxLength = dpLength[i];
            pathEnd = dpPath[i];
        }
    }
    while(pathEnd != -1) {
        contigPath.push_back(pathEnd);
        pathEnd = dpPath[pathEnd];
    }
    for(int i = contigPath.size() - 1; i >= 0; -- i) {
        std::vector<int> distance;
        for(Vertex ps : comps[i]) {
            for(Vertex ns : comps[i + 1]) {
                if(ns.moleId == ps.moleId && ns.startSite - ps.startSite == 1) {
                    distance.push_back(moleSet[boost::lexical_cast<int>(ns.moleId) - 1].getInterval(ps.startSite));
                }
            }
        }
        if(distance.size() != 0) {
            //distance average is not a fair value
            contig.push_back(std::accumulate(distance.begin(), distance.end(), 0.0) / distance.size());
        }
    }
    return true;
}
bool ContigBuilder::build(const std::string& input, const std::string& output, double minScore, int threads) const {
    std::vector<Mole> moleSet;
    if(boost::filesystem::exists(input)) {
        std::ifstream moleInstream(input.c_str());
        std::string line;
        int count = 0;
        while(std::getline(moleInstream, line)) {
            count ++;
            std::string lineName = line;
            std::getline(moleInstream, line);
            std::vector<std::string> lineSplit;
            boost::split(lineSplit, line, boost::is_any_of(" "), boost::token_compress_on);
            std::vector<int> tempMole;
            std::string moleName = std::to_string(count);
            BOOST_FOREACH(std::string interval, lineSplit) {
                tempMole.push_back(boost::lexical_cast< int > (interval));
            }
            Mole m(moleName, tempMole);
            moleSet.push_back(m);
        }
    }
    std::vector<Alignment> alignments;
    alignment(moleSet, alignments, minScore);
    Graph g;
    edgeSet edges;
    constructGraph(alignments, g, edges);
    Components comps;
    bfsSearch(g, edges, 2, comps);
    
}

