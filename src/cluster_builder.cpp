#include "cluster_builder.h"

#include <queue>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <log4cxx/logger.h>


static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.cluster"));

void constructGraph(const std::string& input, double minScore, Graph& graph, edgeSet& edges) {
    std::ifstream overlapIn(input.c_str());
    std::string line;
    while(std::getline(overlapIn, line)) {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        double score = boost::lexical_cast<double>(data[2]);
        if(score < minScore) continue;
        Vertex v1 = data[0] + "|" + data[3];
        Vertex v2 = data[1] + "|" + data[4];
        edges[v1 + v2] = score;
        edges[v2 + v1] = score;
        graph[v1].insert(v2);
        graph[v2].insert(v1);
    }
}
void bfsSearch(const Graph& graph, edgeSet& edges, int minCluster, Components& comps) {
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
        if(comp.size() == 2) {
            std::string edge;
            for(auto it : comp) {
                edge += it;
            }
            if(edges[edge] < 30)
                continue;
        }
        count ++;
        LOG4CXX_DEBUG(logger, boost::format("#%d cluster, cluster size is %d.") % count % comp.size());
        comps.push_back(comp);
    }
}

bool ClusterBuilder::build(const std::string& input, double minScore, int minCluster, const std::string& output) const {
    Graph graph;
    edgeSet edges;
    if(boost::filesystem::exists(input)) {
        constructGraph(input, minScore, graph, edges);
    } else {
        LOG4CXX_WARN(logger, boost::format("%s not is not existed.") % input);
        return false;
    }
    LOG4CXX_DEBUG(logger, boost::format("%d verteies in graph.") % graph.size());
    Components comps;
    bfsSearch(graph, edges, minCluster, comps);
    LOG4CXX_DEBUG(logger, boost::format("%d cluster, minimal cluster size is %d.") % comps.size() % minCluster);
    std::ofstream os(output.c_str());
    for(auto comp : comps) {
        for(auto it = comp.begin(); it != comp.end(); ++ it) {
            os << *it << " ";
        } 
        os << '\n';
    }
    return true;
}
