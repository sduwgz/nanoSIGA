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

bool ClusterBuilder::build(const std::string& input, double minScore, int minCluster, const std::string& output) const {
    Graph graph;
    if(boost::filesystem::exists(input)) {
        constructGraph(input, minScore, graph);
    } else {
        LOG4CXX_WARN(logger, boost::format("%s not is not existed.") % input);
        return false;
    }
    LOG4CXX_DEBUG(logger, boost::format("%d verteies in graph.") % graph.size());
    Components comps;
    bfsSearch(graph, minCluster, comps);
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
void ClusterBuilder::constructGraph(const std::string& input, double minScore, Graph& graph) const {
    std::ifstream overlapIn(input.c_str());
    std::string line;
    while(std::getline(overlapIn, line)) {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        double score = boost::lexical_cast<double>(data[2]);
        if(score < minScore) continue;
        Vertex v1 = data[0] + "|" + data[3] + "|" + data[5];
        Vertex v2 = data[1] + "|" + data[4] + "|" + data[6];
        graph[v1].insert(v2);
        graph[v2].insert(v1);
    }
} 
void ClusterBuilder::bfsSearch(const Graph& graph, int minCluster, Components& comps) const {
    std::map<Vertex, int> colors;
    for(auto it : graph) {
        colors[it.first] = 0;
    }
    std::queue<Vertex> q;
    std::set<Vertex> comp;
    
    for(auto it : graph) {
        if(colors[it.first] != 0) continue;
        comp.clear();
        q.push(it.first);
        while(!q.empty()) {
            Vertex tmp = q.front();
            colors[tmp] = 2;
            comp.insert(tmp);
            for(auto iv : it.second) {
                if(colors[iv] == 0) {
                    q.push(iv);
                    colors[iv] = 1;
                }
            }
        }
        comps.push_back(comp);
    }
}
