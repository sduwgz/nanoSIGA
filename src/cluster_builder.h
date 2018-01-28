#ifndef cluster_builder_h_
#define cluster_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <set>
#include <map>

typedef std::string Vertex;
typedef std::map<Vertex, std::set<Vertex>> Graph;
typedef std::vector<std::set<Vertex>> Components;

class ClusterBuilder {
public:
    ClusterBuilder(const std::string& prefix="default") : _prefix(prefix) {
    }
    bool build(const std::string& input, double minScore, int minCluster, const std::string& output) const;
    void constructGraph(const std::string& input, double minScore, Graph& graph) const;
    void bfsSearch(const Graph& graph, int minCluster, Components& coms) const;
private:
    std::string _prefix;

};
#endif //cluster_builder_h_
