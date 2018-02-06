#ifndef contig_builder_h_
#define contig_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <map>

typedef std::string Vertex;
typedef std::map<Vertex, std::set<Vertex>> Graph;
typedef std::vector<std::set<Vertex>> Components;

class ContigBuilder {
public:
    ContigBuilder(Map& maptool, const::string& prefix="default") : _maptool(maptool), _prefix(prefix) {
    }
    bool build(const std::string& input, const std::string& output, const double minScore, int threads=1) const;
    void alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignmentSet);
    void component(const std::vector<Alignment>& alignmentSet, Components& components) const;
    void bfsSearch(const Graph& graph, edgeSet& edges, int minCluster, Components& coms) const;
    
private:
    Map _maptool;
    std::string _prefix;
};
#endif //contig_builder_h_

