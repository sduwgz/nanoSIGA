#ifndef contig_builder_h_
#define contig_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <boost/lexical_cast.hpp>


struct Vertex {
    //moleId|startSite
    std::string moleId;
    int startSite;
    Vertex(const std::string& mid, const int ssite) : moleId(mid), startSite(ssite) {}
    bool operator < (const Vertex& v) const {
        int moleNum1 = boost::lexical_cast<int>(moleId);
        int moleNum2 = boost::lexical_cast<int>(v.moleId);

        return moleNum1 * 1000 + startSite < moleNum2 * 1000 + v.startSite;
    }
    std::string to_string() const {
        std::string sv = "(" + moleId + ", " + std::to_string(startSite) + ")";
        return sv;
    }
};

typedef std::map<std::string, double> edgeSet;
typedef std::map<Vertex, std::set<Vertex>> Graph;

class ContigBuilder {
public:
    ContigBuilder(Map& maptool, std::string prefix="default") : _maptool(maptool), _prefix(prefix) {
    }
    void start(const std::vector<Mole>* moleSetPtr, std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId) const;
    void alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, int threadNum, double minScore) const;
void constructGraph(const std::vector<Alignment>& alignments, Graph& graph1, Graph& graph2, edgeSet& edges) const;
    bool build(const std::string& input, const std::string& output, const std::string& alignmentFile, double minScore, int threads=16) const;
    
private:
    Map _maptool;
    std::string _prefix;
};
#endif //contig_builder_h_
