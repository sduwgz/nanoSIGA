#ifndef contig_builder_h_
#define contig_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <map>
#include <set>

struct Vertex {
    //moleId|startSite
    std::string moleId;
    int startSite;
    Vertex(const std::string& mid, const int ssite) : moleId(mid), startSite(ssite) {}
    bool operator < (const Vertex& v) const {
        std::string sa = moleId + "|" + std::to_string(startSite);
        std::string sv = v.moleId + "|" + std::to_string(v.startSite);
        return sa < sv;
    }
};


class ContigBuilder {
public:
    ContigBuilder(Map& maptool, std::string prefix="default") : _maptool(maptool), _prefix(prefix) {
    }
    void alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, double minScore) const;
    bool build(const std::string& input, const std::string& output, double minScore, int threads=1) const;
    
private:
    Map _maptool;
    std::string _prefix;
};
#endif //contig_builder_h_
