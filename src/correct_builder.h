#ifndef correct_builder_h_
#define correct_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <set>
#include <map>

typedef std::vector<int> Center;

struct Index {
    std::string moleId;
    int startSite;
    int endSite;
    Index(const std::string  a, int b, int c) : moleId(a), startSite(b), endSite(c) {}
    inline std::string toString() const {
        std::string tmp = moleId + "|" + std::to_string(startSite) + "|" + std::to_string(endSite);
        return tmp;
    }
};
inline bool operator==(const Index& index1, const Index& index2) {
    return index1.moleId == index2.moleId && index1.startSite == index2.startSite && index1.endSite == index2.endSite;
}
inline bool operator!=(const Index& index1, const Index& index2) {
    return !(index1 == index2);
}

class CorrectBuilder {
public:
    CorrectBuilder(Map& maptool, const std::string& prefix="default") : _maptool(maptool), _prefix(prefix) {
    }
    bool build(const std::string& moleFile, const std::string& clusterFile, double minScore, const std::string& output, int thread=1) const;
    void alignment(const std::map<std::string, Mole>& moleSet, const std::vector<Index>& indexes, const Index& voteCenter, std::vector<Alignment>& alignments) const;
    void divide(const std::vector<Alignment>& alignments, std::map<int, std::vector<Alignment>>& dividedAlignments) const;
    Center vote(std::vector<Alignment>& alignments) const;
    Index center(const std::vector<Index>& indexes) const;
private:
    Map _maptool;
    std::string _prefix;

};
#endif //correct_builder_h_
