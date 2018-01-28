#ifndef correct_builder_h_
#define correct_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <set>
#include <map>

typedef std::vector<std::vector<int> > Center;

struct Index {
    int moleId;
    int startSite;
    int endSite;
};

class CorrectBuilder {
public:
    CorrectBuilder(const std::string& prefix="default") : _prefix(prefix) {
    }
    bool build(const std::string& moleFile, const std::string& clusterFile, double minScore, const std::string& output) const;
    void alignment(const std::vector<Mole>& moleSet, const Index& center, std::vector<Alignment>& alignments) const;
    Center vote(std::vector<Alignment>& alignments) const;
    int center(const std::vector<Index>& indexes) const;
private:
    std::string _prefix;

};
#endif //correct_builder_h_
