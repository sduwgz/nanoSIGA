#ifndef overlap_builder_h_
#define overlap_builder_h_

#include <string>
#include <vector>

#include "map.h"
#include "index.h"

class OverlapBuilder {
public:
    OverlapBuilder(const std::string& parameterFile, const std::string& prefix) : _maptool(parameterFile), _prefix(prefix) {
    }
    bool build(const std::string& input, double minScore, const std::string& output, size_t threads, int trim, bool reverse, bool useHash);
void start(std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId, int batch);
void alignment(std::ofstream& overlapOutstream, int threads, double minScore, int trim);
void buildIndex() {
    _index.build(_moleSet);
}
private:
    std::string _prefix;
    std::vector<Mole> _moleSet;
    Map _maptool;
    Index _index;
};
#endif //overlap_builder_h_
