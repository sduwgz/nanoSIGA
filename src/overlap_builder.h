#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "map.h"
#include <string>
#include <vector>

class OverlapBuilder {
public:
    OverlapBuilder(Map& maptool, const std::string& prefix="default") : _maptool(maptool), _prefix(prefix) {
    }
    bool build(const std::string& input, double minScore, const std::string& output, size_t threads=1, int trim=0, int reverseLabel=0) const;

private:
    Map _maptool;
    std::string _prefix;
    void alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, int threads, double minScore) const;
    void start(const std::vector<Mole>* moleSetPtr, std::vector<Alignment>* alignmentsPtr, double minScore, int threads, int threadId) const;
};
#endif //overlap_builder_h_
