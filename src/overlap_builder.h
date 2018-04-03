#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "map.h"
#include <string>
#include <vector>

class OverlapBuilder {
public:
    OverlapBuilder(const std::string& parameterFile, const std::string& prefix="default") : _parameterFile(parameterFile), _prefix(prefix) {
    }
    bool build(const std::string& input, double minScore, const std::string& output, size_t threads=1, int trim=0, int reverseLabel=0) const;

private:
    std::string _parameterFile;
    std::string _prefix;
    void alignment(const std::vector<Mole>& moleSet, std::vector<Alignment>& alignments, int threads, double minScore) const;
};
#endif //overlap_builder_h_
