#ifndef overlap_builder_h_
#define overlap_builder_h_

#include "map.h"
#include "index.h"
#include <string>
#include <vector>

class OverlapBuilder {
public:
    OverlapBuilder(const std::string& parameterFile, const std::string& prefix="default") : _parameterFile(parameterFile), _prefix(prefix) {
    }
    bool build(const std::string& input, double minScore, const std::string& output, size_t threads=1, int trim=1, bool reverseLabel=false, bool useHash=false) const;
private:
    std::string _parameterFile;
    std::string _prefix;
};
#endif //overlap_builder_h_
