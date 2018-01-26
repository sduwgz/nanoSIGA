#ifndef cluster_builder_h_
#define cluster_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>

class ClusterBuilder {
public:
    ClusterBuilder(const std::string& prefix="default") : _prefix(prefix) {
    }
    bool build(const std::string& input, double minScore, int minCluster) const;
private:
    std::string _prefix;

};
#endif //overlap_builder_h_
