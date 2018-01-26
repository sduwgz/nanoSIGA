#include "cluster_builder.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::getLogger("nanoARCS.cluster"));

bool ClusterBuilder::build(const std::string& input, double minScore, int minCluster, const std::string& output) const {
    
}
