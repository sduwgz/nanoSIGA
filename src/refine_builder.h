#ifndef contig_builder_h_
#define contig_builder_h_

#include "map.h"
#include "mole.h"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <boost/lexical_cast.hpp>

class RefineBuilder {
public:
    RefineBuilder(const std::string parameter_file="parameter_file") : _parameter_file(parameter_file) {
    }
    bool build(const std::string& input, const std::string& moleFile, const std::string& output) const;
    
private:
    std::string _parameter_file;
};
#endif //refine_builder_h_
