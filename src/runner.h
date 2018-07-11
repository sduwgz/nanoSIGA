#ifndef runner_h_
#define runner_h_

#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>

typedef boost::property_tree::ptree Properties;
typedef std::vector< std::string > Arguments;

class Runner {
public:
    const std::string options() const {
        return _options;
    }
    std::string transform(const char key) const {
        auto it = _transform.find(key);
        if(it != _transform.end()) {
            return it->second;
        }
        return std::string(1, key);
    }
    virtual int run(const Properties options, const Arguments& arguments) = 0;
protected:
    Runner(const std::string& options = "", const std::map<char, std::string>& table=std::map<char, std::string>()) : _options(options), _transform(table) {}
    std::string _options;
    std::map<char, std::string> _transform;
};

typedef Runner* RunnerPtr;
class RunnerManager {
public:
    virtual ~RunnerManager() {
    }
    static RunnerManager* instance() {
        static RunnerManager mgr;
        return &mgr;
    }
    bool install(const std::string& name, RunnerPtr ptr) {
        if(_runners.find(name) != _runners.end())
            return false;
        _runners.insert(std::make_pair(name, ptr));
        return true;
    }
    bool uninstall(const std::string& name) {
        auto it = _runners.find(name);
        if(it == _runners.end()) {
            return false;
        }
        _runners.erase(it);
        return true;
    }
    RunnerPtr get(const std::string name) const {
        if(_runners.find(name) != _runners.end()) {
            return _runners.find(name)->second;
        }
        return NULL;
    }
    int help(int argc, char* argv[]) const {
        std::cout << "Usage: XXX" << std::endl;
    }
private:
    RunnerManager() {
    }
    std::map<std::string, RunnerPtr> _runners;
};
#endif //runner_h_
