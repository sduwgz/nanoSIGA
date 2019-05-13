#include<unordered_map>
#include<set>
#include<iostream>
#include<string>

#include "mole.h"

typedef std::pair<std::string, int> Value;
typedef std::unordered_map<std::string, std::set<Value>> HashTable;
class Index {
public:
    Index() {}
    Index(int k) : K(k) {}
    void build(const std::vector<Mole>& moleSet);
    std::set<std::string> query(const Mole& m) const;
private:
    int K = 3;
    HashTable htable;
};
