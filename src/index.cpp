#include "index.h"
#include <map>

void Index::build(const std::vector<Mole>& moleSet) {
    for(Mole m : moleSet) {
        std::vector<int> distance = m.getData();
        if(distance.size() < 6) continue;
        for(int i = 0; i < distance.size() - K; ++ i) {
            std::string key;
            for(int j = i; j < i + K; ++ j) {
                key += std::to_string((distance[j] + 500) / 1000);
                key += "-";
            }
            htable[key].insert(std::make_pair(m.getID(), i));
        }
    }
}
/*
std::set<std::string> Index::query(const Mole& m) const {
    std::set<std::string> res;
    std::map<std::string, int> hitCount;
    std::vector<int> distance = m.getData();
    
    for(int i = 0; i < distance.size() - K; ++ i) {
        std::string key;
        for(int j = i; j < i + K; ++ j) {
            key += std::to_string((distance[j] + 500) / 1000);
            key += "-";
        }
        auto it = htable.find(key);
        if(it != htable.end()) {
            for(auto value : it->second) {
                hitCount[value.first] ++;
            }
        }
    }
    for(auto it : hitCount) {
        if(it.second >= 2)
            res.insert(it.first);
    }
    return res;
}
*/
std::set<std::string> Index::query(const Mole& m) const {
    std::set<std::string> res;
    std::vector<int> distance = m.getData();
    std::vector<std::set<Value>> hitList(distance.size(), std::set<Value>());
    for(int i = 0; i < distance.size() - K; ++ i) {
        std::string key;
        for(int j = i; j < i + K; ++ j) {
            key += std::to_string((distance[j] + 500) / 1000);
            key += "-";
        }
        auto it = htable.find(key);
        if(it != htable.end())
            hitList[i] = it->second;
    }
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> hitMolePos;
    for(int i = 0; i < hitList.size(); ++ i) {
        for(auto value : hitList[i]) {
            if(hitMolePos.find(value.first) != hitMolePos.end()) {
                for(auto molePos : hitMolePos[value.first]) {
                    int d1 = i - molePos.first;
                    int d2 = value.second - molePos.second;
                    if(d2 * d1 > 0 && abs(d1 - d2) < 4) {
                        res.insert(value.first);
                    }
                } 
            }
        }
        for(auto value : hitList[i]) {
            hitMolePos[value.first].push_back(std::make_pair(i, value.second));
        }
    }
    return res;
}
