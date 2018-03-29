#ifndef mole_h_
#define mole_h_

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>


typedef std::vector<int> Fragment;
typedef std::vector<Fragment> AlignedMole;

struct Alignment {
    double score;
    std::string mole1Id;
    std::string mole2Id;
    int mole1Start, mole2Start;
    int mole1End, mole2End;
    AlignedMole alignedMole1, alignedMole2;
    std::ostream& print(std::ostream& os, bool distail=false) const {
        os << mole1Id << " " << mole2Id << " " << score << " ";
        os << mole1Start << " " << mole2Start << " " << mole1End << " " << mole2End << '\n';
        if(distail) {
            for(int i = 0; i < alignedMole1.size(); ++ i) {
                for(int j = 0; j < alignedMole1[i].size() - 1; ++ j) {
                    os << alignedMole1[i][j] << ",";
                }
                os << alignedMole1[i][alignedMole1[i].size() - 1] << " ";
            }
            os << '\n';
            for(int i = 0; i < alignedMole2.size(); ++ i) {
                for(int j = 0; j < alignedMole2[i].size() - 1; ++ j) {
                    os << alignedMole2[i][j] << ",";
                }
                os << alignedMole2[i][alignedMole2[i].size() - 1] << " ";
            }
            os << '\n';

        }
        return os;
    }
    void trimHead() {
        while(abs(std::accumulate(alignedMole1[0].begin(), alignedMole1[0].end(), 0) - std::accumulate(alignedMole2[0].begin(), alignedMole2[0].end(), 0)) > 1500 || alignedMole1[0].size() > 1 || alignedMole2[0].size() > 1) {
            mole1Start += alignedMole1[0].size();
            mole2Start += alignedMole2[0].size();
            alignedMole1 = std::vector<Fragment>(alignedMole1.begin() + 1, alignedMole1.end());
            alignedMole2 = std::vector<Fragment>(alignedMole2.begin() + 1, alignedMole2.end());
        }
    }
    void trimTail() {
        while(abs(std::accumulate(alignedMole1.back().begin(), alignedMole1.back().end(), 0) - std::accumulate(alignedMole2.back().begin(), alignedMole2.back().end(), 0)) > 1500 || alignedMole1.back().size() > 1 || alignedMole2.back().size() > 1) {
            mole1End -= alignedMole1.back().size();
            mole2End -= alignedMole2.back().size();
            alignedMole1.pop_back();
            alignedMole2.pop_back();
        }
    }
};

class Mole {
public:
    Mole() {} 
    explicit Mole(std::string& id) : _id(id) {
    }
    Mole(std::string& id, std::vector<int>& dis) : _id(id), _distance(dis) {
    }
    virtual ~Mole() {}
    bool getDistance(); 
    Mole reverseMole();    
    std::vector<int> getDistance() const {
       return _distance;
    } 
    std::string getID() const {
       return _id;
    } 
    void setID(const std::string& s) {
       _id = s;
    } 
    int getInterval(const int index) const {
        if(index > _distance.size()) return 0;
        return _distance[index];
    }
    void shear(const int startSite) {
        _distance = std::vector<int>(_distance.begin() + startSite, _distance.end());
    }
    friend class MoleReader;
private:
    std::string _id;
    std::vector<long long> _position;
    std::vector<int> _distance;
};

class MoleReader {
public:
    MoleReader(std::istream& stream) : _stream(stream) {};
    bool read(Mole& mole);
    void reset(Mole& mole) {
        mole._id = "";
        mole._distance.clear();
        mole._position.clear();
    }

private:
    std::istream& _stream;
};
#endif //mole_h_
