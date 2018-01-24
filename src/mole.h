#ifndef mole_h_
#define mole_h_

#include <iostream>
#include <vector>
#include <string>


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
    void trim() {
        
    }
};

class Mole {
public:
    Mole() {} 
    explicit Mole(std::string id) : _id(id) {
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
    friend class MoleReader;
private:
    std::string _id;
    std::vector<long> _position;
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
