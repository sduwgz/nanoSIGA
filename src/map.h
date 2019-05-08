#ifndef map_h_
#define map_h_

#include "mole.h"
#include "constant.h"

#include <vector>
#include <string>
#include <map>

typedef std::vector<double> PunishScore;
typedef std::map<std::string, double> ParametersList;

class Map {
public:
    Map() {}
    Map(const std::string& parameter_file) : _parameterfile(parameter_file) {
        if(this->initParameters(parameter_file))
            this->initPunishScore();
    };

    static Map* instance(const std::string& parameterFile) {
        static Map maptool;
        if(maptool.initParameters(parameterFile)) {
            maptool.initPunishScore();
        } else {
        }
        return &maptool;
    }

    double wholeDPscore(const std::vector<int>& d1, const std::vector<int>& d2) const;
    Alignment localDPscore(const Mole& m1, const Mole& m2) const;
    double validScore(const Fragment& moleFragment, const Fragment& geneFragment) const;
private:
    bool initParameters(const std::string& parameter_file);
    void initPunishScore();

    double probLaplace(int delta) const;
    double probDeletion(int siteNumber, int moleLength) const;
    double probInsertion(int k) const;
    double probBackground(int delta) const;

    ParametersList _parameters;
    PunishScore _insertionScore;
    PunishScore _deletionScore;
    PunishScore _laplaceScore;
    PunishScore _backgroundScore;
    std::string _parameterfile;
};
#endif  //map_h_
