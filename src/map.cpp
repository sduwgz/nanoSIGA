#include "map.h"
#include "constant.h"

#include <math.h>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("nanoARCS.map"));

double Map::validScore(const Fragment& moleFragment, const Fragment& geneFragment) const {
    int moleLength = 0, geneLength = 0;
    moleLength = std::accumulate(moleFragment.begin(), moleFragment.end(), 0);
    geneLength = std::accumulate(geneFragment.begin(), geneFragment.end(), 0);

    int delta = moleLength - geneLength + 100000;
    if (delta < 0) {
        delta = 0;
    }
    if (delta > 200000) {
        delta = 200000;
    }

    int moleSiteNumber = moleFragment.size() - 1;
    int geneSiteNumber = geneFragment.size() - 1;

    int deleteNumber = static_cast< int > ((geneSiteNumber + 0.0) / moleLength * UNIT_LENGTH + 0.5);
    if (deleteNumber < 1) {
        deleteNumber = 1;
    } else if (deleteNumber > MAX_DELETION) {
        deleteNumber = MAX_DELETION;
    }

    if (moleSiteNumber != 0 || geneSiteNumber != 0) {
        return _laplaceScore[delta] + _deletionScore[deleteNumber] + _insertionScore[moleSiteNumber] - _backgroundScore[delta];
        //return probLaplace(delta) + probDeletion(geneSiteNumber, moleLength) + probInsertion(moleSiteNumber) - probBackground(delta);
    } else {
        return _laplaceScore[delta]  - _backgroundScore[delta];
        //return probLaplace(delta)  - probBackground(delta);
        //return laplace(delta) + pI(0) + pD(0, moleLength) - background(delta);
    }
}

double Map::probDeletion(int siteNumber, int moleLength) const {
    int deleteNumber = static_cast< int > ((siteNumber + 0.0) / moleLength * UNIT_LENGTH + 0.5);
    if (deleteNumber < 1) {
        deleteNumber = 1;
    } else if (deleteNumber > MAX_DELETION) {
        deleteNumber = MAX_DELETION;
    }
    double lambda = _parameters.find("lambda_poisson")->second;
    boost::math::poisson_distribution<> p(lambda);
    return log(boost::math::pdf(p, deleteNumber));
}

double Map::probInsertion(int k) const {
    double lambda = _parameters.find("lambda_exponent")->second;
    boost::math::exponential_distribution<> e(lambda);
    if (k == 0) {
        return log(cdf(e, 0.5) - cdf(e, 0));
    } 
    return log(cdf(e, k + 0.5) - cdf(e, k - 0.5));
    //return log(pow(lambda, k));
}

double Map::probBackground(int delta) const {
    double mu = _parameters.find("mu_background")->second;
    double sigma = _parameters.find("sigma_background")->second;
    boost::math::normal_distribution<> n(0, sigma);
    int distance = delta - mu;
    int d = distance / Interval;
    int interval_left = d * Interval;
    int interval_right = (d + 1) * Interval;
    if (distance < 0) {
        interval_left = (d - 1) * Interval;
        interval_right = d * Interval;
    }
    return log(boost::math::cdf(n, interval_right) - boost::math::cdf(n, interval_left));
}

double Map::probLaplace(int delta) const {
    double mu = _parameters.find("mu_laplace")->second;
    double sigma = _parameters.find("sigma_laplace")->second;
    //boost::math::laplace_distribution<> l(mu, sigma);
    boost::math::laplace_distribution<> l(0, sigma);
    int distance = delta - mu;
    int d = distance / Interval;
    int interval_left = d * Interval;
    int interval_right = (d + 1) * Interval;
    if (distance < 0) {
        interval_left = (d - 1) * Interval;
        interval_right = d * Interval;
    }
    return log(boost::math::cdf(l, interval_right) - boost::math::cdf(l, interval_left));
}

Alignment Map::localDPscore(const Mole& m1, const Mole& m2) const {
    std::vector<int> d1 = m1.getDistance();
    std::vector<int> d2 = m2.getDistance();
    std::vector<Fragment> alignedMole1;
    std::vector<Fragment> alignedMole2;
    int rows = d1.size() + 1, cols = d2.size() + 1;
    LOG4CXX_DEBUG(logger, boost::format("Rows: %d, Cols: %d") %rows %cols);
    double scoreMatrix[rows][cols];
    int trackMatrix[rows][cols];
    for (int i = 0; i < rows; ++ i) {
        for (int j = 0; j < cols; ++ j) {
            scoreMatrix[i][j] = INIT_SCORE;
            trackMatrix[i][j] = -1;
        }
    }
    for (int i = 0; i < rows; ++ i) {
        scoreMatrix[i][0] = 0;
    }
    for (int j = 0; j < cols; ++ j) {
        scoreMatrix[0][j] = 0;
    }
    LOG4CXX_DEBUG(logger, boost::format("Init is successful."));
    for (int i = 1; i < rows; ++ i) {
        for (int j = 1; j < cols; ++ j) {
            Fragment moleFragment, geneFragment;
            //one interval on mole is mapped with (1 -- K) intervals on gene
            moleFragment.push_back(d1[i - 1]);
            for (int k = j - 1; k >= j - MAX_MISS_MATCH && k >= 0; -- k) {
                geneFragment.push_back(d2[k]);
                double temp = scoreMatrix[i - 1][k] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                    trackMatrix[i][j] = (i - 1) * cols + k;
                }
            }

            //(1 -- K) intervals on mole is mapped with one interval on gene
            moleFragment.clear();
            geneFragment.clear();
            geneFragment.push_back(d2[j - 1]);
            for (int k = i - 1; k >= i - MAX_MISS_MATCH && k >= 0; -- k) {
                moleFragment.push_back(d1[k]);
                double temp = scoreMatrix[k][j - 1] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                    trackMatrix[i][j] = k * cols + j - 1;
                }
            }

            //two intervals on mole is mapped with two intervals on gene
            if (i > 2 && j > 2){
                moleFragment.clear();
                geneFragment.clear();
                geneFragment.push_back(d2[j - 1]);
                geneFragment.push_back(d2[j - 2]);
                moleFragment.push_back(d1[i - 1]);
                moleFragment.push_back(d1[i - 2]);

                double temp = scoreMatrix[i - 2][j - 2] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                    trackMatrix[i][j] = (i - 2) * cols + j - 2;
                }
            } 
        }
    }
    LOG4CXX_DEBUG(logger, boost::format("DP finished."));
    
    double maxScore = INIT_SCORE;
    int maxI = 0, maxJ = 0, startI = -1, startJ = -1;
     
    for (int i = 0; i < rows; ++ i) {
        if(scoreMatrix[i][cols - 1] > maxScore) {
            maxScore = scoreMatrix[i][cols - 1];
            maxI = i;
            maxJ = cols - 1;
        }
    }
    for (int j = 0; j < cols; ++ j) {
        if(scoreMatrix[rows - 1][j] > maxScore) {
            maxScore = scoreMatrix[rows - 1][j];
            maxI = rows - 1;
            maxJ = j;
        }
    }
    int trackI = maxI, trackJ = maxJ;
    while(trackI != 0 && trackJ != 0) {
        Fragment f1, f2;
        LOG4CXX_DEBUG(logger, boost::format("Track: %d %d %d") %trackI %trackJ %trackMatrix[trackI][trackJ]);
        int nextI = trackMatrix[trackI][trackJ] / cols;
        int nextJ = trackMatrix[trackI][trackJ] % cols;
        LOG4CXX_DEBUG(logger, boost::format("Next: %d %d") %nextI %nextJ);
        for(int i = nextI; i < trackI; ++ i) {
            //std::cout << d1[i - 1] << std::endl;
            f1.push_back(d1[i]);
        } 
        for(int j = nextJ; j < trackJ; ++ j) {
            f2.push_back(d2[j]);
        } 
        trackI = nextI, trackJ = nextJ;
        reverse(f1.begin(), f1.end());
        reverse(f2.begin(), f2.end());
        alignedMole1.push_back(f1);
        alignedMole2.push_back(f2);
    }
    reverse(alignedMole1.begin(), alignedMole1.end());
    reverse(alignedMole2.begin(), alignedMole2.end());
    LOG4CXX_DEBUG(logger, boost::format("Track finished."));
    int len1 = maxI - trackI;
    int len2 = maxJ - trackJ;
    //len1 compare with len2
    int startSiteI = trackI, startSiteJ = trackJ;
    //save the alignment in ret.
    Alignment ret;
    ret.score = maxScore;
    ret.mole1Id = m1.getID();
    ret.mole2Id = m2.getID();
    ret.mole1Start = trackI;
    ret.mole2Start = trackJ;
    ret.mole1End = maxI;
    ret.mole2End = maxJ;
    ret.alignedMole1 = alignedMole1;
    ret.alignedMole2 = alignedMole2;
    /*
    if(maxScore >= 15) {
        if(trackI == 0 && maxI == rows - 1 && trackJ == 0 && maxJ == cols - 1) {
            //std::cout << "equal " << k1 << " " << k2 << " ";
            std::cout << "equal " << k1 << "|" << startSiteI << "|" << rows << " " << k2 << "|" << startSiteJ << "|" << cols << " ";
        } else {
            if(trackI == 0 && maxI == rows - 1) {
                std::cout << "contain1 " << k1 << "|" << startSiteI << "|" << rows << " " << k2 << "|" << startSiteJ << "|" << cols << " ";
            }else if(trackJ == 0 && maxJ == cols - 1) {
                std::cout << "contain2 " << k1 << "|" << startSiteI << "|" << rows << " " << k2 << "|" << startSiteJ << "|" << cols << " ";
            } else {
                std::cout << "normal " << k1 << "|" << startSiteI << "|" << rows << " " << k2 << "|" << startSiteJ << "|" << cols << " ";
            }
        }
        //std::cout << trackI << " " << maxI << " " << trackJ << " " << maxJ << std::endl;
         * It is necessary to remove bad head and tail in an alignment in cluster and contig procession.
         * It is a good choise to remove bad head and tail in python. So, this part would not run.
         * But in vote, I must jump it to make sure every alignment has a same start site.
         *
        while(abs(accumulate((alignedMole1.back()).begin(), (alignedMole1.back()).end(), 0) - accumulate((alignedMole2.back()).begin(), (alignedMole2.back()).end(), 0)) > 1500 || alignedMole1.back().size() > 1 || alignedMole2.back().size() > 1) {
            startSiteI += alignedMole1.back().size();
            startSiteJ += alignedMole2.back().size();
            alignedMole1.pop_back();
            alignedMole2.pop_back();
        }

        if(trackI == 0) {
            int dis = accumulate(d2.begin(), d2.begin() + trackJ, 0);
            std::cout << k2 << "|" << startSiteJ << "|" << cols << " " << k1 << "|" << startSiteI << "|" << rows << " " << dis << " " << maxScore << std::endl;
        } else {
            int dis = accumulate(d1.begin(), d1.begin() + trackI, 0);
            std::cout << k1 << "|" << startSiteI << "|" << rows << " " << k2 << "|" << startSiteJ << "|" << cols << " " << dis << " " << maxScore << std::endl;
        }
        
    }
    */
    return ret;
}

double Map::wholeDPscore(const std::vector<int>& d1, const std::vector<int>& d2) const {
    int rows = d1.size() + 1, cols = d2.size() + 1;
    double scoreMatrix[rows][cols];
    for (int i = 0; i < rows; ++ i) {
        for (int j = 0; j < cols; ++ j) {
            scoreMatrix[i][j] = INIT_SCORE;
        }
    }
    scoreMatrix[0][0] = 0;
    for (int i = 1; i < rows; ++ i) {
        for (int j = 1; j < cols; ++ j) {
            Fragment moleFragment, geneFragment;
            moleFragment.push_back(d1[i - 1]);
            for (int k = j - 1; k >= j - MAX_MISS_MATCH && k >= 0; -- k) {
                geneFragment.push_back(d2[k]);
                double temp = scoreMatrix[i - 1][k] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                }
            }

            moleFragment.clear();
            geneFragment.clear();
            geneFragment.push_back(d2[j - 1]);
            for (int k = i - 1; k >= i - MAX_MISS_MATCH && k >= 0; -- k) {
                moleFragment.push_back(d1[k]);
                double temp = scoreMatrix[k][j - 1] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                }
            }

            if (i > 2 && j > 2){
                moleFragment.clear();
                geneFragment.clear();
                geneFragment.push_back(d2[j - 1]);
                geneFragment.push_back(d2[j - 2]);
                moleFragment.push_back(d1[i - 1]);
                moleFragment.push_back(d1[i - 2]);

                double temp = scoreMatrix[i - 2][j - 2] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                }
            } 
        }
    }
    return scoreMatrix[rows - 1][cols - 1];
}

bool Map::initParameters(const std::string& parameter_file) {
    boost::property_tree::ptree parameters;
    if (boost::filesystem::exists(parameter_file)) {
        try {
            boost::property_tree::read_ini(parameter_file, parameters);
        } catch (const boost::property_tree::ini_parser_error& e) {
            LOG4CXX_ERROR(logger, boost::format("load %s failed(%s).") % parameter_file % e.what());
            return false;
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % parameter_file);
        return false;
    }
    for (boost::property_tree::ptree::const_iterator it = parameters.begin(); it != parameters.end(); it++) {
        _parameters[it->first] = boost::lexical_cast<double> (it->second.data());
    }
    LOG4CXX_DEBUG(logger, boost::format("parameters mount: %d") % parameters.size());
    for (auto i = _parameters.begin(); i != _parameters.end(); ++ i) {
        LOG4CXX_DEBUG(logger, boost::format("%s: %s") % i->first % i->second);
    }
    return true;
}
void Map::initPunishScore() {
    for (int i = 0; i < MAX_MISS_MATCH + 1; ++ i) {
        _insertionScore.push_back(probInsertion(i));
    }
    for (int i = 0; i < MAX_DELETION + 1; ++ i) {
        _deletionScore.push_back(probDeletion(i, UNIT_LENGTH));
    }
    for (int i = 0; i <= 200000; ++ i) {
        _laplaceScore.push_back(probLaplace(i - 100000));
        _backgroundScore.push_back(probBackground(i - 100000));
    }
}
