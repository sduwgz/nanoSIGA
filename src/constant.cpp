#include "constant.h"

const char* kLogConfig = "log4cxx.properties";
//Parameters in DP
const double PI = 3.1415926;
const int BIN = 1000;
const int MINCLUSTER = 20; 
const int MISS = 3;
const int MINMOLELEN = 2;
const int MAXITER = 100;
const int EPSILON = 1;
const int MAXCLUSTER = 10000;
const int DEC = 1;

//
const int BATCH = 1000;

//Length of integral interval
const int Interval = 100;

//Number of enzyme cleavage site
const int MAX_MISS_MATCH = 3;
const int MIN_MATCH_NUMBER = 3;

//Initial value of scoring matrix
const double INIT_SCORE = -100000.0;

//Parameter in deletion model
const int UNIT_LENGTH = 10000;
const int MAX_DELETION = 20; 

//Chip id
const std::string CHIP_ID = "20245,11541,3/19/2013,0009155";
//Cmp header
const std::string CMP_HEADER = "# CMAP File Version:    0.1\n# Label Channels:   1\n# Nickase Recognition Site 1:   GCTCTTC\n# Enzyme1:  Nt.BspQI\n# Number of Consensus Nanomaps: 1\n#h CMapId   ContigLength    NumSites    SiteID  LabelChannel    Position    StdDev  Coverage    Occurrence\n#f int  float   int int int float   float   int int\n";
