#ifndef constant_h_
#define constant_h_
#include<string>

extern const char* kLogConfig;

//Parameters in DP
extern const double PI;
extern const int BIN;
extern const int MINCLUSTER; 
extern const int MISS;
extern const int MINMOLELEN;
extern const int MAXITER;
extern const int EPSILON;
extern const int MAXCLUSTER;
extern const int DEC;
extern const int BATCH;
//Length of integral interval
extern const int Interval;

//Number of enzyme cleavage site
extern const int MAX_MISS_MATCH;
extern const int MIN_MATCH_NUMBER;

//Initial value of scoring matrix
extern const double INIT_SCORE;

//Parameter in deletion model
extern const int UNIT_LENGTH;
extern const int MAX_DELETION; 

//Chip id
extern const std::string CHIP_ID;
extern const std::string CMP_HEADER;

//k-gram
extern const int K;
#endif //constant_h_

