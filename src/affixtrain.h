#ifndef AFFIXTRAIN_H
#define AFFIXTRAIN_H

#include "settingsaffixtrain.h"

class optionStruct;

class countAndWeight
    {
    int Nnodes;
    double Weight;
    LONG CountByDepth;
    public:
        void setNnodes(int N){Nnodes = N;};
        void setWeight(double W){Weight = W;};
        void setCountByDepth(double C){CountByDepth = C;};
        int getNnodes()const{return Nnodes;}
        double getWeight()const{return Weight;}
        LONG getCountByDepth()const{return CountByDepth;}
    };


//void trainRules(const char * tag, optionStruct * options);
void trainRules(const char * tag, optionStruct * options,countAndWeight * Counts);

extern int openfiles;

#endif
