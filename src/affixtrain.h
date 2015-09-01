#ifndef AFFIXTRAIN_H
#define AFFIXTRAIN_H

#include "settingsaffixtrain.h"

class optionStruct;

class countAndWeight
    {
    int Nnodes;
    double Weight;
    LONG CountByDepth;
    LONG CountBySize;
    public:
        void setNnodes(int N){Nnodes = N;};
        void setWeight(double W){Weight = W;};
        void setCountByDepth(LONG C){CountByDepth = C;};
        void setCountBySize(LONG C){CountBySize = C;};
        int getNnodes()const{return Nnodes;}
        double getWeight()const{return Weight;}
        LONG getCountByDepth()const{return CountByDepth;}
        LONG getCountBySize()const{return CountBySize;}
    };


//void trainRules(const char * tag, optionStruct * options);
void trainRules(const char * tag, optionStruct * options,countAndWeight * Counts);

extern int openfiles;

#endif
