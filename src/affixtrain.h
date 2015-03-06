#ifndef AFFIXTRAIN_H
#define AFFIXTRAIN_H

class optionStruct;

class countAndWeight
    {
    int Nnodes;
    double Weight;
    public:
        void setNnodes(int N){Nnodes = N;};
        void setWeight(double W){Weight = W;};
        int getNnodes()const{return Nnodes;}
        double getWeight()const{return Weight;}
    };


//void trainRules(const char * tag, optionStruct * options);
void trainRules(const char * tag, optionStruct * options,countAndWeight * Counts);

extern int openfiles;

#endif
