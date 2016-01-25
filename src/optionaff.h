/*
AFFIXTRAIN - supervised learning of affix rules for CSTLEMMA

Copyright (C) 2012  Center for Sprogteknologi, University of Copenhagen

This file is part of AFFIXTRAIN.

AFFIXTRAIN is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

AFFIXTRAIN is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AFFIXTRAIN; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef OPTIONAFF_H
#define OPTIONAFF_H

#include <stdio.h> // FILE

typedef enum {GoOn = 0,Leave = 1,Error = 2} OptReturnTp;
typedef enum {econstant = 0, edepth = 1, esupport = 2, eentropy = 3, esize = 4} OptWeightFunction;

void cleanUpOptions();

class optionStruct
    {
    private:
        double * D; // list of penalty parameters
                    // R__R;W__R;R__W;W__W
                    // or
                    // R__R;W__R;R__W;W__W;R__NA;W__NA
                    // Defaults to 0;0;1;0;0;0
        int nD; // set by counting numbers in D, defaults to 6
        int c; // cutoff
        int C; // expected cutoff
        const char * e; // extra
        const char * f; // compfunc
        const char * i; // word list
        const char * j; // temp dir
        const char * n; // columns
        const char * o; // flexrules
        const char * B; // Best parms
        const char * P; // Current parms
        const char * b; // raw rules
        const char * G; // External training program
        const char * E; // External lemmatizer
        const char * X; // rule weight function, defaults to constant (1)
        const char * I; // Input (with option -b: additionally lemmatise a word list - one word per line)
        const char * O; // Output (lemmas of input -I) (with option -b)
        bool ComputeParms;// compute parms
        bool SuffixOnly;// suffix only
        bool ExpensiveInfix;/*to create intermediate between full affix
                            and suffix only (which is better than affix,
                            sometimes, for unclear reasons.)*/
        bool Verbose;   // verbose
        bool Test;
        bool TrainTest;
        bool Remove;    // remove test files after use (default false)
        bool VX; // 10-fold cross validation
        double Minfraction;
        double Maxfraction;
        int K;    // Number of differently sized fractions of trainingdata 
        double M; // # Iterations when training with Maxfraction of input
        double N; // # Iterations when training with Minfraction of input
        bool Redo; /* option R. Set to true if training has to be done once
                   more after the removal of the homographs that will be
                   handled in the next iteration.*/
        int Q; // Max recursion depth when attempting to create candidate rule
  /*      int q;*/ /* Percentage of training pairs to set aside for testing set to
               positive value if you want to set PERC percent of the avaliable
               data aside for testing.*/
        bool F; // Create flexrules. Can be combined with computation (-p) and
                // testing (-t, -T).
        /* End of options. */

        int Blobs; // Number of blobs found in word list
        int Lines; // Number of lines found in word list
        int FracBlobs; // Number of blobs used for training
        int FracLines; // Number of lines used for training
        bool SuffixOnlyParmSeen;
        int Swath;
        int SwathIteration;
        int NumberOfNodes;
        int TrainingPairsLines;
        double TreePenalty;
//        double DepthWeight;
        char * Argstring;
        OptWeightFunction WeightFunction;
        OptReturnTp doSwitch(int c, char * locoptarg, char * progname);
        OptReturnTp readOptsFromFile(char * locoptarg,char * progname);
        void detectFloatingPointNumbers(const char * S);
    public:
        const OptWeightFunction getWeightFunction(){return WeightFunction;}
        const int blobs() const { return Blobs; }
        const int lines() const { return Lines; }
        const int fracblobs() const { return FracBlobs; }
        const int fraclines() const { return FracLines; }
        void setReadBlobs(int bl) { FracBlobs = bl; }
        void setReadLines(int li) { FracLines = li; }
        const int cutoff() const{ return c; }
        const int expectedCutoff() const{return C;}
        const char * extra() const{return e ? e : "";}
        const char * compfunc() const{return f;}
        const char * wordPerLineList() const{return I;}
        const char * wordLemmaList() const{return i;}
        const char * tempDir() const{return j;}
        const char * columns() const{return n;}
        //const char * flexrules() const{return o;}
        const char * flexrules();
        const char * lemmas() const{return O;}
        const char * bestParms() const{return B;}
        const char * currentParms() const{return P;}
        const char * ruleWeightFunction() const{return X;}
        const char * rawRules() const{return b;}
        const char * externalTrainer() const{return G;}
        const char * externalLemmatizar() const{return E;}
        const int maxRecursionDepthRuleCreation()const{return Q;}
//        const int percentageTestPairs()const{return q;}
        const bool computeParms()const{return ComputeParms;}
        const bool suffixOnly()const{return SuffixOnly;}
        const bool expensiveInfix(){ return ExpensiveInfix; }
        const bool verbose()const{return Verbose;}
        const bool remove()const{ return Remove; }
        const bool tenfoldCrossValidation()const{ return VX; }
        const bool test()const{return Test;}
        const bool trainTest()const{return TrainTest;}        
        const double minfraction()const{return Minfraction;}
        const double maxfraction()const{return Maxfraction;}
        const bool redo()const{return Redo;}
        const int swaths()const{return K;}
        const double minIterations()const{return M;}
        const double maxIterations()const{return N;}
        const int numberOfParms()const{return nD;}
        const double parm(int i)const{return (0<=i&&i<nD) ? D[i] : 0.0;}
        void setParm(int i,double v){if(0<=i&&i<nD)D[i]=v;}
        const bool createFlexRules(){return F;}
        void setSwath(int S){Swath = S;}
        void setSwathIteration(int SI){SwathIteration = SI;}
        void setNumberOfNodes(int NoN){NumberOfNodes = NoN;}
        void setTrainingPairsLines(int TPL){TrainingPairsLines = TPL;}
        void setWeight(double W){TreePenalty = W;}

        void setI(const char * WordList);
        void seti(const char * WordLemmaList);
        void setO(const char * Result);
        void seto(const char * Result);
        void sete(const char * Extra);
        void setn(const char * Columns);
        void setf(const char * Compfunc);
        void setP(const char * ParamFile);
        void setc(int cutoff){c = cutoff;}

        const char * argstring() const { return Argstring; }
        const char * argstring_no_path() const;
        void setArgstring();
        void print(FILE * fp) const;
        void printArgFile() const;
        void printEvaluation(const char * introduction,char * evaluation,char * postScriptum) const;        
        optionStruct();
        optionStruct(optionStruct & O);
        ~optionStruct();
        void completeArgs();
        OptReturnTp readArgs(int argc, char * argv[]);
    };
#endif
