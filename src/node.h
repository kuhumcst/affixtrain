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

#ifndef GRAPH_H
#define GRAPH_H

#include "strng.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stddef.h> // gcc ptrdiff_t
#include "settingsaffixtrain.h"

#if PRUNETRAININGPAIRS
extern FILE * fprune;
#endif

extern FILE * flog;
extern int parmsoff;
/*
parmsoff, the idea:

For weighting candidate rules, the program uses a number of rows of parameters,
curently 3. If the first row cannot decide which of two rules is best, the
second row is used, and so on. If, after using the last row of parameters, no
difference is found, that is the answer.

parsmoff is introduced to jump directly to the second row if the number of
candidate rules that still are in the race has fallen below a threshold, or to
the third row if that number is even lower than a second threshold.

The first threshold is 
                                                3/4
    thres34 = (initial number of candidate rules)

The second threshold is 
                                                1/2
    thres12 = (initial number of candidate rules)

thres34 and thres12 are rounded down to the nearest integer.
If the initial number of candidate rules is below 16, thres34 is set to 8 and
thres12 to 4.
*/

extern int visitno;
extern double wght;

#if AMBIGUOUS
extern int ambivalentWords;
extern int alternatives;
#endif
extern int allwords;


extern int NodeCount;
extern int RulePairCount;



class vertex;
class trainingPair;
class vertexCount;
class ruleTemplate;
class optionStruct;

const int b_doublet             = 1 << 0;
//const int b_test                = 1 << 1;
const int b_skip                = 1 << 2;
const int b_ambiguous           = 1 << 3;
const int b_solved              = 1 << 4;
const int b_tentativelysolved   = 1 << 5;


class vertexPointer;


class trainingPairPointer;
class shortRulePair;





class node // node in decision tree
    {
    private:
        bool compatibleSibling(node * sib);
    public:
        static int mostPenalized;
        void Counting(int & nodes,int & pairs,FILE * f);
        vertex * V;
        int Ncandidates;
        int IDXcandidates;
        node * IfPatternSucceeds; 
        node * IfPatternFails;
        trainingPair * Right; // Pattern succeeds and replacement is right.
        void printSep(FILE * f,int level);

        node * cleanup(node * parent);
        int count()
            {
            int ret = 0;
            node * n = this;
            while(n)
                {
                ++ret;
                if(n->IfPatternSucceeds)
                    ret += n->IfPatternSucceeds->count();
                n = n->IfPatternFails;
                }
            return ret;
            }
        LONG countByDepth(int Depth)
            {
            LONG ret = 0;
            for ( node * n = this
                ; n
                ; n = n->IfPatternFails
                )
                {
                if(n->Right)
                    {
                    ret += Depth;
                    }
                
                if(n->IfPatternSucceeds)
                    ret += n->IfPatternSucceeds->countByDepth(Depth+1);
                }
            return ret;
            }
        LONG countBySize();
        double weightedcount();
        double entropy(double Nnodes);
        int prune(int threshold);

        int pruneAll(int threshold)
            {
            int n = prune(threshold);
            return n;
            }

        int print(FILE * fo,int ind,int & Nnodes,int &NnodesR);

        void init(trainingPair ** allRight,trainingPair ** allWrong,int level,optionStruct * options);
        void splitTrainingPairList(trainingPair * all,trainingPair **& pNotApplicable,trainingPair **& pWrong,trainingPair **& pRight,bool SuffixOnly);
#if PRUNETRAININGPAIRS
        vertex ** cleanUpUnusedVertices(vertex ** pvf, vertex ** pvN
#if _NA
            , trainingPair * deletedRightPairs, trainingPair * deletedWrongPairs
#endif
            );
#endif
        node(vertex * V);
        ~node();
    };




class optionStruct;

int printRules
    ( node * nd
#if RULESASTEXTINDENTED
    , FILE * fo
#endif
#if BRACMATOUTPUT
    , FILE * fobra
#endif
    , FILE * folem
    , int ind
    , strng * L
    , strng * R
    , int & nr
    , optionStruct * options
    );


#endif
