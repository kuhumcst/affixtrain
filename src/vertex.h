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

#ifndef VERTEX_H
#define VERTEX_H

#include "settingsaffixtrain.h"
#include "strng.h"
#include <stdio.h>


class trainingPair;
class shortRulePair;
class hashTable;

class vertex;

extern int VertexPointerCount;

struct topScore
    {
    vertex * best;
    long MaxR;
    long MaxW;
    long MaxN;
    };

class vertex
    {
    public:
        strng Pattern;
        strng Replacement;
        double wght;
    private:
        friend class hashTable;
        vertex * Next; // vertices with the same hash key.
        vertex ** Head;
        hashTable * Hash;
    public:
        long R__R;
        long R__W;
        long W__R;
        long W__W;
#if _NA
        long R__NA;
        long W__NA;
#endif
    private:
        int RefCount;
        int Relations;
#if PRUNETRAININGPAIRS
        int RuleLikes;
#endif
    public:
        matchResult apply(trainingPair * trainingpair);
        matchResult applym(trainingPair * pair, char * mask);
#if _NA
        void adjustNotApplicableCountsByRecalculatingR_NA(trainingPair * NotApplicableRight,int total);
        void adjustNotApplicableCountsByRecalculatingW_NA(trainingPair * NotApplicableWrong,int total);
#endif
        void adjustWeight();
#if DOIMPEDANCE
        double impedance; 
        /* Regard the part of the pattern that is enclosed in * as a series of condensors.
        The capacity of each condensor is equal to the length of the pattern string.

        For example, the impedance of ^tvö*k*nur$ is 1/len("k") = 1
        impedance(^*ê*à*ête$) = (1/len("ê")+1/len("à")) = 2
        impedance(^a*ge*etste$) = 1/2
        impedance(^ge*etste$) = 0
        Hypothesis: Rules with low impedance are better (less fluffy) than rules with high impedance.
        */
        void computeImpedance();
#endif
        void incRelations(){++Relations;}
        void decRelations(){--Relations;}
        int relations(){return Relations;}
        void print(FILE * f,int level);
        void printRule(FILE * f,int level,int nr);
        void print1(FILE * f)const;
        vertex * getNext(){return Next;}
        void setNext(vertex * next){Next = next;}
        void setHead(vertex ** head){Head = head;}
        char * itsTxt(){return Pattern.itsTxt();}
        vertex * findReplacement(vertex * Rule)
            {
            vertex * p = this;
            while (p && (!p->Replacement.eq(Rule->Replacement) || !p->Pattern.eq(Rule->Pattern)))
                p = p->Next;
            return p;
            }
        void destroy();
        void incRefCnt(){++RefCount;}
        void decRefCnt(){--RefCount;}
        int refCount(){return RefCount;}
        ~vertex();
        void deleteThis();
        void construct(const char * pat, const char * rep);
        vertex(vertex * Rule,hashTable * Hash);
        vertex(shortRulePair * Rule);
        vertex(const char * pat, const char * rep);
        matchResult lemmatise(trainingPair * pair);
        int goodness(trainingPair * pairs, topScore * Top);
        void nlemmatiseStart();
#if PRUNETRAININGPAIRS
        int ruleLikes()const{return RuleLikes;}
#endif
#if SMALLMEMORY
        int nlemmatise(trainingPair * pairs,int n,bool InputRight,int skip);
#else
        int nlemmatise(trainingPair * pairs,bool InputRight,int skip);
#endif
#if AMBIGUOUS
        void markAmbiguousForNextRound(trainingPair * pair);
#endif
    };

extern bool building;
class vertexPointer
    {
    private:
        vertex * V;
        vertexPointer * Next;
        bool InputRight;
        bool Right;
    public:
        vertexPointer(vertex * V,vertexPointer * Next,bool InputRight,bool Right);
        ~vertexPointer();
#if PRUNETRAININGPAIRS
        bool fewerLikesThan(int thresh) const;
#endif
        void deleteAll()
            {
            vertexPointer * p = this;
            while(p)
                {
                vertexPointer * n = p->Next;
                delete p;
                p = n;
                }
            }
        bool has(vertex * V)
            {
            for(vertexPointer * p = this; p; p = p->Next)
                if(p->V == V)
                    return true;
            return false;
            }

    };

#endif