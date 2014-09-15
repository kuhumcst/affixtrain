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

extern int ambivalentWords;
extern int alternatives;
extern int allwords;


extern int NodeCount;

extern int VertexPointerCount;
extern int VertexCount;
extern int TrainingPairCount;
extern int HashCount;
extern int RulePairCount;
extern int RuleTemplateCount;
extern int ShortRulePairCount;
extern int FullRulePairCount;
//extern int XX;



edif dif(char * Txt, char * s_Txt);

class vertex;
class trainingPair;
class hash;
class vertexCount;
class ruleTemplate;


const int b_ambiguous     = 1 << 0;
const int b_doublet       = 1 << 1;
const int b_test          = 1 << 2;
const int b_ok            = 1 << 3; // if AMBIGUOUS is set TRUE, b_ok and 
const int b_wrong         = 1 << 4; // b_wrong are used during build-up as well
const int b_skip          = 1 << 5;
const int b_oksibling     = 1 << 6;
const int b_bench         = 1 << 7; // Homographs that are not lemmatised 
                                    // correctly but that have a correctly
                                    // lemmatised sibling get this flag and 
                                    // are ignored for the remainder of the
                                    // training.


class vertexPointer;

class trainingPair
    {
    private:
        const char * Word;
        const char * LemmaHead;
#if WORDCLASS
        const char * WordClass;
#endif
#if LEMMACLASS
        const char * LemmaClass;
#endif
        size_t wordlength;
        size_t lemmalength;
#if WORDCLASS
        size_t wordclasslength;
#endif
#if LEMMACLASS
        size_t lemmaclasslength;
#endif
#if LEMMAINL
        int Inl;
        int Lemma_Inl;
#endif
        trainingPair * Next;
        char * Mask; /* Initially 0. When trainingpair becomes member of a
                     node's Wrong of Right list, it borrows Mask from the
                     trainingPairPointer pointing to this trainingPair. */
        char * Lemma; // as computed
        vertex * V; // the rule that made Lemma
        int bits;
        vertexPointer * applicableRules;
        void deleteRules();
    public:
#if AMBIGUOUS
        trainingPair * Alt; // forms closed loop of alternatives (of lemmas)
#endif
#if PESSIMISTIC
        trainingPair * AltLemma; // forms closed loop of inflected forms (of same lemma)
#endif
#if AMBIGUOUS
//        void makeLoop(trainingPair ** first,trainingPair ** ceiling);
        void makeWrongAmbiguousIfRightPresent(trainingPair *& Wrong,trainingPair *& Ambiguous);
#endif
        void addRule(vertex * V,bool InputRight,bool Right);
        void allDeleteRules();
        trainingPair * nth(int n);
        void makeCandidateRules(hash * Hash,vertex * parent,bool alreadyRight);
        int makeRuleEx(hash * Hash,vertex * parent,bool alreadyRight);
        int makeCorrectRules(hash * Hash,ruleTemplate * Template,const char * similar,vertex * parent, int mlow, int recurse);
        // if recurse = 0, don't recurse. Else recurse. If a succeeding rule
        // already has been found, decrement recurse and then test whether a
        // recursive call must be made.
        trainingPair ** pNext()
            {
            return &Next;
            }
        void setNext(trainingPair * n)
            {
            Next = n;
            }
        const char * getMask() const
            {
            return Mask;
            }
        void setLemma(char * lemma)
            {
            if(Lemma)
                delete [] Lemma;
            Lemma = dup(lemma);
            }
        void unsetLemma()
            {
            delete [] Lemma;
            Lemma = 0;
            }
        void setVertex(vertex * v)
            {
            V = v;
            }
        int count()
            {
            int n = 0;
            trainingPair * t = this;
            while(t)
                {
                ++n;
                t = t->Next;
                }
            return n;
            }
        void takeMask(char * mask)
            {
            if(Mask)
                delete [] Mask;
            Mask = dup(mask);
            }
        trainingPair * next() const
            {
            return Next;
            }
        static void makeChains(int allPairs,trainingPair * TrainingPair,trainingPair ** train,trainingPair ** test);
        static int makeNextTrainingSet(int allPairs,trainingPair * TrainingPair,FILE * train,FILE * done,FILE * combined, FILE * disamb);
        void set(int bit)
            {
            bits |= bit;
            }
        void unset(int bit)
            {
            bits &= ~bit;
            }
        int isset(int bit)
            {
            return bit & bits;
            }
        void fprint(FILE * fp);
        void fprintAll(FILE * fp);
        void fprintTraining(FILE * fp
#if WRITEINVERTED
            ,FILE * fpinverted
#endif
            );
        const char * itsWord(){return Word;}
        const char * itsLemmaHead(){return LemmaHead;}
#if LEMMACLASS
        const char * itsLemmaClass(){return LemmaClass;}
        size_t itsLemmaClassLength(){return lemmaclasslength;}
#endif
        size_t itsWordlength(){return wordlength;}
        size_t itsLemmalength(){return lemmalength;}
        void init(const char * Word
            ,size_t wordlength
            ,const char * LemmaHead
            ,size_t lemmalength
#if WORDCLASS
            ,const char * WordClass
            ,size_t wordclasslength
#endif
#if LEMMACLASS
            ,const char * LemmaClass
            ,size_t lemmaclasslength
#endif
#if LEMMAINL
            ,int Inl,int Lemma_Inl
#endif
            )
            {
            this->Word = Word;
            this->wordlength = wordlength;
            this->LemmaHead = LemmaHead;
            this->lemmalength = lemmalength;
#if WORDCLASS
            this->WordClass = WordClass;
            this->wordclasslength = wordclasslength;
#endif
#if LEMMACLASS
            this->LemmaClass = LemmaClass;
            this->lemmaclasslength = lemmaclasslength;
#endif

            this->bits = 0;
#if LEMMAINL
            this->Inl = Inl;
            this->Lemma_Inl = Lemma_Inl;
#endif
#if AMBIGUOUS
            this->Alt = this;
#endif
#if PESSIMISTIC
            this->AltLemma = this;
#endif
            }
        bool isCorrect(const char * lemma) const;
        int printAll(FILE * f,const char * h,int s);
        void print(FILE * f);
        void printMore(FILE * f);
        void printSep(FILE * f);
        trainingPair():Next(0),Mask(0),Lemma(0),V(0),applicableRules(0)
            {
            ++TrainingPairCount;
            }
        ~trainingPair();
        int cmpWord(const trainingPair * B) const;
        int cmpLemma(const trainingPair * B) const;
#if WORDCLASS
        int cmpWordClass(const trainingPair * B) const;
#endif
#if LEMMACLASS
        int cmpLemmaClass(const trainingPair * B) const;
#endif
#if LEMMAINL
        int cmpFreq(const trainingPair * B) const;
#endif
        trainingPair * append(trainingPair * list)
            {
            trainingPair * p = this;
            while(p->Next)
                {
                p = p->Next;
                }
            p->Next = list;
            return this;
            }
    };

class rulePair
    {
    protected:
        virtual char * itsPattern() = 0;
        virtual char * itsReplacement() = 0;
    public:
        const char * pattern()
            {
            return itsPattern();
            }
        const char * replacement()
            {
            return itsReplacement();
            }
        bool apply(trainingPair * trainingpair,size_t lemmalength,char * lemma,char * mask);
        edif dif(rulePair * other);
        rulePair(){++RulePairCount;}
        virtual ~rulePair(){--RulePairCount;}
    };

class trainingPairPointer;

typedef enum {failure,wrong,right} matchResult;

struct topScore
    {
    vertex * best;
    //bool cut;
    long MaxR;
    long MaxW;
    long MaxN;
    };

class vertex : public rulePair
    {
    private:
        friend class hash;
        strng * Pattern;
        strng * Replacement;
        vertex * Next; // vertices with the same hash key.
        vertex ** Head;
        int RefCount;
        int Relations;
        hash * Hash;
    public:
        long R__R;
        long R__W;
        long W__R;
        long W__W;
#if _NA
        long R__NA;
        long W__NA;
        /* checkPV
        long R(){return R__R+R__W+R__NA;}
        long W(){return W__R+W__W+W__NA;}
        */
#endif
    public:
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
        virtual char * itsPattern(){return Pattern->itsTxt();}
        virtual char * itsReplacement(){return Replacement->itsTxt();}
        virtual const strng * itsstrngPattern()const{return Pattern;}
        virtual const strng * itsstrngReplacement()const{return Replacement;}
        void incRelations(){++Relations;}
        void decRelations(){--Relations;}
        int relations(){return Relations;}
        void print(FILE * f,int level);
        void printRule(FILE * f,int level,int nr);
        void print1(FILE * f);
        vertex * getNext()
            {
            return Next;
            }
        void setNext(vertex * next)
            {
            Next = next;
            }
        void setHead(vertex ** head)
            {
            Head = head;
            }
        char * itsTxt()
            {
            return Pattern->itsTxt();
            }
        strng * replacement()
            {
            return Replacement;
            }
        vertex * findReplacement(rulePair * Rule)
            {
            vertex * p = this;
            while(p && (!p->Replacement->eq(Rule->replacement()) || !p->Pattern->eq(Rule->pattern())))
                p = p->Next;
            return p;
            }
        void destroy();
        void incRefCnt()
            {
            ++RefCount;
            }
        void decRefCnt()
            {
            --RefCount;
            }
        int refCount()
            {
            return RefCount;
            }
        ~vertex();
        void deleteThis();
        vertex(rulePair * Rule,hash * Hash);
        matchResult lemmatise(trainingPair * pair,char ** mask,char ** plemma);
        int goodness(trainingPair * pairs,topScore * Top/*,bool maycut*/);
        void nlemmatiseStart();
        int nlemmatise(trainingPair * pairs,int n,bool InputRight);
    };

extern bool building;
class vertexPointer
    {
    private:
        vertex * V;
        vertexPointer * Next;
        bool InputRight:1;
        bool Right:1;
    public:
        vertexPointer(vertex * V,vertexPointer * Next,bool InputRight,bool Right)
            : V(V),Next(Next),InputRight(InputRight),Right(Right)
            {
            if(InputRight)
                {
                if(Right)
                    ++(V->R__R);
                else
                    ++(V->R__W);
                }
            else
                {
                if(Right)
                    ++(V->W__R);
                else
                    ++(V->W__W);
                }
            ++VertexPointerCount;
            }
        ~vertexPointer()
            {
            if(building)
                {
                if(InputRight)
                    {
                    if(Right)
                        --(V->R__R);
                    else
                        --(V->R__W);
                    }
                else
                    {
                    if(Right)
                        --(V->W__R);
                    else
                        --(V->W__W);
                    }
                }
            --VertexPointerCount;
            }
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
    };


    class node // node in decision tree
        {
        private:
            bool compatibleSibling(node * sib);
        public:
            void Counting(int & nodes,int & pairs,FILE * f)
                {
                int Nodes = 0;
                int Pairs = 0;
                ++Nodes;
                if(this->Right)
                    Pairs += this->Right->count();
                if(this->IfPatternSucceeds)
                    this->IfPatternSucceeds->Counting(Nodes,Pairs,f);
                if(this->IfPatternFails)
                    this->IfPatternFails->Counting(Nodes,Pairs,f);
                fprintf(f,"%d\t%d\t%f\t%f\n",Nodes,Pairs,log((double)Nodes),log((double)Pairs));
                nodes += Nodes;
                pairs += Pairs;
                }
            vertex * V;
           // vertex ** pv; // array of candidates of the IfPatternFails list, more or less sorted in descending goodness order
            int Ncandidates;
            int IDXcandidates;
            node * IfPatternSucceeds; 
            node * IfPatternFails;
            trainingPair * Right; // Pattern succeeds and replacement is right.
            void printSep(FILE * f,int level)
                {
                this->V->print(f,level);
                if(Right)
                    Right->printSep(f);
                if(this->IfPatternSucceeds)
                    this->IfPatternSucceeds->printSep(f,level+1);
                if(this->IfPatternFails)
                    this->IfPatternFails->printSep(f,level);
                }
#if AMBIGUOUS
            //trainingPair * Ambiguous; // Pattern succeeds, but replacement is wrong. However, the same word has more lemmas, one of which was produced by the rule.
#endif

            matchResult lemmatise(trainingPair * pair);
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
            double weightedcount()
                {
                double ret = 0.0;
                node * n = this;
                while(n)
                    {
                    if(this->Right)
                        {
                        double rcount = (double)this->Right->count();
                        assert(rcount >= 0.9);
                        ret += 5.0*(rcount-1.0)*exp(-0.9*rcount)+1; 
                        /* penalizes as follows:
                        rcount >> 1 :  1
                        rcount ~ 2 : 1.75
                        rcount == 1 : 1

                        So it penalizes the rules with the lowest number of 
                        supporting word/lemma examples that we aren't going to
                        cut away. The rules that are going to be cut away are
                        less extremely penalized, but enough to make it 
                        disadvantageous to split the heaviest penalized rules
                        in even lower grade rules.
                        */
                        }
                    if(n->IfPatternSucceeds)
                        ret += n->IfPatternSucceeds->weightedcount();
                    n = n->IfPatternFails;
                    }
                return ret;
                }
            int countWords()
                {
                int ret = 0;
                node * n = this;
                while(n)
                    {
                    if(n->Right)
                        ret += n->Right->count();
                    if(n->IfPatternSucceeds)
                        ret += n->IfPatternSucceeds->countWords();
                    n = n->IfPatternFails;
                    }
                return ret;
                }
            int prune(int threshold)
                {
                int N = 0;
                while(IfPatternFails)
                    {
                    N = IfPatternFails->prune(threshold);
                    if(0 <= N && N <= threshold)
                        {
                        node * child = IfPatternFails;
                        IfPatternFails = child->IfPatternFails;
                        child->IfPatternFails = 0;
                        delete child;
                        continue;
                        }
                    break;
                    }
                while(IfPatternSucceeds)
                    {
                    N = IfPatternSucceeds->prune(threshold);
                    if(0 <= N && N < threshold)
                        {
                        node * child = IfPatternSucceeds;
                        IfPatternSucceeds = child->IfPatternFails;
                        child->IfPatternFails = 0;
                        delete child;
                        continue;
                        }
                    break;
                    }
                if(IfPatternSucceeds)
                    return -1;
                else
                    N = Right->count();
                return N;
                }

            int pruneAll(int threshold)
                {
                //printf("prune\n");
                int n = prune(threshold);
                //printf("prune DONE\n");
                return n;
                }

            int print(FILE * fo,int ind,int & Nnodes,int &NnodesR)
                {
                V->print(fo,ind);
                int n = 0;
                ++Nnodes;
                if(Right)
                    {
                    ++NnodesR;
                    if(IfPatternSucceeds)
                        n += Right->printAll(fo,"\t\tR+\t",',');
                    else
                        n += Right->printAll(fo,"\t\tR\t",',');
                    }
                if(IfPatternSucceeds)
                    n += IfPatternSucceeds->print(fo,ind+1,Nnodes,NnodesR);
                if(IfPatternFails)
                    n += IfPatternFails->print(fo,ind,Nnodes,NnodesR);
                return n;
                }
           

#if AMBIGUOUS
            void init(trainingPair ** allRight,trainingPair ** allWrong/*,trainingPair ** allAmbiguous*/,int level);
#else
            void init(trainingPair ** allRight,trainingPair ** allWrong,int level/*,vertex ** pvp,int Np*/);
#endif
            void splitTrainingPairList(trainingPair * all,trainingPair **& pNotApplicable,trainingPair **& pWrong,trainingPair **& pRight);
            node(vertex * V):V(V),IfPatternSucceeds(0),IfPatternFails(0),Right(0)
                {
                ++NodeCount;
                V->incRefCnt();
                }
            ~node()
                {
                --NodeCount;
                delete IfPatternSucceeds;
                delete IfPatternFails;
//                delete [] pv;
               // V->decRefCnt();
                V->destroy();
                }
        };

typedef  void (vertex::*forallfunc)();

class hash
    {
    private:
        long hash_size;
        long record_count;
        vertex ** hash_table;
        hash * next;
    public:
        vertex ** convertToList(int & N);
        hash * getNext()
            {
            return next;
            }
        hash *& getNextR()
            {
            return next;
            }
        void inc()
            {
            ++record_count;
            }
        void dec()
            {
            --record_count;
            }
        void rehash(int loadFactor/*1-100*/);
        long loadfactor()
            {
            if(record_count < 10000000L)
                return (100 * record_count) / hash_size;
            else
                return record_count / (hash_size/100);
            }
        hash(long size);
        ~hash();
        long key(const char * ckey);
        vertex * find(const char * ckey,vertex **& head);
        int forall(forallfunc fnc);
        vertex * getVertex(rulePair * Rule,bool & New);
        bool deleteVertex(rulePair * Rule);
    };



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
    );


#endif
