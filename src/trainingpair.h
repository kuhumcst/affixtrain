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
#ifndef TRAININGPAIR_H
#define TRAININGPAIR_H

#include "settingsaffixtrain.h"
#include <stdio.h>

class vertex;
class vertexPointer;
class hashTable;
class optionStruct;
class ruleTemplate;

extern FILE * fpmourn;

enum eResult {undecided, yes, no
#if AMBIGUOUS
                                , notme
#endif
    };

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
    public:
        trainingPair * Next;
#if AMBIGUOUS
        trainingPair * Alt; // forms closed loop of alternatives (of lemmas)
#endif
#if PESSIMISTIC
        trainingPair * AltLemma; // forms closed loop of inflected forms (of same lemma)
#endif
    private:
        char * Mask; /* Initially 0. When trainingpair becomes member of a
                     node's Wrong of Right list, it borrows Mask from the
                     trainingPairPointer pointing to this trainingPair. */
        char * Lemma; // as computed
        vertex * V; // the rule that made Lemma
        vertexPointer * applicableRules;
        unsigned int bits:8;
        unsigned int ambs:3;
        unsigned int tentativeAmbs:3;
        // For prepruning, pairs with fewer than 3 or 2 likes can removed from
        // training.
#if PRUNETRAININGPAIRS
        unsigned int likes:2; // 0:do not train with this pair
                              // 1:this pair has a unique inflection (initial, pessimistic value)
                              // 2:this pair has an almost unique inflection
                              // 3:this pair has not an almost unique inflection
#endif
        void deleteRules();
    public:
#if AMBIGUOUS
        trainingPair * makeWrongAmbiguousIfRightPresent(trainingPair *& Ambiguous);
#endif
        void addRule(vertex * V,bool InputRight,bool Right);
        void allDeleteRules();
        trainingPair * nth(int n);
        void makeCandidateRules(hashTable * Hash,vertex * parent,bool alreadyRight,optionStruct * options);
        int makeCorrectRules(hashTable * Hash,ruleTemplate * Template,const char * similar,vertex * parent, int mlow, int recurse,optionStruct * options);
        // if recurse = 0, don't recurse. Else recurse. If a succeeding rule
        // already has been found, decrement recurse and then test whether a
        // recursive call must be made.
#if PRUNETRAININGPAIRS
        void like(unsigned int L)
            {
            if(0 < likes && likes < L)
                likes = L > 3 ? 3 : L;
            }
        bool fewerLikesThan(int thresh) const;
        int notLemmatizedBy(vertex * V);
        void mourn();
#endif
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
        void setLemma(char * lemma);
        void unsetLemma()
            {
            delete [] Lemma;
            Lemma = 0;
            }
        void setVertex(vertex * v)
            {
            V = v;
            }
        unsigned long count()
            {
            unsigned long n = 0;
            trainingPair * t = this;
            while(t)
                {
                ++n;
                t = t->Next;
                }
            return n;
            }
        void takeMask(char * mask);
        trainingPair * next() const
            {
            return Next;
            }
        inline void set(int bit)
            {
            bits |= bit;
            }
        inline void unset(int bit)
            {
            bits &= ~bit;
            }
        inline unsigned int isset(unsigned int bit)
            {
            return bit & bits;
            }
        inline void setRes(eResult r)
            {
            ambs = r;
            }
        inline unsigned int getRes()
            {
            return ambs;
            }
        inline void setTentativeRes(eResult r)
            {
            tentativeAmbs = r;
            }
        inline unsigned int getTentativeRes()
            {
            return tentativeAmbs;
            }
        bool checkResAll(eResult what)
            {
            bool ret = false;
            trainingPair * p;
            for(p = this;p;p = p->next())
                {
                if(p->getRes() == (unsigned int)what)
                    {
                    ret = true;
                    printf("CHECK:");
                    p->print(stdout);
                    printf("\n");
                    }
                }
            return ret;
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
        void printMore(FILE * f)const;
        void printSep(FILE * f);
        trainingPair();
        ~trainingPair();
        void deepDelete();
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

#endif

