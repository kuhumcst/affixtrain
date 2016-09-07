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

#include "vertex.h"
#include "node.h"
#include "trainingpair.h"
#include "shortrulepair.h"
#include "ruletemplate.h"


int VertexPointerCount = 0;
int VertexCount = 0;


#if DOIMPEDANCE
void vertex::computeImpedance()
    {
    /* Regard the part of the pattern that is enclosed in * as a series of condensors.
    The capacity of each condensor is equal to the length of the pattern string.

    For example, the impedance of ^tvö*k*nur$ is 1/len("k") = 1
    impedance(^*ê*à*ête$) = (1/len("ê")+1/len("à")) = 2
    impedance(^a*ge*etste$) = 1/2
    impedance(^ge*etste$) = 0
    Hypothesis: Rules with low impedance are better (less fluffy) than rules with high impedance.
    */
    impedance = 0.0;
    
    const char * p = Pattern.itsTxt();
    while(*p && *p != ANY)
        {
        //++p;
        p += skipUTF8char(p);
        --impedance; // characters at start (or end) diminish impedance.
        }
    if(*p)
        { // *p == ANY, the first star
        double len = 0.0;
        //while(*++p)
        int k;
        while((k = getUTF8char(p,UTF8)) != 0)
            {
            if(k == ANY) // the second, third, ... star
                {
                impedance += 1.0/len;
                len = 0.0;
                }
            else
                ++len;
            }
        impedance -= len; // characters at end diminish impedance.
        }
        
    }
#endif

void vertex::print1(FILE * f)const
    {
    fprintf(f,"{%s\t%s}",Pattern.itsTxt(),Replacement.itsTxt());
    }

void vertex::print(FILE * f,int level)
    {
        fprintf(f,"%*s|\t{%s\t%s}\n",level,"",Pattern.itsTxt(),Replacement.itsTxt());
    }

void vertex::printRule(FILE * f,int level,int nr)
    {
    fprintf(f,"%d\t%d\t%.*s\t%.*s\n",
        nr,
        level+1, // 1-based
        (int)(strlen(Pattern.itsTxt()) - 2),Pattern.itsTxt()+1,
        (int)(strlen(Replacement.itsTxt()) - 2),Replacement.itsTxt()+1);
    }


void vertex::destroy()
    {
    decRefCnt();
    int rc = refCount();
    if(rc == 0)
        {
        deleteThis();
        }
    }

vertex::~vertex()
    {
    --VertexCount;
    }

vertex::vertex(vertex * Rule,hashTable * Hash):
             Head(0)
            ,Hash(Hash)
            ,RefCount(0)
            ,Relations(0)
#if PRUNETRAININGPAIRS
            ,RuleLikes(0)
#endif
        ,R__R(0)
        ,R__W(0)
        ,W__R(0)
        ,W__W(0)
#if _NA
        ,R__NA(0)
        ,W__NA(0)
#endif

        {
        assert(Hash);
        Pattern.dup(Rule->Pattern.itsTxt());
        Replacement.dup(Rule->Replacement.itsTxt());
        ++VertexCount;
        }

matchResult vertex::lemmatise(trainingPair * pair)
    {
    matchResult ret = apply(pair);
    switch (ret)
        {
        case right:
            {
#if AMBIGUOUS
            for(trainingPair * q = pair->Alt;q != pair;q = q->Alt)
                {
                q->setTentativeRes(notme);
                q->set(b_tentativelysolved);
                }
#endif
            pair->setTentativeRes(yes);
            pair->set(b_tentativelysolved);
            break;
            }
        case wrong:
            switch(pair->getTentativeRes())
                {
                case yes:
                case undecided:
                case no:
                    {
#if AMBIGUOUS
                    for(trainingPair * q = pair->Alt;q != pair;q = q->Alt)
                        {
                        q->setTentativeRes(undecided);
                        }
#endif
                    pair->setTentativeRes(no);
                    break;
                    }
#if AMBIGUOUS
                case notme:
                    if(!pair->isset(b_tentativelysolved))
                        pair->setTentativeRes(no);
                    break;
#endif
                default:
                    //printf("?????????????????\n");
                    ;
                }
            assert(pair->getTentativeRes() != yes);
            break;
        default:
            ;
        }
    return ret;
    }

void vertex::nlemmatiseStart()
    {
    R__W = R__R = W__W = W__R = 0;
#if _NA
    R__NA = W__NA = 0;
#endif
    }

#if AMBIGUOUS
void vertex::markAmbiguousForNextRound(trainingPair * pair)
    {
    trainingPair * p;
    for(p=pair;p;p = p->next())
        {
        p->setRes(undecided);
        p->unset(b_tentativelysolved);
        }
    for(p=pair;p;p = p->next())
        {
        lemmatise(p);
        }
    }
#endif

int vertex::nlemmatise ( trainingPair * pair
#if SMALLMEMORY
                       , int n
#endif
                       , bool InputRight
                       , ptrdiff_t skip
                       )
    {
    int ret = 0;
    trainingPair * p;
    int skipped = 0;
#if SMALLMEMORY
    int m;
    for(p = pair, m = n; p && (m != 0);p = p->next(), --m)
#else
    for(p = pair; p; p = p->next())
#endif
        {
        p->setTentativeRes(undecided);
        p->unset(b_tentativelysolved);
        }
#if SMALLMEMORY
    for(p = pair, m = n; p && (m != 0);p = p->next(), --m)
#else
    for(p = pair; p; p = p->next())
#endif
        {
        ++skipped;
        if(skipped > skip)
            {
            skipped = 0;
            ++ret;
            switch(lemmatise(p))
                {
                case wrong:
#if AMBIGUOUS
                    assert(p->getTentativeRes() == no || p->getTentativeRes() == notme);
                    if(p->getTentativeRes() != notme)
                        { // don't count homographs that have an ok sibling
                        p->addRule(this,InputRight,false);
                        }
#else
                    p->addRule(this,InputRight,false);
#endif
                    break;
                case right:
                    assert(p->getTentativeRes() == yes);
                    // do opportunistically count homographs that are ok
                    p->addRule(this,InputRight,true);
#if PRUNETRAININGPAIRS
                    ++RuleLikes;
#endif
                    break;
#if _NA
                default:
                    if(InputRight)
                        ++R__NA;
                    else
                        ++W__NA;
#else
                default:
                    ;
#endif
                }
            }
        }
    assert(!p);
    return ret;
    }

void vertex::deleteThis()
    {
    if(Hash)
        {
        while(*Head && *Head != this)
            Head = &(*Head)->Next;
        if(*Head)
            {
            *Head = Next;
            Next = 0;
            delete this;
            }
        }
    else
        delete this;
    }

void vertex::adjustWeight()
    {    
    double i = 1.0;
    for (const char * p = this->Pattern.itsTxt(); *p; ++p)
        if (*p == ANY)
            i += 1.0;
    this->wght *= i;
    }

#if _NA
void vertex::adjustNotApplicableCountsByRecalculatingR_NA(trainingPair * pair,int total)
    {
    R__NA = 0;
    while(pair)
        {
        if (pair->notLemmatizedBy(this))
            ++R__NA;
        /*
        switch (lemmatise(pair))
            {
            case wrong:
            case right:
                break;
            default:
                ++R__NA;
            }*/
        pair = pair->next();
        }
    W__NA = total - (R__NA + R__R + R__W + W__R + W__W);
    }

void vertex::adjustNotApplicableCountsByRecalculatingW_NA(trainingPair * pair,int total)
    {
    W__NA = 0;
    while (pair)
        {
        if (pair->notLemmatizedBy(this))
            ++W__NA;
        /*        switch (lemmatise(pair))
            {
            case wrong:
            case right:
                break;
            default:
                ++W__NA;
            }*/
        pair = pair->next();
        }
    R__NA = total - (R__R + R__W + W__NA + W__R + W__W);
    }
#endif

vertex::vertex(shortRulePair * Rule) :
      Head(0)
    , Hash(0)
    , R__R(0)
    , R__W(0)
    , W__R(0)
    , W__W(0)
#if _NA
    , R__NA(0)
    , W__NA(0)
#endif
    , RefCount(0)
    , Relations(0)
    {
    construct((const char*)Rule->itsPatternArray(), (const char*)Rule->itsReplacementArray());
    }

vertex::vertex(const char * pat, const char * rep) :
Head(0), RefCount(0), Relations(0), Hash(0)
, R__R(0)
, R__W(0)
, W__R(0)
, W__W(0)
#if _NA
, R__NA(0)
, W__NA(0)
#endif
    {
    construct(pat, rep);
    }


void vertex::construct(const char * pat, const char * rep)
    {
    char lpattern[1000];
    char lreplacement[1000];
    char * ppat = lpattern;
    char * prep = lreplacement;

    size_t patlen = strlen(pat);
    size_t replen = strlen(rep);
    if (*pat != START)
        {
        if (*rep != START)
            {
            strcpy(ppat, StartAny);
            strcpy(prep, StartAny);
            ppat += 2;
            prep += 2;
            }
        else
            {
            strcpy(ppat, Start);
            ++ppat;
            }
        }
    else if (*rep != START)
        {
        strcpy(prep, Start);
        ++prep;
        }
    strcpy(ppat, pat);
    strcpy(prep, rep);
    ppat += patlen - 1;
    prep += replen - 1;
    if (*ppat == ANY && *prep == ANY)
        {
        strcpy(++ppat, End);
        strcpy(++prep, End);
        }
    else
        {
        if (*ppat != END)
            strcpy(++ppat, AnyEnd);
        if (*prep != END)
            strcpy(++prep, AnyEnd);
        }
    Pattern.dup(lpattern);
    Replacement.dup(lreplacement);
    ++VertexCount;
    }

static matchResult wrongOrFailure(char * p,char *w)
    {
    while (*p)
        {
        if (*p == *w)
            {
            do
                {
                ++p;
                ++w;
                } while (*p && *p != ANY && *p == *w);

            if (*p != ANY && *p != *w)
                {
                return failure;
                }
            }
        else if (*p == ANY)
            {
            ++p;
            char * ep = strchr(p, ANY);
            if (ep)
                {
                *ep = '\0';
                char * sub = strstr(w, p);
                if (sub)
                    {
                    w = sub;
                    }
                else
                    {
                    *ep = ANY;
                    return failure;
                    }
                *ep = ANY;
                }
            else
                {
                char * sub = strstr(w, p);
                if (sub)
                    {
                    w = sub;
                    }
                else
                    {
                    return failure;
                    }
                }
            }
        else
            {
            return failure;
            }
        }
    return wrong;
    }

matchResult vertex::apply(trainingPair * trainingpair)
    {
    static char wrd[100];
    size_t L1 = trainingpair->itsWordlength();
    if (L1 + 3 > sizeof(wrd))
        {
        printf("vertex::apply too small buffer");
        exit(1);
        }
    wrd[0] = START;
    strncpy(wrd + 1, trainingpair->itsWord(), L1);
    wrd[L1 + 1] = END;
    wrd[L1 + 2] = 0;
    char * p = Pattern.itsTxt();
    char * r = Replacement.itsTxt();
    if (!*p)
        {
        static char lStartAnyEnd[4] = { START, ANY, END, 0 };
        p = lStartAnyEnd/*"^*$"*/;
        r = p;
        }
    char * w = wrd;
    const char * lh = trainingpair->itsLemmaHead();
    const char * last = lh + trainingpair->itsLemmalength();
    
    ++r;
    while (*p)
        {
        if (*p == *w)
            {
            while (*r && *r != ANY)
                {
                if (*r != END && (lh == last || *lh++ != *r))
                    {
                    while (*r && *r != ANY)
                        {
                        ++r;
                        }

                    do
                        {
                        ++p;
                        ++w;
                        } while (*p && *p != ANY && *p == *w);

                    if (*p == *r)
                        {
                        return wrongOrFailure(p, w);
                        }
                    else
                        {
                        return failure;
                        }
                    }
                ++r;
                }

            do
                {
                ++p;
                ++w;
                } while (*p && *p != ANY && *p == *w);

            if (*p != *r)
                {
                return failure;
                }
            }
        else if (*p == ANY)
            {
            ++p;
            ++r;
            char * ep = strchr(p, ANY);
            if (ep)
                {
                *ep = '\0';
                char * sub = strstr(w, p);
                if (sub)
                    {
                    while (w < sub)
                        {
                        if (lh == last || *lh++ != *w)
                            {
                            w = sub;
                            *ep = ANY;
                            return wrongOrFailure(p, w);
                            }
                        ++w;
                        }
                    }
                else
                    {
                    *ep = ANY;
                    return failure;
                    }
                *ep = ANY;
                }
            else
                {
                char * sub = strstr(w, p);
                if (sub)
                    {
                    while (w < sub)
                        {
                        if (lh == last || *lh++ != *w)
                            {
                            w = sub;
                            return wrongOrFailure(p, w);
                            }
                        ++w;
                        }
                    }
                else
                    {
                    return failure;
                    }
                }
            }
        else
            {
            return failure;
            }
        }
    if (lh != last)
        return wrong;
    return right;
    }

static int utfchar(char * p, int & U) /* int is big enough for all UTF-8 bytes */
    {
    if (!*p)
        {
        U = 0;
        return 0;
        }
    char * q = p;
    U = *q++;
    while ((*q & 0xC0) == 0x80) /* Assuming valid UTF-8, such bytes can only
                                occur after an initial byte with highest two
                                bits set, so we don't need to check for the
                                existence of that initial byte. */
        {
        U <<= 8;
        U += *q++;
        } 
    return (int)(q - p);
    }

matchResult vertex::applym(trainingPair * pair, char * mask)
    {
    char wrd[100];
    char lemma[100];
    size_t L1 = pair->itsWordlength();
    if (L1 + 3 > sizeof(wrd))
        {
        printf("vertex::apply too small buffer");
        exit(1);
        }
    wrd[0] = START;
    strncpy(wrd + 1, pair->itsWord(), L1);
    wrd[L1 + 1] = END;
    wrd[L1 + 2] = 0;
    char * p = Pattern.itsTxt();
    char * r = Replacement.itsTxt();
    if (!*p) // root
        {
        static char lStartAnyEnd[4] = { START, ANY, END, 0 };
        p = lStartAnyEnd/*"^*$"*/;
        r = p;
        }
    char * w = wrd;
    char * d = lemma;
    char * last = lemma + sizeof(lemma) - 7;
    char * m = mask;
    int P,W,R;
    int pinc = utfchar(p, P);
    int winc = utfchar(w, W);
    int rinc = utfchar(r, R);
    while (P && R)
        {
        if (P == W)
            {
            *m = unequal;
            while (R && R != ANY && d < last)
                {
                strncpy(d, r, rinc);
                r += rinc;
                d += rinc;
                rinc = utfchar(r, R);
                }

            do
                {
                *++m = unequal;
                p += pinc;
                w += winc;
                pinc = utfchar(p, P);
                winc = utfchar(w, W);
                } while (P && P != ANY && P == W);

            if (P != R)
                {
                return failure;
                }
            }
        else if (R == ANY)
            {
            p += pinc;
            r += rinc;
            pinc = utfchar(p, P);
            rinc = utfchar(r, R);
            char * ep = strchr(p, ANY);
            if (ep)
                {
                *ep = '\0';
                char * sub = strstr(w, p);
                if (sub)
                    {
                    while (w < sub && d < last)
                        {
                        *m++ = equal;
                        *d++ = *w++;
                        while ((*w & 0xC0) == 0x80)
                            *d++ = *w++;
                        }
                    winc = utfchar(w, W);
                    }
                else
                    {
                    *ep = ANY;
                    d = lemma;
                    break;
                    }
                *ep = ANY;
                }
            else
                {
                char * sub = strstr(w, p);
                if (sub)
                    {
                    while (w < sub && d < last)
                        {
                        *m++ = equal;
                        *d++ = *w++;
                        while ((*w & 0xC0) == 0x80)
                            *d++ = *w++;
                        }
                    winc = utfchar(w, W);
                    }
                else
                    {
                    d = lemma;
                    break;
                    }
                }
            }
        else
            {
            d = lemma;
            break;
            }
        }
    if (P || R)
        {
        return failure;
        }
    *--m = '\0';
    for (m = mask; *m; ++m)
        m[0] = m[1];
    *d = '\0';
    d = lemma;
    char * oldd = d;
    while ((*++d & 0xC0) == 0x80)
        ;
    while (*d)
        {
        *oldd++ = *d++;
        }
    *oldd = '\0';
    if (oldd > lemma + 1)
        oldd[-1] = '\0';

    if (pair->isCorrect(lemma))
        {
#if AMBIGUOUS
        for(trainingPair * q = pair->Alt;q != pair;q = q->Alt)
            {
            q->set(b_solved);
            q->setRes(notme);                
            }
        pair->set(b_solved);
        pair->setRes(yes);
#endif
        return right;
        }
    else
        {
#if AMBIGUOUS
        switch(pair->getRes())
            {
            case yes:
            case undecided:
            case no:
                {
                for(trainingPair * q = pair->Alt;q != pair;q = q->Alt)
                    {
                    q->setRes(undecided);
                    }
                }
                pair->setRes(no);
                break;
            case notme:
                if(!pair->isset(b_solved))
                    {
                    pair->setRes(no);
                    }
                break;
            }
#endif
        return wrong;
        }
    }

vertexPointer::vertexPointer(vertex * V,vertexPointer * Next,bool InputRight,bool Right)
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

vertexPointer::~vertexPointer()
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

#if PRUNETRAININGPAIRS
bool vertexPointer::fewerLikesThan(int thresh) const
    {
    const vertexPointer * a;
    for(a = this;a;a = a->Next)
        {
        if(a->Right)
            {
            if(a->V->ruleLikes() >= thresh)
                return false;
            }
        }
        
    return true;
    }
#endif