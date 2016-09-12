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

#define _CRT_SECURE_NO_DEPRECATE

#include "node.h"
#include "comp.h"
#include "optionaff.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include "utf8func.h"
#include "isofunc.h"
#include "trainingpair.h"
#include "vertex.h"
#include "hashtable.h"
#include "ruletemplate.h"

//int (*comp)(const vertex * a,const vertex * b) = 0; 
// returns b > a ? 1 : b < a ? -1 : 0
// (Chosen like this to let qsort sort in descending order.)
int visitno;

//FILE * fprune = stdout;

#if AMBIGUOUS
int ambivalentWords;
int alternatives;
#endif
int allwords;

int NodeCount = 0;

int RulePairCount = 0;
int StrngCount = 0;

int VCount = 0;

int node::mostPenalized = 1;

bool building = true;


        /*
    & ( translatePat
      =   p,n,f
        .     !arg:(?n.?arg)
            & @(!arg:?p "*" ?arg)
            & chr$!n:?f
            &   !p
                ( glf$('(?.$f)):(=?f)
                & !f
                )
                translatePat$(!n+1.!arg)
          | !arg:?RR&
      )
      */


#if BRACMATOUTPUT
strng * protect(strng * p)
    {
    const char * txt = p->itsTxt();
    if(strcspn(txt,"()\"!?' \t{}[]-#@%/\\;+*^$`&|=_:,.><") < strlen(txt))
        {
        strng * ret = new strng("\"");
        strng slash("\\\\");
        strng quot("\\\"");
        strng qu("\"");
        const char * pq;
        while((pq = strpbrk(txt,"\\\"")) != NULL)
            {
            strng stxt(txt,(size_t)(pq - txt));
            if(*pq == '\\')
                ret->cat(&stxt,&slash,(const strng *)0);
            else
                ret->cat(&stxt,&quot,(const strng *)0);
            txt = pq+1;
            }
        strng stxt2(txt);
        ret->cat(&stxt2,&qu,(const strng *)0);
        return ret;
        }
    return NULL;
    }
#endif

#if BRACMATOUTPUT
static const strng * translatePat(int f,const strng * arg,strng ** RR)
#else
static void translatePat(const strng * arg,strng ** RR)
#endif
    {
    ptrdiff_t star = arg->pos(ANY);
#if BRACMATOUTPUT
    strng * ret;
#endif
    if(star >= 0)
        {
        strng * npat = arg->substr(star+1);
#if BRACMATOUTPUT
        strng * P = arg->substr(0,star);
        strng * p = protect(P);if(p) delete P; else p = P;
        const strng * rem = translatePat(f+1,npat,RR);
        strng blank(" ");
        strng question("?");
        strng var(f);
        ret = p;
        ret->cat(&blank,&question,&var,&blank,rem,(const strng *)0);
        delete rem;
#else
        translatePat(npat,RR);
#endif
        delete npat;
        }
    else
        {
        if(RR)
            {
            delete *RR;
            *RR = new strng(arg);
            }
#if BRACMATOUTPUT
        ret = new strng("");
#endif
        }
#if BRACMATOUTPUT
    return ret;
#endif
    }

#if BRACMATOUTPUT
static const strng * translatePat2(int f,const strng * arg,strng ** patreps,int index)
#else
static void translatePat2(const strng * arg,strng ** patreps,int index)
#endif
    {
    ptrdiff_t star = arg->pos(ANY);
#if BRACMATOUTPUT
    strng * ret;
#endif
    if(star >= 0)
        {
        delete patreps[index];
        patreps[index] = arg->substr(0,star);
        strng * npat = arg->substr(star+1);
#if BRACMATOUTPUT
        strng * P = arg->substr(0,star);
        strng * p = protect(P);if(p) delete P; else p = P;
        const strng * rem = translatePat2(f+1,npat,patreps,index+2);
        strng blank(" ");
        strng question("?");
        strng var(f);
        ret = p;
        ret->cat(&blank,&question,&var,&blank,rem,(const strng *)0);
        delete rem;
#else
        translatePat2(npat,patreps,index+2);
#endif
        delete npat;
        }
    else
        {
        delete patreps[index];
        patreps[index] = new strng(arg);
#if BRACMATOUTPUT
        ret = new strng("");
#endif
        }
#if BRACMATOUTPUT
    return ret;
#endif
    }
        /*
    & ( translateRep
      =   p,n,f,rem
        .     !arg:(?n.?arg)
            & @(!arg:?p "*" ?arg)
            & translateRep$(!n+1.!arg):(=?rem)
            & chr$!n:?f
            & glf$('(!.$f)):(=?f)
            & (   !p:~
                & (   '$rem:(=)
                    & '($p ()$f)
                  | '($p ()$f ()$rem)
                  )
              | '$rem:(=)&'$f
              | '($f ()$rem)
              )
          | '$arg
      )
      */
#if BRACMATOUTPUT
static const strng * translateRep(int f,const strng * arg,strng ** patreps,int index)
#else
static void translateRep(const strng * arg,strng ** patreps,int index)
#endif
    {
    ptrdiff_t star = arg->pos(ANY);
#if BRACMATOUTPUT
    strng * ret;
#endif
    if(star >= 0)
        {
        patreps[index] = arg->substr(0,star);
        strng * narg = arg->substr(star+1);
#if BRACMATOUTPUT
        strng * P = arg->substr(0,star);
        strng * p = protect(P);if(p) delete P; else p = P;
        const strng * rem = translateRep(f+1,narg,patreps,index+2);
        strng blank(" ");
        strng bang("!");
        strng var(f);
        ret = p;
        if(star)
            {
            if(rem->length() == 0)
                ret->cat(&blank,&bang,&var,(const strng *)0);
            else
                ret->cat(&blank,&bang,&var,&blank,rem,(const strng *)0);
            }
        else if(rem->length() == 0)
            {
            ret->cat(&bang,&var,(const strng *)0);
            }
        else
            {
            ret->cat(&bang,&var,&blank,rem,(const strng *)0);
            }
        delete rem;
#else
        translateRep(narg,patreps,index+2);
#endif
        delete narg;
        }
    else
        {
        delete patreps[index];
        patreps[index] = new strng(arg);
#if BRACMATOUTPUT
        strng * Ret = new strng(arg);
        ret = protect(Ret);if(ret) delete Ret; else ret = Ret;
#endif
        }
#if BRACMATOUTPUT
    return ret;
#endif
    }
        /*
    & ( makeNode
      =   L R pat rep
        .   !arg:((?pat.?rep).?L.?R)
          & @(!pat:!L ?pat !R)
          & ( @(!pat:?LL ("*" ?:?pat))
            | :?LL
            )
          & :?RR
          & translatePat$(asc$A.!pat):?pat
          & ( (     !LL:~
                  & (   !RR:~
                      &   
                        ' ( $LL ?W ()$RR
                          & @(!W:$pat)
                          )
                    | '($LL ?W&@(!W:$pat))
                    )
                |   !RR:~
                  & '(?W ()$RR&@(!W:$pat))
                | '(?W:$pat)
              . translateRep$(asc$A.!rep)
              )
            . !nr+1:?nr
            )
      )
    & (makeNodeAndTrim=.makeNode$!arg.!LL.!RR)

*/

#if BRACMATOUTPUT
strng * makeNode(strng ** patreps,int & nr,const strng * pat,const strng * rep,const strng * L,const strng * R,strng ** pLL,strng ** pRR)
#else
void makeNode(strng ** patreps,const strng * pat,const strng * rep,const strng * L,const strng * R,strng ** pLL,strng ** pRR)
#endif
    {
    strng * LL = 0;
    strng * RR = 0;
    const strng * npat = pat->salad(L,R);
    ptrdiff_t star = npat->pos(ANY);
    delete patreps[0];
    if(star >= 0)
        {
        LL = npat->substr(0,star);
#if BRACMATOUTPUT
        const strng * dummy = translatePat2('A',npat,patreps,0);
        delete dummy;
#else
        translatePat2(npat,patreps,0);
#endif
        const strng * nnpat = npat->substr(star);
        delete npat;
        npat = nnpat;
        }
    else
        {
        LL = new strng("");
        patreps[0] = new strng("");
#if BRACMATOUTPUT
        const strng * dummy = translatePat2('A',npat,patreps,0);
        delete dummy;
#else
        translatePat2(npat,patreps,0);
#endif
        }
#if BRACMATOUTPUT
    const strng * nnnpat = translatePat('A',npat,&RR);
#else
    translatePat(npat,&RR);
#endif
    delete npat;
    LL->trim();
    RR->trim();
#if BRACMATOUTPUT
    strng blank(" ");
    strng question("?");
    strng bang("!");
    strng questionW("?W");
    strng * ret = new strng("(((=");

    strng part2("&@(!W:");
    strng part3(").(=");
    strng part2bis("?W:");
    strng rpar(")");

    if(LL->length() > 0)
        {
        strng * ll = protect(LL);if(!ll) ll = LL;
        if(RR->length() > 0)
            {
            strng * rr = protect(RR);if(!rr) rr = RR;
            ret->cat(ll,&blank,&questionW,&blank,rr,&part2,nnnpat,&rpar,(const strng *)0);
            if(rr != RR)
                delete rr;
            }
        else
            {
            ret->cat(ll,&blank,&questionW,&part2,nnnpat,&rpar,(const strng *)0);
            }
        if(ll != LL)
            delete ll;
        }
    else if(RR->length() > 0)
        {
        strng * rr = protect(RR);if(!rr) rr = RR;
        ret->cat(&questionW,&blank,rr,&part2,nnnpat,&rpar,(const strng *)0);
        if(rr != RR)
            delete rr;
        }
    else
        {
        ret->cat(&part2bis,nnnpat,(const strng *)0);
        }
    delete nnnpat;
    const strng * nrep = translateRep('A',rep,patreps,1);
    ++nr;
    char buf[12];
    sprintf(buf,"%d",nr);
    strng Nr(buf);
    strng pardot(")).");
//    strng pardot(".");
    ret->cat(&part3,nrep,&pardot,&Nr,&rpar,(const strng *)0);
    delete nrep;
#else
    translateRep(rep,patreps,1);
#endif
    if(pLL)
        *pLL = LL;
    else
        delete LL;

    if(pRR)
        *pRR = RR;
    else
        delete RR;
#if BRACMATOUTPUT
    return ret;
#endif
    }

char * dup(const char * buf)
    {
    char * Txt = new char[strlen(buf)+1];
    strcpy(Txt,buf);
    return Txt;
    }

char * dupn(const char * buf,size_t n)
    {
    char * Txt = new char[n+1];
    strncpy(Txt,buf,n);
    Txt[n] = 0;
    return Txt;
    }



#if WRITEINVERTED
char * strnrev(const char * s,size_t n)
    {
    if(n == 0)
        n = strlen(s);
    static char * ret = 0;
    static size_t N = 0;
    if(N < n+1)
        {
        if(N > 0)
            delete [] ret;
        N = n+1;
        ret = new char[N];
        }
    size_t i = n;
    char * d = ret;
    if(UTF8)
        {
        while(i > 0)
            {
            --i;
            if(s[i] & 0x80)
                {
                while(i > 0 && (s[--i] & 0xC0) != 0xC0)
                    ;
                size_t ii = i;
                do 
                    *d++ = s[i];
                while((s[++i] & 0xC0) != 0xC0 && s[i] & 0x80);
                i = ii;
                }
            else
                *d++ = s[i];
            }
        }
    else
        while(i > 0)
            *d++ = s[--i];
    *d = '\0';
    return ret;
    }

char * Strrev(char * s)
    {
    return strnrev(s,strlen(s));
    }
#endif


void node::Counting(int & nodes,int & pairs,FILE * f)
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

double node::weightedcount()
    {
    double ret = 0.0;
    node * n = this;
    while(n)
        {
        if(n->Right)
            {
            double rcount = (double)n->Right->count();
            assert(rcount >= 0.9);
            /* 
            ret += 5.0*(rcount-1.0)*exp(-0.9*rcount)+1; 
            penalizes as follows:
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
            /*
            Purpose: for rule sets that have rules with support
            from fewer than 3 word/lemma examples removed.
            Less pronounced difference between values for max 
            penalty at 
            rcount == 3 (1.1036) and for 
            rcount == 2 (1.0268) and
            rcount == 1 (0.7165).

            Weight goes to zero for large rcount.
            rcount == 30 (0.00136)
            */
            ret += rcount*rcount*exp(-2.0*rcount/mostPenalized);
            /* Squared version of previous function. Gives more 
               pronounced penalty at cutoff+1, while penalties
               for rules with fewer examples are just barely
               enough penalized to not let cutoff+1 rules split
               in rules that are later cut away. In other words:
               tends to create rules with more, not less examples
               than cutoff+1.
            */
            }
        if(n->IfPatternSucceeds)
            ret += n->IfPatternSucceeds->weightedcount();
        n = n->IfPatternFails;
        }
    return ret;
    }

double node::entropy(double Nnodes)
    {
    double ret = 0.0;
    node * n = this;
    while (n)
        {
        if (n->Right)
            {
            double rcount = (double)n->Right->count();
            rcount /= Nnodes;
//                        assert(rcount >= 0.9);
            ret += rcount*log(rcount);
            }
        if (n->IfPatternSucceeds)
            ret += n->IfPatternSucceeds->entropy(Nnodes);
        n = n->IfPatternFails;
        }
    return ret;
    }

int node::prune(int threshold)
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


void node::printSep(FILE * f,int level)
    {
    this->V->print(f,level);
    if(Right)
        Right->printSep(f);
    if(this->IfPatternSucceeds)
        this->IfPatternSucceeds->printSep(f,level+1);
    if(this->IfPatternFails)
        this->IfPatternFails->printSep(f,level);
    }

int node::print(FILE * fo,int ind,int & Nnodes,int &NnodesR)
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

node::node(vertex * V):V(V),IfPatternSucceeds(0),IfPatternFails(0),Right(0)
    {
    ++NodeCount;
    V->incRefCnt();
    }

node::~node()
    {
    --NodeCount;
    delete IfPatternSucceeds;
    delete IfPatternFails;
    V->destroy();
    }


LONG node::countBySize()
    {
    LONG ret = 0;
    for ( node * n = this
        ; n
        ; n = n->IfPatternFails
        )
        {
        if(n->Right)
            {
            ret += V->Pattern.length();
            }
        
        if(n->IfPatternSucceeds)
            ret += n->IfPatternSucceeds->countBySize();
        }
    return ret;
    }

void node::splitTrainingPairList(trainingPair * pair,trainingPair **& pNotApplicable,trainingPair **& pWrong,trainingPair **& pRight,bool SuffixOnly)
    {
    while(pair)
        {
        trainingPair * nxt = pair->next();
        pair->setNext(0);
        static char mask[100];
        matchResult res = V->applym(pair,mask);
        if(res == failure)
            {
            *pNotApplicable = pair;
            pNotApplicable = pair->pNext();
            }
        else
            {
            if(SuffixOnly)
                suffix(mask);
            pair->takeMask(mask);
            if(res == wrong)
                {
                *pWrong = pair;
                pWrong = pair->pNext();
                }
            else
                {
                *pRight = pair;
                pRight = pair->pNext();
                }
            }
        pair = nxt;
        }
    }

bool node::compatibleSibling(node * sib)
    {
    if(this->V->Pattern.dif(&sib->V->Pattern) != dif_incompatible)
        return true;
    else if(IfPatternFails)
        return IfPatternFails->compatibleSibling(sib);
    else
        return false;
    }

node * node::cleanup(node * parent)
    {
    if(this->IfPatternSucceeds) 
        {
        this->IfPatternSucceeds = this->IfPatternSucceeds->cleanup(this);
        }
    if(IfPatternFails)
        {
        IfPatternFails = IfPatternFails->cleanup(parent);
        }
    if(  IfPatternFails 
      && (  this->IfPatternSucceeds 
         || IfPatternFails->compatibleSibling(this)
         )
      )
        {
        return this;
        }
    else if(parent)
        {
        if(this->Right)
            {
            if(!parent->IfPatternFails && !parent->Right)
                {
                // remove parent, keep this.
                return this;
                }
            trainingPair * R = this->Right;
            for(;;)
                {
                matchResult res = parent->V->lemmatise(R);

                if(res != right)
                    {
                    return this;
                    }
                if(R->next())
                    R = R->next();
                else
                    break;
                }
            R->setNext(parent->Right);
            parent->Right = this->Right;
            this->Right = 0;
            }
        node * ret = IfPatternFails;
        if(ret)
            {
            IfPatternFails = 0;
            }
        else
            {
            ret = this->IfPatternSucceeds;
            this->IfPatternSucceeds = 0;
            }
        delete this;
        return ret;
        }
    else
        {
        return this;
        }
    }

#if PRUNETRAININGPAIRS
static int deletedRules = 0;
static trainingPair * pruneTrainingPairs(trainingPair * pairs,trainingPair *& deletedPairs)
    {            
    trainingPair * tp = pairs;
    while(tp && tp->fewerLikesThan(2))
        {
        ++deletedRules;
        trainingPair * nxt = tp->next();
        //delete tp; Impossible. Training pairs are allocated as an array.
        tp->Next = deletedPairs;
        deletedPairs = tp;
        tp = nxt;
        }
    trainingPair * ret = tp;
    if(tp)
        {
        trainingPair * lastGood = tp;
        tp = tp->next();
        lastGood->setNext(0);
        while(tp)
            {
            if(tp->fewerLikesThan(2))
                {
                ++deletedRules;
                trainingPair * nxt = tp->next();
                //delete tp; Impossible. Training pairs are allocated as an array.
                tp->Next = deletedPairs;
                deletedPairs = tp;
                tp = nxt;
                }
            else
                {
                lastGood->setNext(tp);
                lastGood = tp;
                tp = tp->next();
                lastGood->setNext(0);
                }
            }
        }
    return ret;
    }

vertex ** cleanUpUnusedVertices(vertex ** pvf, vertex ** pvN
#if _NA
                                      , trainingPair * deletedRightPairs, trainingPair * deletedWrongPairs
#endif
                                      )
    {
#if _NA
#endif
    for(vertex ** pvi = pvf;pvi < pvN;)
        {
        vertex * V = *pvi;
        if(  V->R__R == 0 
          && V->W__R == 0 
          && V->R__W == 0
          && V->W__W == 0
          )
            { // rule has become irrelevant
            V->destroy();
            --pvN;
            *pvi = *pvN;
            }
        else
            {
#if _NA
            int M;
            M = deletedRightPairs->notLemmatizedBy(V);
            V->R__NA -= M;
            M = deletedWrongPairs->notLemmatizedBy(V);
            V->W__NA -= M;
#endif
            ++pvi;
            }
        }
    return pvN;
    }
#endif

int skipping(unsigned long N,trainingPair * pairs)
    {
    if(pairs)
        {
//#define PROD 3ULL
//#define MINTRAINING 2
#define PROD 10000000ULL
#define MINTRAINING 1000
        unsigned long M = pairs->count();
        unsigned long long prod = (unsigned long long)M*(unsigned long long)N;
        if(prod > PROD)
            {
            unsigned long long m = PROD / (unsigned long long)N;
            if(m < MINTRAINING)
                m = MINTRAINING;
            int skip = (int)((M / m) - 1);
            return skip;
            }
        }
    return 0;
    }

void node::printNumbers(trainingPair ** allRight,trainingPair ** allWrong)
    {
    int r = 0;
    int w = 0;
    for (trainingPair * tp = *allRight; tp; tp = tp->next())
        {
        ++r;
        }
    for (trainingPair * tp = *allWrong; tp; tp = tp->next())
        {
        ++w;
        }
    printf("node(right=%d,wrong=%d) rule:",r,w);
    this->V->print1(stdout);
    printf("\n");
    }

void node::unsetSolved(trainingPair ** allPairs)
    {
    for(trainingPair * tp = *allPairs;tp;tp = tp->next())
        {
        tp->unset(b_solved);
        }
    }

void ifThereAreUnmatchedWordsThenTellUsAndExit(vertex ** pv,trainingPair * Wrong, unsigned long first, unsigned long lastN)
    {
    fprintf(stderr,"***** first:%lu >= lastN:%lu\n",first,lastN);
    fprintf(stderr,"***** (This happens if one or more wrongly lemmatised words aren't handled by any candidate rules.)\n");
    fprintf(stderr,"***** (List of candidate rules:)\n");
    fprintf(stderr,"***** pv:\n");
    for(unsigned long i = 0;i < lastN;++i)
        {
        fprintf(stderr,"i=%lu:",i);
        pv[i]->print1(stderr);
        fprintf(stderr,"\n");
        }
    fprintf(stderr,"***** (List of unmatched words that need better lemmatisation rule(s):)\n");
    fprintf(stderr,"\n***** Wrong:\n");
    Wrong->printAll(stderr,"unmatched words that need better lemmatisation rule(s)\n",'\n');
    printf("\nAFFIXTRAIN failed\n");
    fprintf(stderr,"\n***** AFFIXTRAIN failed\n");
    exit(-1);
    }

void letTheBestRuleBeFirst(vertex ** pvf, vertex ** pvN,optionStruct * options)
    {
    if (comp == comp_parms)
        {
        for (vertex ** pvi = pvf; pvi < pvN; ++pvi)
            {
            computeWeight(*pvi);
            }
        if (options->verbose())
            {
            printf("Weights computed     %c",STARTLINE);
            }
        if (options->expensiveInfix())
            {
            for (vertex ** pvi = pvf; pvi < pvN; ++pvi)
                {
                (*pvi)->adjustWeight();
                }
            if (options->verbose())
                {
                printf("Weights adjusted     %c",STARTLINE);
                }
            }
        for (vertex ** pvi = pvf + 1; pvi < pvN; ++pvi)
            {
            if (comp_parms(*pvf, *pvi) < 0)
                {
                vertex * tmp = *pvf;
                *pvf = *pvi;
                *pvi = tmp;
                }
            }
        if (options->verbose())
            {
            printf("candidates sorted     %c",STARTLINE);
            }
        }
    else
        {
        for (vertex ** pvi = pvf + 1; pvi < pvN; ++pvi)
            {
            if (comp(*pvf, *pvi) < 0)
                {
                vertex * tmp = *pvf;
                *pvf = *pvi;
                *pvi = tmp;
                }
            }
        }
    }


vertex ** adjustCounts(vertex ** pvf, vertex ** pvN,optionStruct * options)
    {
    if (options->verbose())
        {
        printf("delete the rules that have become irrelevant and adjust counts for the good ones.    %c",STARTLINE);
        }
    for (vertex ** pvi = pvf; pvi < pvN; )
        {
        vertex * V = *pvi;
        if(  V->R__R == 0 
          && V->W__R == 0 
          && V->R__W == 0
          && V->W__W == 0
          )
            { /* The rule has become irrelevant. Delete it and put the last in 
              the list in its place. (The list is unsorted, so that is safe to
              do.) The list becomes a shorter. */
            V->destroy();
            --pvN;
            *pvi = *pvN;
            /* Do not increment pvi! */
            }
        else
            {
#if _NA
            if(outputR < outputW)
                {
                V->adjustNotApplicableCountsByRecalculatingR_NA(this->Right, outputR + outputW);
                }
            else
                {
                V->adjustNotApplicableCountsByRecalculatingW_NA(Wrong, outputR + outputW);
                }
#endif
            ++pvi;
            }
        }
    if (options->verbose())
        {
        printf("delete the rules that have become irrelevant and adjust counts for the good ones. DONE                                            %c",STARTLINE);
        }
    return pvN;
    }

#if PRUNETRAININGPAIRS
unsigned long pruneTheWeirdWords(trainingPair *& Right, trainingPair *& Wrong, vertex ** pv, unsigned long N, optionStruct * options)
    {
    trainingPair * deletedRightPairs = 0;
    if (options->verbose())
        {
        printf("pruning right pairs %c",STARTLINE);
        }
    Right = pruneTrainingPairs(Right, deletedRightPairs);
    trainingPair * deletedWrongPairs = 0;
    if (options->verbose())
        {
        printf("pruning wrong pairs %c",STARTLINE);
        }
    Wrong = pruneTrainingPairs(Wrong, deletedWrongPairs);
    if (options->verbose())
        {
        printf("clean up unused candidate rules               <-- %lu%c",N,STARTLINE);
        }
    N = (unsigned long)(cleanUpUnusedVertices(pv, pv + N
#if _NA
        , deletedRightPairs, deletedWrongPairs
#endif
        ) - pv);
    if (options->verbose())
        {
        printf("clean up unused candidate rules  %lu%c", N,STARTLINE);
        }
    if (deletedRightPairs)
        deletedRightPairs->mourn();
    if(deletedWrongPairs)
        deletedWrongPairs->mourn();
    if (options->verbose())
        {
        printf("clean up pruned training pairs%c",STARTLINE);
        }
    return N;
    }
#endif


void node::init(trainingPair ** allRight,trainingPair ** allWrong,int level,optionStruct * options)
    {
    /* At entry, this node has a rule (a pattern and a replacement).
       The parent's two sets of training pairs are given as arguments, so that
       this node, this node's siblings or this node's children may "steal" 
       them.
       The return value are the leftovers, the trainingpairs for which no
       rules could be made. (This can only be the case if the parent's rule's
       pattern matches all of a word. The leftovers consist of at most one
       training pair.)
    */
    if (options->verbose())
        printNumbers(allRight,allWrong);

    trainingPair * NotApplicableRight = NULL;
    trainingPair * NotApplicableWrong = NULL;
    trainingPair * Wrong = NULL;

    trainingPair ** pNotApplicableRight = &NotApplicableRight; // Pattern fails.
    trainingPair ** pNotApplicableWrong = &NotApplicableWrong; // Pattern fails.
    trainingPair ** pWrong = &Wrong; // Pattern succeeds and replacement is wrong.
    trainingPair ** pRight = &this->Right; // Pattern succeeds and replacement is right.
    assert(this->Right == NULL);

    unsetSolved(allRight);
    unsetSolved(allWrong);

    this->splitTrainingPairList(*allRight,pNotApplicableRight,pWrong,pRight,options->suffixOnly());
    this->splitTrainingPairList(*allWrong,pNotApplicableWrong,pWrong,pRight,options->suffixOnly());
    *pNotApplicableRight = 0;
    *pNotApplicableWrong = 0;
    *pWrong = 0;
    *pRight = 0; // Pattern succeeds and replacement is right.
    *allRight = NotApplicableRight;
    *allWrong = NotApplicableWrong;

    this->Right->allDeleteRules();
#if AMBIGUOUS
    trainingPair * Ambiguous = 0;
#endif

    if(Wrong)
        {
        Wrong->allDeleteRules();
#if AMBIGUOUS
        Wrong = Wrong->makeWrongAmbiguousIfRightPresent(Ambiguous);
#endif
        }

    if(Wrong)
        {
        /* Create candidate rules for all inputs and make pointers to each of
           the rules from a list. */
        unsigned long N;
        hashTable VertexHash(1000);
        Wrong->makeCandidateRules(&VertexHash,this->V,false,options);
        if(this->Right)
            this->Right->makeCandidateRules(&VertexHash,this->V,true,options);
        N = VertexHash.recordCount();
        vertex ** pv = VertexHash.convertToList();
        if(N == 0)
            {
            delete [] pv;
            return;
            }
        assert(N > 0);
        if(options->verbose())
            {
            printf("             %lu candidate rules %c",N,STARTLINE);
            }

        /* Test all/some pairs (up to upper bound) on all candidate rules. */
        int skipW = skipping(N,Wrong);
        int skipR = skipping(N,this->Right);
        for(unsigned long i = 0;i < N;++i)
            {
            if(options->verbose())
                {
                printf("%lu %c",i,STARTLINE);
                }
            /* Reset all counters in all candidate rules:
                    R__W = R__R = W__W = W__R = 0 */
            pv[i]->nlemmatiseStart();

            pv[i]->nlemmatise(Wrong,false,skipW);
            pv[i]->nlemmatise(this->Right,true,skipR);
            }
#if PRUNETRAININGPAIRS
        N = pruneTheWeirdWords(this->Right, Wrong, pv, N, options);
#endif
        /* Point to the pointer to the head of the sequence of children. */
        node ** pnode = &this->IfPatternSucceeds;
        unsigned long first = 0; /* Index of the first in the list of
                                 candidates. Incremented for each candidate
                                 that made it to a new node. */
        unsigned long lastN = N; /* Index of the last in the list of
                                 candidates. Decremented as candidates become
                                 irrelevant. */
        /* Wrong is the list of training pairs that are not correctly
           lemmatised by this rule. We make them disappear by creating
           child nodes that absorb them. */
        while(Wrong)
            {
            if (options->verbose())
                {
                printf("Finding the best candidate child rule%c",STARTLINE);
                }
            /* Rules having position on list beyond lastN do not apply to any
               of the remaining pairs, except perhaps if pairs were skipped for
               saving processing time. */
            if (first >= lastN)
                {
                ifThereAreUnmatchedWordsThenTellUsAndExit(pv, Wrong, first, lastN);
                }
            assert(first < lastN);
            vertex ** pvf = pv+first;
            vertex ** pvN = pv+lastN;
            letTheBestRuleBeFirst(pvf, pvN, options);
            *pnode = new node(*pvf++);
            (*pnode)->init(&this->Right,&Wrong,level+1,options);
#if _NA
            int outputR = (this->Right ? this->Right->count() : 0); 
            int outputW = (Wrong ? Wrong->count() : 0);
#endif
            if(skipW == 0 && skipR == 0)
                pvN = adjustCounts(pvf, pvN, options);

            lastN = (unsigned long)(pvN - pv);

            ++first;
            if (options->verbose()) printf("Treat the unmatched inputs with this node's next sibling. \n");
            pnode = &(*pnode)->IfPatternFails;
            } // while(Wrong)

        if (options->verbose())
            {
            printf("All input is handled as it should. Clean up before returning \n");
            }

        *pnode = 0;
        for(unsigned long i = 0;i < lastN;++i)
            {
            pv[i]->destroy();
            }
        delete [] pv;
        if (options->verbose())
            {
            printf("Cleaned up. Return. \n");
            }
        }
    }

