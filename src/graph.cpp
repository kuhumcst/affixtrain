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

#include "graph.h"
#include "comp.h"
#include "optionaff.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include "utf8func.h"
#include "isofunc.h"

//int (*comp)(const vertex * a,const vertex * b) = 0; 
// returns b > a ? 1 : b < a ? -1 : 0
// (Chosen like this to let qsort sort in descending order.)
int visitno;

#if AMBIGUOUS
int ambivalentWords;
int alternatives;
#endif
int allwords;

int NodeCount = 0;

int VertexPointerCount = 0;
int VertexCount = 0;
int TrainingPairCount = 0;
int HashCount = 0;
int RulePairCount = 0;
int StrngCount = 0;
int RuleTemplateCount = 0;
int ShortRulePairCount = 0;

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
            strng stxt(txt,pq - txt);
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

bool strongerCondition(char ** A,char ** B)
    {
    // Return true if all B are contained in an A.
    // The loop exits when there are no more A's.
    // If there still is a B, then this B must be unmatched, and so the
    // function must return false.
    // If there is no B left, it may still be the case that the last
    // B is not contained in an A, so we also have to check for the
    // result from the last call to strstr.
    char * sub = 0;
    char * a = *A;
    char * b = *B;
    while(b && a)
        {
        sub = 0;
        while(a)
            {
            sub = strstr(a,b);
            if(sub)
                {
                a = sub + strlen(b);
                if(!*a)
                    a = *++A;
                break;
                }
            a = *++A;
            }
        b = *++B;
        }
    return !*B && sub;
    }

void cut(char * s,char ** a)
    {
    *a = s;
    while(*s && *s != ANY)
        ++s;
    if(*s)
        {
        *s++ = '\0';
        cut(s,++a);
        }
    else
        *++a = 0;
    }

void uncut(char ** a)
    {
    if(*++a)
        {
        (*a)[-1]= ANY;
        uncut(a);
        }
    }

bool compatibleCondition(char ** A,char ** B)
    {
    char * a = *A;
    char * b = *B;
    // Check compatibility at start of pattern
    if(a[0] == START && b[0] == START)
        {
        int i;
        for(i = 1;a[i] && b[i];++i)
            if(a[i] != b[i])
                return false;
        }
    while(A[1])
        ++A;
    while(B[1])
        ++B;
    if(a != *A || b != *B)
        {
        a = *A;
        b = *B;

        // Check compatibility at end of pattern
        size_t alast = strlen(a) - 1;
        size_t blast = strlen(b) - 1;
        if(a[alast] == END && b[blast] == END)
            {
            int i,j;
            for(i = (int)alast - 1,j = (int)blast - 1;i >= 0 && j >= 0;--i,--j)
                if(a[i] != b[j])
                    return false;
            }
        }
    return true;
    }


edif dif(char * Txt, char * s_Txt)
    {
    /* returns 1 if Txt is a stronger condition (further from the root) 
       than s_Txt.
       It returns -1 if Txt is weaker than s_Txt
       It retuns 0 if the conditions are incommensurable.
    */
    char * A[100];
    char * B[100];
    cut(Txt,A);
    cut(s_Txt,B);
    if(compatibleCondition(A,B))
        {
        bool AB = strongerCondition(A,B);
        bool BA = strongerCondition(B,A);

        edif res = (AB && !BA) ? dif_bigger : (BA && !AB) ? dif_smaller : (AB && BA) ? dif_equal : dif_incommensurable;
        uncut(A);
        uncut(B);
        assert(res == dif_incommensurable || res == dif_smaller || res == dif_bigger || res == dif_equal);
        return res;
        }
    else
        {
        uncut(A);
        uncut(B);
        return dif_incompatible;
        }
    }


edif vertex::dif(vertex * other)
    {
    return ::dif(Pattern.itsTxt(),other->Pattern.itsTxt());
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

strng::strng(const char * buf)
    {
    Txt = ::dup(buf);
    ++StrngCount;
    }

strng::strng(const char * buf,size_t length)
    {
    Txt = dupn(buf,length);
    ++StrngCount;
    }




edif strng::dif(strng * s) // returns Txt - s
    {
    return ::dif(Txt,s->Txt);
    }

void trainingPair::addRule(vertex * V,bool InputRight,bool Right)
    {
    this->applicableRules = new vertexPointer(V,this->applicableRules,InputRight,Right);
    }

void trainingPair::deleteRules()
    {
    if(this->applicableRules)
        {
        this->applicableRules->deleteAll();
        this->applicableRules = 0;
        }
    }

void trainingPair::allDeleteRules()
    {
    trainingPair * x = this;
    while(x)
        {
        trainingPair * n = x->Next;
        x->deleteRules();
        x = n;
        }
    }


trainingPair::~trainingPair()
    {
    deleteRules();
    delete [] Lemma;
    Lemma = NULL;
    delete [] Mask;
    Mask = NULL;
    --TrainingPairCount;
    }

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

void vertex::print1(FILE * f)
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

vertex::vertex(vertex * Rule,hash * Hash):
        Head(0),RefCount(0),Relations(0),Hash(Hash)
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

bool trainingPair::isCorrect(const char * lemma) const
    {
    const char * a = lemma;
    const char * b = LemmaHead;
    const char * e = LemmaHead + lemmalength;
    while(b < e && *a == *b)
        {
        ++a;
        ++b;
        }
    return !*a && b == e;
    }

int trainingPair::printAll(FILE * f,const char * h,int s)
    {
    trainingPair * p = this;
    int n = 0;
    if(p)
        {
        fprintf(f,"%s",h);
        while(p)
            {
            ++n;
            p->print(f);
            if(p->Next)
                fprintf(f,"%c",s);
            p = p->Next;
            }
        fprintf(f,"\n");
        }
    return n;
    }

void trainingPair::print(FILE * f)
    {
#if AMBIGUOUS
    fprintf(f,"%.*s\t%.*s\t%s ",(int)(wordlength),Word,(int)(lemmalength),LemmaHead,getRes() == undecided ? "undecided" : getRes() == yes ? "yes" : getRes() == no ? "no" : getRes() == notme ? "notme" : "???" );
#else
    fprintf(f,"%.*s\t%.*s\t%s ",(int)(wordlength),Word,(int)(lemmalength),LemmaHead,getRes() == undecided ? "undecided" : getRes() == yes ? "yes" : getRes() == no ? "no" : "???" );
#endif
    }

void trainingPair::printMore(FILE * f)
    {
    if(!f)
        return;
    fprintf(f,"%.*s\t%.*s",(int)(wordlength),Word,(int)(lemmalength),LemmaHead);
#if WORDCLASS
    fprintf(f,"\t%.*s",(int)(wordclasslength),WordClass);
#endif
#if LEMMACLASS
    fprintf(f,"\t%.*s",(int)(lemmaclasslength),LemmaClass);
#endif
    fprintf(f,"\n");
    }

void trainingPair::printSep(FILE * f)
    {
    trainingPair * tp = this;
    do
        {
        ++allwords;
        fprintf(f,"%.*s\t%.*s\t\n",(int)(tp->wordlength),tp->Word,(int)(tp->lemmalength),tp->LemmaHead);
#if AMBIGUOUS
        if(tp->Alt && tp->Alt != tp)
            {
            ++ambivalentWords;
            for(trainingPair * p = tp->Alt;p != tp;p = p->Alt)
                {
                ++alternatives;
//                fprintf(f,"%.*s\t%.*s\t@@%s\n",(int)(p->wordlength),p->Word,(int)(p->lemmalength),p->LemmaHead,p->isset(b_test) ? "t" : "");
                fprintf(f,"%.*s\t%.*s\n",(int)(p->wordlength),p->Word,(int)(p->lemmalength),p->LemmaHead);
                }
            }
#endif
        tp = tp->Next;
        }
    while(tp);
    }

matchResult vertex::lemmatisem(trainingPair * pair, char ** pmask, char ** plemma, optionStruct * options)
    {
    static char lemma[100];
    static char mask[100];
    if (applym(pair, sizeof(lemma), lemma, mask, options))
        {
        if (plemma)
            *plemma = lemma;
        if (pmask)
            *pmask = mask;
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
    else
        {
        return failure;
        }
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

trainingPair * trainingPair::nth(int n)
    {
    trainingPair * pairs = this;
    while(pairs && n > 0)
        {
        pairs = pairs->next();
        --n;
        }
    return pairs;
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

int vertex::nlemmatise(trainingPair * pair,int n,bool InputRight)
    {
    int ret = 0;
    trainingPair * p = pair;
    int m = n;
    while(p && m != 0)
        {
        p->setTentativeRes(undecided);
        p->unset(b_tentativelysolved);
        p = p->next();
        --m;
        }
    while(pair && n != 0)
        {
        ++ret;
        switch(lemmatise(pair))
            {
            case wrong:
#if AMBIGUOUS
                assert(pair->getTentativeRes() == no || pair->getTentativeRes() == notme);
                if(pair->getTentativeRes() != notme)
                    { // don't count homographs that have an ok sibling
                    pair->addRule(this,InputRight,false);
                    }
                else
                    {
                    }
#else
                pair->addRule(this,InputRight,false);
#endif
                break;
            case right:
                assert(pair->getTentativeRes() == yes);
                // do opportunistically count homographs that are ok
                pair->addRule(this,InputRight,true);
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
        pair = pair->next();
        --n;
        }
    assert(!pair);
    return ret;;
    }

long nextprime(unsigned long g)
    {
    int i;
    unsigned long smalldivisor;
    static int byte[12]=
        {1,  2,  2,  4,    2,    4,    2,    4,    6,    2,  6};
    /*2-3,3-5,5-7,7-11,11-13,13-17,17-19,19-23,23-29,29-1,1-7*/
    unsigned long bigdivisor;
    if(g & 1)
        ++g; /* even -> uneven */
    smalldivisor = 2;
    i = 0;
    while((bigdivisor = g / smalldivisor) >= smalldivisor)
        {
        if(bigdivisor * smalldivisor == g)
            {
            ++g;
            smalldivisor = 2;
            i = 0;
            }
        else
            {
            smalldivisor += byte[i];
            if(++i > 10)
                i = 3;
            }
        }
    return g;
    }

long casesensitivehash(const char * cp)
    {
    long hash_temp = 0;
    while (*cp != '\0')
        {
        if(hash_temp < 0)
            hash_temp = (hash_temp << 1) +1;
        else
            hash_temp = hash_temp << 1;
        hash_temp ^= *cp;
        ++cp;
        }
    return hash_temp;
    }


void hash::rehash(int loadFactor/*1-100*/)
    {
    long oldsize = hash_size;
    hash_size = nextprime((100 * record_count)/loadFactor);
    vertex ** new_hash_table = new vertex * [hash_size];
    long i;
    for(i = 0;i < hash_size;++i)
        new_hash_table[i] = 0;
    if(hash_table)
        {
        for(i = oldsize;i > 0;)
            {
            vertex * r = hash_table[--i];
            while(r)
                {
                long Key = key(r->itsTxt());
                vertex ** head = new_hash_table + Key;
                vertex * Next = r->getNext();
                r->setNext(*head);
                r->setHead(head);
                *head = r;
                r = Next;
                }
            }
        }
    delete [] hash_table;
    hash_table = new_hash_table;
    }

hash::hash(long size):record_count(0),next(0)
    {
    hash_size = nextprime(size);
    hash_table = new vertex * [hash_size];
    long i;
    for(i = 0;i < hash_size;++i)
        hash_table[i] = 0;
    ++HashCount;
    }

hash::~hash()
    {
    int i;
    int n = 0;
    for(i = 0;i < hash_size;++i)
        {
        vertex * ps = hash_table[i];
        while(ps)
            {
            ++n;
            ps->Hash = 0;
            ps = ps->getNext();
            }
        }
    delete [] hash_table;
    // do we want to delete the strngs as well?
    delete next;
    --HashCount;
    }

long hash::key(const char * ckey)
    {
    return ((unsigned long)casesensitivehash(ckey)) % hash_size;
    }

vertex * hash::find(const char * ckey,vertex **& head)
    {
    long Key = key(ckey);
    head = hash_table + Key;
    vertex * p;
    for(p = *head;p && strcmp(p->itsTxt(),ckey);p = p->getNext())
        ;
    return p;
    }

vertex ** hash::convertToList(int & N)
    {
    N = record_count;
    vertex ** ret = new vertex * [record_count];
    int i;
    int n = 0;
    for(i = 0;i < hash_size;++i)
        {
        vertex * ps = hash_table[i];
        while(ps)
            {
            ps->incRefCnt();
            ret[n++] = ps;
            ps = ps->getNext();
            }
        }
    return ret;
    }

int hash::forall(forallfunc fnc)
    {
    int i;
    int n = 0;
    for(i = 0;i < hash_size;++i)
        {
        vertex * ps = hash_table[i];
        while(ps)
            {
            ++n;
            (ps->*fnc)();
            ps = ps->getNext();
            }
        }
    return n;
    }

void strng::checkIntegrity()
    {
    assert(Txt);
    }

bool strng::hasWildCard()
    {
    return Txt[0] != START || Txt[strlen(Txt)-1] != END || strchr(Txt,ANY);
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

vertex * hash::getVertex(vertex * Rule,bool & New)
    {
    vertex * ret = 0;
    New = false;
    vertex ** head;
    long lf = loadfactor();
    if(lf > 100)
        {
        rehash(60);
        }
    vertex * p = find(Rule->Pattern.itsTxt(),head);
    if(p)
        p = p->findReplacement(Rule);
    if(p)
        {
        ret = p;
        }
    else
        {
        ret = new vertex(Rule,this);
        ret->setNext(*head);
        *head = ret;
        ret->setHead(head);
        New = true;
        inc();
        }
    return ret;
    }


int trainingPair::cmpWord(const trainingPair * B) const
    {
    int ret = (int)this->wordlength - (int)B->wordlength;
    if(!ret)
        ret = strncmp(this->Word,B->Word,this->wordlength);
    return ret;
    }

int trainingPair::cmpLemma(const trainingPair * B) const
    {
    int ret = (int)this->lemmalength - (int)B->lemmalength;
    if(!ret)
        ret = strncmp(this->LemmaHead,B->LemmaHead,this->lemmalength);
    return ret;
    }

#if WORDCLASS
int trainingPair::cmpWordClass(const trainingPair * B) const
    {
    int ret = (int)this->wordclasslength - (int)B->wordclasslength;
    if(!ret)
        ret = strncmp(this->WordClass,B->WordClass,this->wordclasslength);
    return ret;
    }
#endif
#if LEMMACLASS
int trainingPair::cmpLemmaClass(const trainingPair * B) const
    {
    int ret = (int)this->lemmaclasslength - (int)B->lemmaclasslength;
    if(!ret)
        ret = strncmp(this->LemmaClass,B->LemmaClass,this->lemmaclasslength);
    return ret;
    }
#endif

#if LEMMAINL
int trainingPair::cmpFreq(const trainingPair * B) const
    {
    int ret = (int)this->Inl - (int)B->Inl;
    if(!ret)
        ret = (int)this->Lemma_Inl - (int)B->Lemma_Inl;
    return ret > 0 ? 1 : ret < 0 ? -1 : 0;
    }
#endif

void trainingPair::fprint(FILE * famb)
    {
    fprintf(famb,"%.*s\t%.*s",(int)(this->wordlength),this->Word,(int)(this->lemmalength),this->LemmaHead);
    //if(this->Lemma && this->isset(b_wrong))
    if(this->Lemma && this->getRes() != yes)
        {
        fprintf(famb,"\t%s",this->Lemma);
        if(V)
            {
            fprintf(famb,"\t");
            V->print1(famb);
            }
        }
#if WORDCLASS
    fprintf(famb,"\t%.*s",this->wordclasslength,this->WordClass);
#endif
#if LEMMACLASS
    fprintf(famb,"\t%.*s",this->lemmaclasslength,this->LemmaClass);
#endif
    fprintf(famb,"\n");
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

void trainingPair::fprintTraining(FILE * fp
#if WRITEINVERTED
                                  ,FILE * fpinverted
#endif
                                  )
    {
    if(UTF8 
      ? isUpper(UTF8char(this->Word,UTF8)) && !isUpper(UTF8char(this->LemmaHead,UTF8)) 
      : isUpperISO(this->Word) && !isUpperISO(this->LemmaHead)
      )
        {
        char * lemma;
        if(isAllUpper(this->Word,this->wordlength))
            {
            size_t length = this->lemmalength;
            const char * upper = changeCase(this->LemmaHead,false,length);
            lemma = new char[length+1];
            strcpy(lemma,upper);
            }
        else
            {
            lemma = new char[this->lemmalength+1];
            strncpy(lemma,this->LemmaHead,this->lemmalength);
            lemma[this->lemmalength] = '\0';
            char * plemma = lemma;
            const char * pword = this->Word;
            size_t len = this->wordlength;
            if(UTF8)
                {
                while(len-- > 0 && *plemma)
                    {
                    char * lemma = plemma;
                    int LEMMA = upperEquivalent(getUTF8char((const char *&)plemma,UTF8));
                    int WORD =  getUTF8char(pword,UTF8);
                    if(LEMMA == WORD)
                        {
                        char NLEMMA[7];
                        int n = UnicodeToUtf8(LEMMA,NLEMMA,sizeof(NLEMMA)-1);
                        if(n == plemma - lemma)
                            {
                            strncpy(lemma,NLEMMA,n);
                            // Notice: no '\0'
                            }
                        else
                            break;
                        }
                    else
                        break;
                    }
                }
            else
                {
                while(len-- > 0 && *plemma && isUpperISO(pword))
                    {
                    toUpperISO(plemma);
                    ++plemma;
                    ++pword;
                    }
                }
            }
        fprintf(fp,"%.*s\t%s",(int)(this->wordlength),this->Word,lemma);
#if WRITEINVERTED
        fprintf(fpinverted,"%s\t%s",strnrev(this->Word,this->wordlength),Strrev(lemma));         
#endif
        delete [] lemma;
        }
    else
        {
        fprintf(fp,"%.*s\t%.*s",(int)(this->wordlength),this->Word,(int)(this->lemmalength),this->LemmaHead);
        char * lemmainverted = new char[this->lemmalength + 1];
#if WRITEINVERTED
        strcpy(lemmainverted,strnrev(this->LemmaHead,this->lemmalength));
        fprintf(fpinverted,"%s\t%s",strnrev(this->Word,this->wordlength),lemmainverted);
#endif
        delete [] lemmainverted;
        }
#if WORDCLASS
    fprintf(fp,"\t%.*s",this->wordclasslength,this->WordClass);
#endif
#if LEMMACLASS
    fprintf(fp,"\t%.*s",this->lemmaclasslength,this->LemmaClass);
#endif
#if LEMMAINL
    fprintf(fp,"\t%d\t%d",this->Inl,this->Lemma_Inl);
#endif
    fprintf(fp,"\n");
#if WRITEINVERTED
#if WORDCLASS
    fprintf(fpinverted,"\t%.*s",this->wordclasslength,this->WordClass);
#endif
#if LEMMACLASS
    fprintf(fpinverted,"\t%.*s",this->lemmaclasslength,this->LemmaClass);
#endif
#if LEMMAINL
    fprintf(fpinverted,"\t%d\t%d",this->Inl,this->Lemma_Inl);
#endif
    fprintf(fpinverted,"\n");
#endif
    }

void trainingPair::fprintAll(FILE * famb)
    {
    fprintf(famb,"%.*s\t%.*s",(int)(this->wordlength),this->Word,(int)(this->lemmalength),this->LemmaHead);
    if(this->Lemma && this->getRes() != yes)
        {
        fprintf(famb,"\t%s",this->Lemma);
        if(V)
            {
            fprintf(famb,"\t");
            V->print1(famb);
            }
        }
#if WORDCLASS
    fprintf(famb,"%.*s",this->wordclasslength,this->WordClass);
#endif
#if LEMMACLASS
    fprintf(famb,"%.*s",this->lemmaclasslength,this->LemmaClass);
#endif
#if LEMMAINL
    fprintf(famb,"\t%d\t%d",this->Inl,this->Lemma_Inl);
#endif
    if(this->isset(b_ambiguous))
        fprintf(famb," ambiguous");
    if(this->isset(b_skip))
        fprintf(famb," skip");
    fprintf(famb,"\n");
    }

void trainingPair::makeCandidateRules(hash * Hash,vertex * parent,bool alreadyRight /*20101206 used to calm makeRuleEx*/,optionStruct * options)
    {
    trainingPair * tp = this;
    while(tp)
        {
        tp->makeRuleEx(Hash,parent,alreadyRight,options);
        tp = tp->next();
        }
    }

int printRules(node * nd
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
             )
    {
    int n = 0;
    while(nd)
        {
        strng * nLL = 0;
        strng * nRR = 0;
#if BRACMATOUTPUT
        strng * snode;
#endif
        // pat is the complete pattern, not an optimized pattern that only 
        // specifies the differences with its parent pattern.
        // pat starts with '^' and ends with '$'.
        // Wildcards are the character in the constant ANY (currently ':')
        if(options->verbose())
            printf("pat [%s] rep [%s]\n",nd->V->Pattern.itsTxt(),nd->V->Replacement.itsTxt());
        strng * pat = nd->V->Pattern.substr(1,strlen(nd->V->Pattern.itsTxt()) - 2);
        strng * rep = nd->V->Replacement.substr(1,strlen(nd->V->Replacement.itsTxt()) - 2);
        strng * patreps[100];
        unsigned int i;
        for(i = 0;i < sizeof(patreps)/sizeof(patreps[0]);++i)
            patreps[i] = 0;
        if(nd->IfPatternSucceeds)
            {
#if BRACMATOUTPUT
            if(fobra)
                {
                snode = makeNode(patreps,nr,pat,rep,L,R,&nLL,&nRR);
                fprintf(fobra,"(");
                }
            else
                makeNode(patreps,nr,pat,rep,L,R,&nLL,&nRR);
#else
            makeNode(patreps,pat,rep,L,R,&nLL,&nRR);
#endif
            }
        else
            {
#if BRACMATOUTPUT
            snode = makeNode(patreps,nr,pat,rep,L,R,NULL,NULL);
#else
            makeNode(patreps,pat,rep,L,R,NULL,NULL);
#endif
            }
        for(i = 0;patreps[i];i+=2)
            ;
        //*
        if(i > 4)
            {
            i -= 2;
            strng * tmppat = patreps[i];
            strng * tmprep = patreps[i+1];
            while(i > 2)
                {
                patreps[i] = patreps[i - 2];
                patreps[i+1] = patreps[i - 1];
                i -= 2;
                }
            patreps[2] = tmppat;
            patreps[3] = tmprep;
            }
        //*/
        if(folem)
            {
            fprintf(folem,"%d\t",ind);
            if(patreps[0])
                {
                for(i = 0;patreps[i+2];i+=2)
                    {
                    fprintf(folem,"%s\t%s\t",patreps[i]->itsTxt(),patreps[i+1]->itsTxt());
                    }
                fprintf(folem,"%s\t%s",patreps[i]->itsTxt(),patreps[i+1]->itsTxt());
                }
            fprintf(folem,"\n");
            }
        for(i = 0;patreps[i];++i)
            delete patreps[i];

        delete pat;
        delete rep;
#if RULESASTEXTINDENTED
        nd->V->printRule(fo,ind,nr);
#endif
#if BRACMATOUTPUT
        if(fobra)
            {
            fprintf(fobra,"%s\n",snode->itsTxt());
            delete snode;
            }
#endif
        if(nd->IfPatternSucceeds)
            {
#if BRACMATOUTPUT
            if(fobra)
                fprintf(fobra,",");
#endif
            strng nL(L);
            strng nR(nRR);
            nL.cat(nLL,(const strng *)0);
            nR.cat(R,(const strng *)0);
            n += printRules
                ( nd->IfPatternSucceeds
#if RULESASTEXTINDENTED
                , fo
#endif
#if BRACMATOUTPUT
                , fobra
#endif
                , folem
                , ind+1
                , &nL
                , &nR
                , nr
                , options
                );
#if BRACMATOUTPUT
            if(fobra)
                fprintf(fobra,")\n");
#endif
            }
        delete nLL;
        delete nRR;
        nd = nd->IfPatternFails;
        }
    return n;
    }

void node::splitTrainingPairList(trainingPair * all,trainingPair **& pNotApplicable,trainingPair **& pWrong,trainingPair **& pRight,optionStruct * options)
    {
    while(all)
        {
        trainingPair * nxt = all->next();
        all->setNext(0);
        char * mask = 0;
        matchResult res = V->lemmatisem(all,&mask,0,options);
        if(mask)
            {
            all->takeMask(mask);
            assert(res != failure);
            }
        else
            {
            assert(res == failure);
            }
        switch(res)
            {
            case failure:
                *pNotApplicable = all;
                pNotApplicable = all->pNext();
                break;
            case wrong:
                *pWrong = all;
                pWrong = all->pNext();
                break;
            default:
                *pRight = all;
                pRight = all->pNext();
            }
        all = nxt;
        }
    }

bool node::compatibleSibling(node * sib)
    {
    if(this->V->dif(sib->V) != dif_incompatible)
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

#if AMBIGUOUS
trainingPair * trainingPair::makeWrongAmbiguousIfRightPresent(trainingPair *& Ambiguous)
    {
    trainingPair * tp;
    // Are there any unambiguously wrongly lemmatized training pairs?
    for(tp = this;tp;tp = tp->Next) 
        {
        if(tp->Alt == tp)
            {
            break;
            }
        else
            { // Check whether any of the alternatives was lemmatized correctly
            trainingPair * q = tp->Alt;
            while(q != tp && q->getRes() != yes)
                {
                q = q->Alt;
                }
            if(q == tp)
                { // None of these homographs are lemmatized correctly.
                break;
                }
            }
        }
    if(tp)
        { // There is at least one training pair that is wrongly lemmatized
          // and that has no homograph. 
        Ambiguous = 0;
        return this;
        }
    else
        {
      //  printf("All ambiguous\n");
        Ambiguous = this;
        return 0;
        }
    }
#endif

/*
static long checkPV(vertex ** pv,int N)
    {
    long R = pv[0]->R();
    long W = pv[0]->W();
    long res = R + W;
    for(int i = 1;i < N;++i)
        {
        assert(pv[i]->R() == R && pv[i]->W() == W);
        }
    return res;
    }

static void qqsort(vertex ** pv,size_t num,size_t width,int (__cdecl *compare )(const void *, const void *))
    {
    for(int i = 1;i < num;++i)
        {
        if(comp(*pv,pv[i]) > 0)
            {
            vertex * tmp = pv[0];
            pv[0] = pv[i];
            pv[i] = tmp;
            }
        }
    }
*/


#if AMBIGUOUS
void node::init(trainingPair ** allRight,trainingPair ** allWrong/*,trainingPair ** allAmbiguous*/,int level,optionStruct * options)
#else
void node::init(trainingPair ** allRight,trainingPair ** allWrong,int level/*,vertex ** pvp,int Np*/,optionStruct * options)
#endif
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

    trainingPair * NotApplicableRight = NULL;
    trainingPair * NotApplicableWrong = NULL;
    trainingPair * Wrong = NULL;

    trainingPair ** pNotApplicableRight = &NotApplicableRight; // Pattern fails.
    trainingPair ** pNotApplicableWrong = &NotApplicableWrong; // Pattern fails.
    trainingPair ** pWrong = &Wrong; // Pattern succeeds and replacement is wrong.
    trainingPair ** pRight = &this->Right; // Pattern succeeds and replacement is right.
    assert(this->Right == NULL);

    for(trainingPair * tp = *allRight;tp;tp = tp->next())
        {
        tp->unset(b_solved);
        }
    for(trainingPair * tp = *allWrong;tp;tp = tp->next())
        {
        tp->unset(b_solved);
        }

    this->splitTrainingPairList(*allRight,pNotApplicableRight,pWrong,pRight,options);
    this->splitTrainingPairList(*allWrong,pNotApplicableWrong,pWrong,pRight,options);

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
        int N;
        hash Hash(1000);
        Wrong->makeCandidateRules(&Hash,this->V,false,options);
        if(this->Right)
            this->Right->makeCandidateRules(&Hash,this->V,true,options);

        vertex ** pv = Hash.convertToList(N);
        if(N == 0)
            {
            delete [] pv;
            return;
            }
        assert(N > 0);
        int maxN = N;
        if(maxN < 16)
            maxN = 16;
        int wpart = -1; // -1 indicating: no upper bound
        int rpart = -1;
        /* If wpart and rpart are positive, the rules are applied to just a 
           subset of training pairs. This saves time and memory.
        */

        // hack if number of pairs is very big
        int ntestpairWrong = Wrong->count();
        int ntestpairRight = this->Right ? this->Right->count() : 0;
        double fraction = (double)MAXPAIRS / (double)(ntestpairWrong + ntestpairRight);
        if(fraction < 1.0)
            { // Define upper bound to make it managable memory-wise
            if(options->verbose())
                {
                printf("%d < %d N=   %d\n",MAXPAIRS,(ntestpairWrong + ntestpairRight),N);
                }
            fprintf(stderr,"***** FAIRLY BIG THIS TRAINING SET IS! %d < %d N=   %d\n",MAXPAIRS,(ntestpairWrong + ntestpairRight),N);
            wpart = (int)(fraction * (double)wpart);
            rpart = MAXPAIRS - wpart;
            }
        //:hack

        // Test all pairs (up to upper bound) on all candidate rules.
        //printf("Test all pairs (up to upper bound) on all candidate rules. %d\n",N);
        if(options->verbose())
            {
            printf("             %d  \r",N);
            }
        for(int i = 0;i < N;++i)
            {
            if(options->verbose())
                {
                printf("%d \r",i);
                }
            // Reset all counters in all candidate rules.
            pv[i]->nlemmatiseStart();
            pv[i]->nlemmatise(Wrong,wpart,false);
            pv[i]->nlemmatise(this->Right,rpart,true);
            }
        node ** pnode = &this->IfPatternSucceeds;
        ptrdiff_t first = 0;
        ptrdiff_t lastN = N; // Rules having position on list beyond lastN do not
                       // apply to any of the remaining pairs.
        do
            {
            if(first >= lastN)
                {
                fprintf(stderr,"***** first:%ld >= lastN:%ld\n",(long)first,(long)lastN);
                fprintf(stderr,"***** (This happens if one or more wrongly lemmatised words aren't handled by any candidate rules.)\n");
                fprintf(stderr,"***** (List of candidate rules:)\n");
                fprintf(stderr,"***** pv:");
                for(int i = 0;i < lastN;++i)
                    {
                    fprintf(stderr,"i=%d:",i);
                    pv[i]->print1(stderr);
                    fprintf(stderr,"\n");
                    }
                fprintf(stderr,"***** (List of unmatched words that need better lemmatisation rule(s):)\n");
                fprintf(stderr,"\n***** Wrong:\n");
                Wrong->printAll(stderr,"unmatched words that need better lemmatisation rule(s)",'\n');
                if(options->verbose())
                    {
                    printf("\nAFFIXTRAIN failed\n");
                    }
                fprintf(stderr,"\n***** AFFIXTRAIN failed\n");
                exit(-1);
                }
            assert(first < lastN);
            vertex ** pvf = pv+first;
            vertex ** pvN = pv+lastN;
            if (comp == comp_parms)
                {
                for (vertex ** pvi = pvf; pvi < pvN; ++pvi)
                    {
                    computeWeight(*pvi);
                    }
                if (options->expensiveInfix())
                    {
                    for (vertex ** pvi = pvf; pvi < pvN; ++pvi)
                        {
                        (*pvi)->adjustWeight();
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
            *pnode = new node(*pvf++);
#if AMBIGUOUS
            (*pnode)->init(&this->Right,&Wrong,level+1,options);
#else
            (*pnode)->init(&this->Right,&Wrong,level+1/*,pv+first,lastN-first*/,options);
#endif
#if _NA
            int outputR = (this->Right ? this->Right->count() : 0); 
            int outputW = (Wrong ? Wrong->count() : 0);
#endif
            if(wpart < 0)
                {
                for(vertex ** pvi = pvf;pvi < pvN;++pvi)
                    {
                    if(  (*pvi)->R__R == 0 
                      && (*pvi)->W__R == 0 
                      && (*pvi)->R__W == 0
                      && (*pvi)->W__W == 0
                      )
                        { // rule has become irrelevant
                        (*pvi)->destroy();
                        --pvN;
                        *pvi = *pvN;
                        }
#if _NA
                    else
                        {
                        if(outputR < outputW)
                            {
                            (*pvi)->adjustNotApplicableCountsByRecalculatingR_NA(this->Right,outputR+outputW);
                            }
                        else
                            {
                            (*pvi)->adjustNotApplicableCountsByRecalculatingW_NA(Wrong,outputR+outputW);
                            }
                        }
#endif
                    }
                if(pvf == pvN)
                    {
                    if(options->verbose()) /* This is not really an error. */
                        fprintf(stderr,"\n***** destroyed all remaining untested rule candidates because they have become inapplicable\n");
                    }
                }
            if(options->verbose() && wpart >= 0)
                {
                if(lastN > pvN - pv)
                    {
                    printf("Deleted %ld of original %d candidates (%f%%)\n",(long int)(lastN - (pvN - pv)),N,(100.0*(lastN - (pvN - pv))/N));
                    lastN = pvN - pv;
                    }
                }
            else
                lastN = pvN - pv;

            ++first;
            if(first >= lastN)
                {
                if(options->verbose()) /* This is not really an error. */
                    fprintf(stderr,"\n***** DANGER first %ld >= lastN %ld\n",(long)first,(long)lastN);
                }
            // hack:

            int CnT;
#if AMBIGUOUS
            // Take apart the wrongly lemmatised homographs that have a correctly lemmatised sibling. 
            Wrong = Wrong->makeWrongAmbiguousIfRightPresent(Ambiguous);
#endif
            if(  wpart >= 0 
              && Wrong 
              &&    ( CnT = 
                      Wrong->count() 
                    + (this->Right ? this->Right->count() : 0) 
                    )
                 <= MAXPAIRS
              )
                {
                if(options->verbose())
                    {
                    printf("%d > %d lastN=%ld\n",MAXPAIRS,CnT,(long)lastN);
                    }
                fprintf(stderr,"***** FAIRLY BIG THIS TRAINING SET IS! BUT NOW %d > %d lastN=%ld\n",MAXPAIRS,CnT,(long)lastN);
                // Test all remaining pairs on all remaining candidate rules.
                for(ptrdiff_t i = first;i < lastN;++i)
                    {
                    // Reset all counters in remaining candidate rules.
                    pv[i]->nlemmatiseStart();
                    pv[i]->nlemmatise(Wrong,-1,false);
                    pv[i]->nlemmatise(this->Right,-1,true);
                    }
                wpart = -1;
                }
            // :hack
            pnode = &(*pnode)->IfPatternFails;
            }
        while(Wrong);

        *pnode = 0;
        for(int i = 0;i < lastN;++i)
            {
            pv[i]->destroy();
            }
        delete [] pv;
        }
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
        switch(lemmatise(pair))
            {
            case wrong:
            case right:
                break;
            default:
                ++R__NA;
            }
        pair = pair->next();
        }
    W__NA = total - (R__NA + R__R + R__W + W__R + W__W);
    }

void vertex::adjustNotApplicableCountsByRecalculatingW_NA(trainingPair * pair,int total)
    {
    W__NA = 0;
    while(pair)
        {
        switch(lemmatise(pair))
            {
            case wrong:
            case right:
                break;
            default:
                ++W__NA;
            }
        pair = pair->next();
        }
    R__NA = total - (R__R + R__W + W__NA + W__R + W__W);
    }

#endif
