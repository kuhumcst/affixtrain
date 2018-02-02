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
#include "trainingpair.h"
#include "ruletemplate.h"
#include "strng.h"
#include "node.h"
#include "utf8func.h"
#include "isofunc.h"
#include "optionaff.h"
#include "shortrulepair.h"
#include "vertex.h"
#include "hashtable.h"


#if PRUNETRAININGPAIRS
FILE * fpmourn = 0;
#endif

int TrainingPairCount = 0;
trainingPair::trainingPair():Next(0),Mask(0),Lemma(0),V(0),applicableRules(0),ambs(undecided),tentativeAmbs(undecided)
#if PRUNETRAININGPAIRS
    ,likes(1)
#endif
    {
    ++TrainingPairCount;
    }

void trainingPair::deepDelete()
    {
    deleteRules();
    delete [] Lemma;
    Lemma = NULL;
    delete [] Mask;
    Mask = NULL;
    Next = 0;
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

void trainingPair::setLemma(char * lemma)
    {
    if(Lemma)
        delete [] Lemma;
    Lemma = dup(lemma);
    }

void trainingPair::takeMask(char * mask)
    {
    if(Mask)
        delete [] Mask;
    Mask = dup(mask);
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

void trainingPair::printMore(FILE * f)const
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

void trainingPair::fprintTraining(FILE * fp
#if WRITEINVERTED
                                  ,FILE * fpinverted
#endif
                                  )
    {
    if( globUTF8 
      ? isUpper(UTF8char(this->Word,globUTF8)) && !isUpper(UTF8char(this->LemmaHead,globUTF8)) 
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
            if(globUTF8)
                {
                while(len-- > 0 && *plemma)
                    {
                    char * lemma = plemma;
                    int LEMMA = (int)upperEquivalent(getUTF8char((const char *&)plemma,globUTF8));
                    int WORD =  getUTF8char(pword,globUTF8);
                    if(LEMMA == WORD)
                        {
                        char NLEMMA[7];
                        int n = UnicodeToUtf8(LEMMA,NLEMMA,sizeof(NLEMMA)-1);
                        if(n == plemma - lemma)
                            {
                            strncpy(lemma,NLEMMA,(size_t)n);
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

static int storeRule(hashTable * Hash, shortRulePair * Rule/*, vertex *& V*/)
    {
    vertex FullRule(Rule);
    //vertex * V;
    bool New;
    /*V =*/ Hash->getVertex(&FullRule, New);
    if (New)
        {
        return 1;
        }
    else
        {
        return 0;
        }
    }

void trainingPair::makeCandidateRules(hashTable * Hash,vertex * parent,bool alreadyRight,optionStruct * options)
    {
    trainingPair * tp = this;
    while(tp)
        {
        char similarArray[1000];
        similData SimilData(similarArray);
        SimilData.match(tp);
        const char * predef = tp->getMask();
        SimilData.mergeTemplates(predef, options);
        shortRulePair Rule(tp, &SimilData);
        if (Rule.checkRule(tp, parent))
            {
            //vertex * e;
            storeRule(Hash, &Rule/*, e*/);
            }
        else
            {
            int nr = tp->makeCorrectRules(Hash, &SimilData, similarArray, parent, 1, options->maxRecursionDepthRuleCreation(), options);
            if (nr == 0 && !alreadyRight)
                {
                if (options->verbose() > 4)
                    { // This is by design. Increasing RECURSE will eventually help.
                    fprintf(stderr, "Error (makeCandidateRules): Cannot construct rule for training pair ");
                    tp->print(stderr);
                    fprintf(stderr, " Based on parent ");
                    parent->print1(stderr);
                    fprintf(stderr, "\n");
                    }
                }
            }
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
    if (options->verbose() > 6)
        printf("ind==(%d)\n",ind);
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
        if (options->verbose() > 6)
            printf("pat [%s] rep [%s]\n",nd->V->Pattern.itsTxt(),nd->V->Replacement.itsTxt());
        strng * pat = nd->V->Pattern.substr(1,(ptrdiff_t)strlen(nd->V->Pattern.itsTxt()) - 2);
        strng * rep = nd->V->Replacement.substr(1, (ptrdiff_t)strlen(nd->V->Replacement.itsTxt()) - 2);
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
    if (options->verbose() > 6)
        printf("ind==(%d)DONE n==(%d)\n",ind,n);
    return n;
    }

int trainingPair::makeCorrectRules(hashTable * Hash, ruleTemplate * Template, const char * similar, vertex * parent, int mlow, int recurse, optionStruct * options)
    {
    /*In the template, replace one = with a # and construct a rule based on
    this new template. If the rule succeeds, store it. Repeat this for all
    places were replacement is allowed. (What should be allowed is a matter
    of heuristics.)*/
    int ret = 0;
    bool different = true;
    int anihilatedGuards;
    ruleTemplate locTemplate;
    locTemplate.copy(Template);

    mlow = 1;

    for (int m = mlow
         ; locTemplate.makebigger(m, anihilatedGuards, options)
         ;   ++m
         , locTemplate.copy(Template)
         )
        {
        shortRulePair Rule(this, &locTemplate);
        if (Rule.checkRule(this, parent))
            {
            different = false;
            //vertex * e;
            storeRule(Hash, &Rule/*, e*/);
            ++ret; // Increment for each new or already existing rule that fits
            }
        }
    /*If no rule with one more character in the search pattern succeeds,
    repeat the procedure with two instead of one additional character,
    and so on.*/
    if (different || --recurse > 0)
        {
        locTemplate.copy(Template);
        for (int m = mlow
             ; locTemplate.makebigger(m, anihilatedGuards, options)
             ;   ++m
             , locTemplate.copy(Template)
             )
            {
            /*Recurse with template that has one more #*/
            ret += makeCorrectRules(Hash, &locTemplate, similar, parent, m + 1 - anihilatedGuards, recurse, options);
            }
        }
    return ret;
    }

#if PRUNETRAININGPAIRS
bool trainingPair::fewerLikesThan(int thresh) const
    {
    bool ret = !applicableRules || applicableRules->fewerLikesThan(thresh);
    return ret;
    }
#endif
#if _NA
int trainingPair::notLemmatizedBy(vertex * V)
    {
    int ret = 0;
    for(trainingPair * p = this;p;p = p->Next)
        {
        if(!p->applicableRules->has(V))
            ++ret;
        }
    return ret;
    }
#endif
#if PRUNETRAININGPAIRS
void trainingPair::mourn()
    {
    for(trainingPair * p = this;p;p = p->Next)
        {
        p->printMore(fpmourn);
        p->deepDelete();
        }
    }

#endif
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
