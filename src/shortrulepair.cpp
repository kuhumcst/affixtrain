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
#include "shortrulepair.h"
#include "ruletemplate.h"
#include "utf8func.h"
#include "isofunc.h"
#include "trainingpair.h"
#include "vertex.h"


int ShortRulePairCount = 0;

shortRulePair::shortRulePair(trainingPair * trainingpair, ruleTemplate * Template)
    {
    bool star = false;
    int i = 0;
    int i2 = 0;
    int j = 0;
    int j2 = 0;
    size_t k = 0;
    size_t k2 = 0;
    int inc;
    while (Template->itsNPatternArray()[i] && Template->itsNReplacementArray()[i2])
        {
        while (UTF8 && Template->itsNPatternArray()[i] && Template->itsNPatternArray()[i] != equal)
            {
            star = false;
            inc = copyUTF8char(trainingpair->itsWord() + k, patternArray + j);
            j += inc;
            k += inc;
            ++i;
            }
        while (Template->itsNPatternArray()[i] && Template->itsNPatternArray()[i] != equal)
            {
            star = false;
            patternArray[j] = trainingpair->itsWord()[k];
            ++j;
            ++k;
            ++i;
            }
        while (UTF8 && Template->itsNReplacementArray()[i2] && Template->itsNReplacementArray()[i2] != equal)
            {
            star = false;
            inc = copyUTF8char(trainingpair->itsLemmaHead() + k2, replacementArray + j2);
            j2 += inc;
            k2 += inc;
            ++i2;
            }
        while (Template->itsNReplacementArray()[i2] && Template->itsNReplacementArray()[i2] != equal)
            {
            star = false;
            replacementArray[j2] = trainingpair->itsLemmaHead()[k2];
            ++j2;
            ++k2;
            ++i2;
            }
        while (Template->itsNPatternArray()[i] == equal && Template->itsNReplacementArray()[i2] == equal)
            {
            if (!star)
                {
                patternArray[j++] = ANY;
                replacementArray[j2++] = ANY;
                patternArray[j] = '\0';
                replacementArray[j2] = '\0';
                star = true;
                }
            k += skipUTF8char(trainingpair->itsWord() + k);
            k2 += skipUTF8char(trainingpair->itsLemmaHead() + k2);
            i++;
            i2++;
            }
        }
    while (UTF8 && Template->itsNPatternArray()[i])
        {
        inc = copyUTF8char(trainingpair->itsWord() + k, patternArray + j);
        j += inc;
        k += inc;
        ++i;
        }
    while (Template->itsNPatternArray()[i])
        {
        patternArray[j] = trainingpair->itsWord()[k];
        ++j;
        ++k;
        ++i;
        }
    while (UTF8 && Template->itsNReplacementArray()[i2])
        {
        inc = copyUTF8char(trainingpair->itsLemmaHead() + k2, replacementArray + j2);
        j2 += inc;
        k2 += inc;
        ++i2;
        }
    while (Template->itsNReplacementArray()[i2])
        {
        replacementArray[j2] = trainingpair->itsLemmaHead()[k2];
        ++j2;
        ++k2;
        ++i2;
        }
    patternArray[j] = '\0';
    replacementArray[j2] = '\0';
    trim();
    ++ShortRulePairCount;
    }

static int cnt(char * s, int k)
    {
    int ret = 0;
    while (*s)
        {
        if (*s == k)
            ++ret;
        ++s;
        }
    return ret;
    }

static void shiftleft(char * a)
    {
    while (*a)
        {
        a[0] = a[1];
        ++a;
        }
    }

static void shiftright(char * A)
    {
    char * a = A + strlen(A);
    a[1] = '\0';
    while (a > A)
        {
        a[0] = a[-1];
        --a;
        }
    }



void shortRulePair::trim()
    {
    char * A = patternArray;
    char * B = replacementArray;
    if (A[0] == ANY && B[0] == ANY)
        {
        int ca = cnt(A, ANY);
        int cb = cnt(B, ANY);
        if (ca == cb)
            {
            shiftleft(A);
            shiftleft(B);
            }
        else if (ca < cb)
            shiftleft(B);
        else
            shiftleft(A);

        }
    else
        {
        if (A[0] != ANY && A[0] != START)
            {
            shiftright(A);
            A[0] = START;
            }
        if (B[0] != ANY && B[0] != START)
            {
            shiftright(B);
            B[0] = START;
            }
        }
    size_t la = strlen(A);
    size_t lb = strlen(B);
    if (la && A[la - 1] == ANY && lb && B[lb - 1] == ANY)
        {
        A[la - 1] = '\0';
        B[lb - 1] = '\0';
        }
    else
        {
        if (la)
            {
            if (A[la - 1] == END)
                {
                if (!lb || (lb && B[lb - 1] == ANY))
                    {
                    B[lb] = END;
                    B[lb + 1] = '\0';
                    }
                }
            else if (A[la - 1] != ANY)
                {
                A[la] = END;
                A[la + 1] = '\0';
                if (!lb || (lb && B[lb - 1] == ANY))
                    {
                    B[lb] = END;
                    B[lb + 1] = '\0';
                    }
                }
            }
        if (lb)
            {
            if (B[lb - 1] == END)
                {
                if (!la || (la && A[la - 1] == ANY))
                    {
                    A[la] = END;
                    A[la + 1] = '\0';
                    }
                }
            else if (B[lb - 1] != ANY)
                {
                B[lb] = END;
                B[lb + 1] = '\0';
                if (!la || (la && A[la - 1] == ANY))
                    {
                    A[la] = END;
                    A[la + 1] = '\0';
                    }
                }
            }
        }
    }

bool shortRulePair::checkRule(trainingPair * trainingpair, vertex * parentPat)
    {
    vertex FullRule(this);
    if (  (  !parentPat
          || FullRule.Pattern.dif(&parentPat->Pattern) == dif_bigger
          )
       && FullRule.apply(trainingpair) == right
       )
        {
        return true;
        }
    else
        {
        return false;
        }
    }



shortRulePair::~shortRulePair(){ --ShortRulePairCount; }
