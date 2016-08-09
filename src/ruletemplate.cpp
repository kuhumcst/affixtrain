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
#include "ruletemplate.h"
#include "trainingpair.h"
#include "utf8func.h"
#include "isofunc.h"
#include "optionaff.h"

#include <string.h>
#include <assert.h>

int RuleTemplateCount;
const int unequal = '#';
const int equal = '=';
static const char * nil = "";
const char Start[2] = { START, 0 };
const char End[2] = { END, 0 };
static const int inil[1] = { 0 };
static const int iStart[2] = { START, 0 };
static const int iEnd[2] = { END, 0 };
const char StartAny[3] = { START, ANY, 0 };
const char AnyEnd[3] = { ANY, END, 0 };
const char StartAnyEnd[4] = { START, ANY, END, 0 };
static bool glob_wildcards = false;

void ruleTemplate::copy(ruleTemplate * Template)
    {
    strcpy(npatternArray, Template->npatternArray);
    strcpy(nreplacementArray, Template->nreplacementArray);
    }

ruleTemplate::ruleTemplate(){ ++RuleTemplateCount; }
ruleTemplate::~ruleTemplate(){ --RuleTemplateCount; }

void similData::match(trainingPair * trainingpair)
    {
    if (UTF8)
        {
        size_t wordLen = trainingpair->itsWordlength();
        size_t lemmaLen = trainingpair->itsLemmalength();
        int * iWord = new int[wordLen + 1];
        int * iLemmaHead = new int[lemmaLen + 1];
        wordLen = Utf8ToUnicode(iWord, trainingpair->itsWord(), wordLen);
        assert(lemmaLen > 0);
        lemmaLen = Utf8ToUnicode(iLemmaHead, trainingpair->itsLemmaHead(), lemmaLen);
        isimil(iWord, iWord + wordLen, iLemmaHead, iLemmaHead + lemmaLen, iStart, iEnd);
        delete[] iWord;
        delete[] iLemmaHead;
        }
    else
        {
        simil(trainingpair->itsWord(), trainingpair->itsWord() + trainingpair->itsWordlength(), trainingpair->itsLemmaHead(), trainingpair->itsLemmaHead() + trainingpair->itsLemmalength(), Start, End);
        }
    *psimilar = *ppattern = *preplacement = '\0';
    *nppattern = *npreplacement = '\0';
    }


bool similData::mergeTemplates(const char * predef, optionStruct * options)
    {
    char * n = npatternArray;
    char * r = nreplacementArray - 1;
    const char * p = predef;
    int eqs = 0;
    int reqs = 0;
    bool changed = false;
    while (*n)
        {
        if (*n == equal)
            {
            ++eqs;
            do
                {
                if (*++r == equal)
                    ++reqs;
                } while (reqs < eqs);
            if (*p == unequal)
                {
                *n = unequal;
                *r = unequal;
                changed = true;
                }
            }
        ++n;
        ++p;
        }
    if (options->suffixOnly())
        {
        suffix(npatternArray);
        suffix(nreplacementArray);
        }
    return changed;
    }


ptrdiff_t similData::simil(
    const char * const s1,
    const char * const n1,
    const char * const s2,
    const char * const n2,
    const char * const start,
    const char * const end)
    {
    const char * ls1;
    const char * s1l = NULL;
    const char * s1r = NULL;
    const char * s2l = NULL;
    const char * s2r = NULL;
    ptrdiff_t max;
    /*

    d o o d l a c h t e n            d o o d l a c h e n
    ^                   ^            ^                 ^
    s1                  n1           s2                n2



    d o o d l a c h t e n            d o o d l a c h e n
    ^                   ^            ^                 ^
    s1                  n1           s2                n2
    ^               ^                ^               ^
    ls1             lls1             ls2             lls2
    ^---------------^
    dif
    max



    d o o d l a c h t e n            d o o d l a c h e n
    ^                   ^            ^                 ^
    s1                  n1           s2                n2
    ^               ^                ^               ^
    ls1             lls1             ls2             lls2
    ^---------------^
    dif
    max
    ^               ^                ^               ^
    s1l             s1r              s2l             s2r
    */



    for (max = 0, ls1 = s1; ls1 + max < n1; ls1++)
        /* consider each character of s1 as candidate start character for a match */
        {
        const char * ls2;
        /* vergelijk met s2 */
        for (ls2 = s2; ls2 + max < n2; ls2++)
            {
            const char * lls1;
            const char * lls2;
            /* determine lenght of equal parts */
            for (lls1 = ls1, lls2 = ls2
                 ;     lls1 < n1
                 &&  lls2 < n2
                 && *lls1 == *lls2
                 ; lls1++, lls2++
                 )
                 ;
            /* adjust score, if needed */
            ptrdiff_t dif = lls1 - ls1;
            if (dif > max)
                {
                max = dif;
                /* remember end positions of left strings
                   and start positions of right strings */
                s1l = ls1;    /* start of longest common substring in s1 */
                s1r = lls1;   /*   end of longest common substring in s1 */
                s2l = ls2;    /* start of longest common substring in s2 */
                s2r = lls2;   /*   end of longest common substring in s2 */
                }
            }
        }
    if (max)
        {
        /*
        The longest sequence of characters occurring in s1 and s2 is found.
        The number of characters in the common substring is <max>.
        */
        if (s1 == s1l)
            {
            /*Potential problematic pattern: a replacement for no pattern.
            (Insertion).
            That only works if the replacement is at the start or the end
            of the word. (prefix or suffix)*/
            if (!*start)
                {
                /*Problem: empty pattern somewhere in the middle.
                Solution: borrow to the left or to the right,
                'reducing' the common substring.*/
                s1l++;
                s2l++;
                /*Caveat: By borrowing the left and right may touch.
                (max == 1 or 2)*/
                --max;
                }
            }
        }
    if (max)
        {
        if (s1r == n1)
            {
            if (!*end)
                {
                s1r--;
                s2r--;
                --max;
                }
            }
        }

    if (max)
        {
        if (s1 != s1l || s2 != s2l)
            {
            /* Recurse if the common substring does not start at the beginning of
            at least one of s1 and s2. */
            if (max)
                max += simil(s1, s1l, s2, s2l, start, nil);

            }
        else
            {
            if (*start && glob_wildcards)
                {
                *psimilar++ = *start;
                }
            }


        for (const char * s = s1l; s < s1r; ++s)
            /* This is the longest common substring. */
            {
            *nppattern++ = equal;
            *npreplacement++ = equal;
            }


        if (s1r != n1 || s2r != n2)
            {
            /* Recurse if the common substring does not end at the end of
            at least one of s1 and s2. */
            max += simil(s1r, n1, s2r, n2, nil, end);
            }
        else
            {
            if (*end && glob_wildcards)
                {
                *psimilar++ = *end;
                }
            }
        }
    else
        {
        /* The strings s1 and s2 are completely different */
        if (*start)
            {
            if (s1 < n1)
                *ppattern++ = *start; // ^ge
            if (s2 < n2)
                *preplacement++ = *start; // ^over
            }
        else if (!(ppattern == pattern && preplacement == replacement))
            {
            *ppattern++ = ANY;
            *preplacement++ = ANY;
            }

        if (s1 < n1) // pattern not nothing
            {
            for (const char * s = s1; s < n1; ++s)
                {
                *ppattern++ = *s;
                *nppattern++ = unequal;
                }
            }

        if (s2 < n2)
            {
            for (const char * s = s2; s < n2; ++s)
                {
                *preplacement++ = *s;
                *npreplacement++ = unequal;
                }
            }
        if (!*start && !*end && glob_wildcards)
            {
            *psimilar++ = '?';
            }

        if (*end)
            {
            *ppattern++ = *end;
            *preplacement++ = *end;
            }
        }
    return max;
    }

ptrdiff_t similData::isimil(
    const int * const s1,
    const int * const n1,
    const int * const s2,
    const int * const n2,
    const int * const start,
    const int * const end)
    {
    const int * ls1;
    const int * s1l = NULL;
    const int * s1r = NULL;
    const int * s2l = NULL;
    const int * s2r = NULL;
    ptrdiff_t max;
    /*

    d o o d l a c h t e n            d o o d l a c h e n
    ^                   ^            ^                 ^
    s1                  n1           s2                n2



    d o o d l a c h t e n            d o o d l a c h e n
    ^                   ^            ^                 ^
    s1                  n1           s2                n2
    ^               ^                ^               ^
    ls1             lls1             ls2             lls2
    ^---------------^
    dif
    max



    d o o d l a c h t e n            d o o d l a c h e n
    ^                   ^            ^                 ^
    s1                  n1           s2                n2
    ^               ^                ^               ^
    ls1             lls1             ls2             lls2
    ^---------------^
    dif
    max
    ^               ^                ^               ^
    s1l             s1r              s2l             s2r
    */



    for (max = 0, ls1 = s1; ls1 + max < n1; ls1++)
        /* consider each character of s1 as candidate start character for a match */
        {
        const int * ls2;
        /* compare with s2 */
        for (ls2 = s2; ls2 + max < n2; ls2++)
            {
            const int * lls1;
            const int * lls2;
            /* determine lenght of equal parts */
            for (lls1 = ls1, lls2 = ls2
                 ;    lls1 < n1
                 && lls2 < n2
                 && *lls1 == *lls2
                 ; lls1++, lls2++
                 )
                 ;
            /* adjust score, if needed */
            ptrdiff_t dif = lls1 - ls1;
            if (dif > max)
                {
                max = dif;
                /* remember end positions of left strings
                   and start positions of right strings */
                s1l = ls1;    /* start of longest common substring in s1 */
                s1r = lls1;   /*   end of longest common substring in s1 */
                s2l = ls2;    /* start of longest common substring in s2 */
                s2r = lls2;   /*   end of longest common substring in s2 */
                }
            }
        }
    if (max)
        {
        /*
        The longest sequence of characters occurring in s1 and s2 is found.
        The number of characters in the common substring is <max>.
        */
        if (s1 == s1l)
            {
            /*Potential problematic pattern: a replacement for no pattern.
            (Insertion).
            That only works if the replacement is at the start or the end
            of the word. (prefix or suffix)*/
            if (!*start)
                {
                /*Problem: empty pattern somewhere in the middle.
                Solution: borrow to the left or to the right,
                'reducing' the common substring.*/
                s1l++;
                s2l++;
                /*Caveat: By borrowing the left and right may touch.
                (max == 1 or 2)*/
                --max;
                }
            }
        }
    if (max)
        {
        if (s1r == n1)
            {
            if (!*end)
                {
                s1r--;
                s2r--;
                --max;
                }
            }
        }

    if (max)
        {
        if (s1 != s1l || s2 != s2l)
            {
            /* Recurse if the common substring does not start at the beginning of
            at least one of s1 and s2. */
            if (max)
                max += isimil(s1, s1l, s2, s2l, start, inil);

            }
        else
            {
            if (*start && glob_wildcards)
                {
                psimilar += UnicodeToUtf8(*start, psimilar, 1000);
                }
            }


        for (const int * s = s1l; s < s1r; ++s)
            /* This is the longest common substring. */
            {
            nppattern += UnicodeToUtf8(equal, nppattern, 1000);
            npreplacement += UnicodeToUtf8(equal, npreplacement, 1000);
            }


        if (s1r != n1 || s2r != n2)
            {
            /* Recurse if the common substring does not end at the end of
            at least one of s1 and s2. */
            max += isimil(s1r, n1, s2r, n2, inil, end);
            }
        else
            {
            if (*end && glob_wildcards)
                {
                psimilar += UnicodeToUtf8(*end, psimilar, 1000);
                }
            }
        }
    else
        {
        /*The strings s1 and s2 are completely different*/
        if (*start)
            {
            if (s1 < n1)
                {
                ppattern += UnicodeToUtf8(*start, ppattern, 1000);// ^ge
                }
            if (s2 < n2)
                {
                preplacement += UnicodeToUtf8(*start, preplacement, 1000);// ^over
                }
            }
        else if (!(ppattern == pattern && preplacement == replacement))
            {
            ppattern += UnicodeToUtf8(ANY, ppattern, 1000);
            preplacement += UnicodeToUtf8(ANY, preplacement, 1000);// ^over
            }

        if (s1 < n1) // pattern not nothing
            {
            for (const int * s = s1; s < n1; ++s)
                {
                ppattern += UnicodeToUtf8(*s, ppattern, 1000);
                nppattern += UnicodeToUtf8(unequal, nppattern, 1000);
                }
            }

        if (s2 < n2)
            {
            for (const int * s = s2; s < n2; ++s)
                {
                preplacement += UnicodeToUtf8(*s, preplacement, 1000);
                npreplacement += UnicodeToUtf8(unequal, npreplacement, 1000);
                }
            }
        if (!*start && !*end && glob_wildcards)
            {
            psimilar += UnicodeToUtf8('?', psimilar, 1000);
            }

        if (*end)
            {
            ppattern += UnicodeToUtf8(*end, ppattern, 1000);
            preplacement += UnicodeToUtf8(*end, preplacement, 1000);
            }
        }
    return max;
    }



void suffix(char * mask)
    {
    char * m;
    for (m = mask; *m && *m == equal; ++m)
        ;
    for (; *m; ++m)
        *m = unequal;
    }

/*
Reykjarnesbjarga
============#==#
============#==
Reykjarnesbjörg
Notice the substring "jar", which occurs twice in the inflected form.
The second "jar" changes to "jör" in the lemma. So these rules would all fail:
============#==#
*     a *a
============#==
*     ö *


===========##==#
*    ja *a
===========##==
*    jö *


============##=#
*     ar*a
============##=
*     ör*

===========###=#
*    jar*a
===========###=
*    jör*

In this case, the best solution is to sacrifice the last two =
============####
*     arga
============###
*     örg

*/

bool ruleTemplate::makebigger(int countdown, int & anihilatedGuards, optionStruct * options)
    {
    /*
    Change a = to a # if
    1) The next byte in replacement is # or
    2) The previous byte in pattern is = and the current byte is #
    */
    if (countdown == 0)
        return false;
    char * pattern = npatternArray;
    char * replacement = nreplacementArray;
    anihilatedGuards = 0;
    while (countdown)
        {
        if (*pattern == equal) // =
            {
            // scan replacement until corresponding equal byte found.
            while (*replacement == unequal)
                {
                ++replacement;
                }
            // Corresponding equal byte found in replacement.
            assert(*replacement);
            if (pattern[1] == unequal
                || replacement[1] == unequal
                || (pattern > npatternArray && pattern[-1] == unequal)
                || (replacement > nreplacementArray && replacement[-1] == unequal)
                || (!options->suffixOnly()
                // Start and end are also transition points:
                // (This gives a markedly better result, at least for Dutch)
                && (pattern[1] == '\0'
                || replacement[1] == '\0'
                || pattern == npatternArray
                || replacement == nreplacementArray
                )
                )
                /*
                // 20101206:
                || (  suffixonly
                // And is also transition point:
                && pattern[1] == '\0'
                && replacement[1] == '\0'
                )
                // :20101206
                */
                )
                {
                if (countdown == 1)
                    {
                    *pattern = unequal; // #
                    *replacement = unequal; // Replace the =s with #s and return
                    /* Ways to get rid of a guard
                        =#  -> ##           (1)
                        #=  -> ##           (2)
                        #=# -> ###          (3)
                        ==# -> =##          (4)
                        #== -> ##=          (5)
                        #==# -> ##=#        (6)
                        #==# -> #=##        (7)
                        */
                    if ((pattern == npatternArray // 1
                        || pattern[-1] == unequal   // 2, 3
                        )
                        && (!pattern[1]              // 2
                        || pattern[1] == unequal    // 1, 3
                        )
                        )
                        anihilatedGuards = 1;
                    else if (pattern[1] == equal       // 5, 6
                             && (!pattern[2]            // 5
                             || pattern[2] == unequal  // 6
                             )
                             )
                             anihilatedGuards = 1;
                    else if (pattern > npatternArray
                             && pattern[-1] == equal              // 4, 7
                             && (pattern == npatternArray + 1   // 4
                             || pattern[-2] == unequal         // 7
                             )
                             )
                             anihilatedGuards = 1;
                    return true;
                    }
                --countdown;
                }
            ++replacement; // Go past this = in replacement.
            }
        ++pattern;
        if (!*pattern)
            {
            return false; // signal to caller that countdown has reached maxvalue.
            }
        }
    return false;
    }

