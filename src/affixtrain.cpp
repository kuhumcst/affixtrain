/*
AFFIXTRAIN - supervised learning of affix rules for CSTLEMMA

Copyright (C) 2014  Center for Sprogteknologi, University of Copenhagen

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

#define VERSION "1.73"

#include "flexcombi.h"
#include "optionaff.h"
#include "applyaffrules.h"
#include "graph.h"
#include "comp.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "utf8func.h"

static const char * nil = "";
static const char Start[2] = { START, 0 };
static const char End[2] = { END, 0 };
static const int inil[1] = { 0 };
static const int iStart[2] = { START, 0 };
static const int iEnd[2] = { END, 0 };
static const char StartAny[3] = { START, ANY, 0 };
static const char AnyEnd[3] = { ANY, END, 0 };
static const char StartAnyEnd[4] = { START, ANY, END, 0 };
static const char * SCUT = ".cutoff";
#define CHECK(a)

const char * tempDir(const char * filename, optionStruct * options)
    {
    static char * fullname = NULL;
    size_t length;
    int written;
    CHECK("AglobTempDir");
    if (options->tempDir())
        {
        if (fullname)
            delete[] fullname;
        length = strlen(options->tempDir()) + strlen(filename) + 2;
        fullname = new char[length];
        written = sprintf(fullname, "%s%c%s", options->tempDir(), DIRSEP, filename);
        if ((size_t)written >= length)
            {
            printf("tempDir: Buffer overrun");
            exit(-1);
            }
        return fullname;
        }
    CHECK("BglobTempDir");
    return filename;
    }

union pointers
    {
    unsigned char * uchars;
    char * chars;
    const char * cchars;
    const short int * shorts;
    const long int * longs;
    };

class tagClass
    {
    public:
        char * name;
        tagClass * next;
        tagClass(const char * name, ptrdiff_t length, tagClass * next)
            {
            this->name = new char[length + 1];
            strncpy(this->name, name, length);
            this->name[length] = '\0';
            this->next = next;
            }
        ~tagClass()
            {
            delete[] name;
            delete next;
            }
        bool has(const char * str, ptrdiff_t length)
            {
            return ((length == (ptrdiff_t)strlen(name))
                    && !strncmp(str, name, length)
                    )
                    || (next && next->has(str, length));
            }
    };


class shortRulePair;
class fullRulePair : public rulePair
    {
    private:
        char Pattern[1000];
        char Replacement[1000];
    protected:
        virtual char * itsPattern(){ return Pattern; }
        virtual char * itsReplacement(){ return Replacement; }
    public:
        fullRulePair(const char * pat, const char * rep);
        fullRulePair(shortRulePair * Rule);
        virtual ~fullRulePair(){ --FullRulePairCount; }
        void copy(fullRulePair * Rule);
        const char * pattern()
            {
            return Pattern;
            }
        const char * replacement()
            {
            return Replacement;
            }
    };

FILE * fopenOrExit(const char * name, const char * mode, const char * descriptionOfFile) // must not contain format specifiers
    {
    assert(name);
    assert(mode);
    assert(descriptionOfFile);
    FILE * fp = fopen(name, mode);
    if (!fp)
        {
        fprintf
            (stderr
            , "%s: Cannot open file %s for %s (%s). Exiting.\n"
            , descriptionOfFile
            , name
            , strchr(mode, 'r') ? "reading" : strchr(mode, 'a') ? "appending" : "writing"
            , mode
            );
        exit(-1);
        }
    return fp;
    }

fullRulePair::fullRulePair(const char * pat, const char * rep)
    {
    strcpy(Pattern, pat);
    strcpy(Replacement, rep);
    ++FullRulePairCount;
    }

void fullRulePair::copy(fullRulePair * Rule)
    {
    strcpy(Pattern, Rule->Pattern);
    strcpy(Replacement, Rule->Replacement);
    }

class ruleTemplate;
class shortRulePair
    {
    private:
        char patternArray[1000];
        char replacementArray[1000];
        void trim();
    public:
        const char * itsPatternArray(){ return patternArray; }
        const char * itsReplacementArray(){ return replacementArray; }
        void copy(shortRulePair * Rule)
            {
            strcpy(patternArray, Rule->patternArray);
            strcpy(replacementArray, Rule->replacementArray);
            }
        shortRulePair(trainingPair * trainingpair, ruleTemplate * Template);
        bool checkRule(/*ruleTemplate * Template,*/trainingPair * trainingpair, rulePair * parentPat, optionStruct * options);
        ~shortRulePair(){ --ShortRulePairCount; }
    };

class ruleTemplate
    {
    protected:
        char npatternArray[1000];
        char nreplacementArray[1000];
    public:
        const char * itsNPatternArray(){ return npatternArray; }
        const char * itsNReplacementArray(){ return nreplacementArray; }
        void copy(ruleTemplate * Template)
            {
            strcpy(npatternArray, Template->npatternArray);
            strcpy(nreplacementArray, Template->nreplacementArray);
            }
        bool makebigger(int countdown, int & anihilatedGuards, optionStruct * options);
        ruleTemplate(){ ++RuleTemplateCount; }
        ~ruleTemplate(){ --RuleTemplateCount; }
    };

static const int unequal = '#';
static const int equal = '=';

void suffix(char * mask)
    {
    char * m;
    for (m = mask; *m && *m == equal; ++m)
        ;
    for (; *m; ++m)
        *m = unequal;
    }




class similData : public ruleTemplate
    {
    private:
        char pattern[1000];
        char replacement[1000];
        char * psimilar;
        char * ppattern;
        char * preplacement;
        char * nppattern;
        char * npreplacement;
    public:
        similData(char * similarArray)
            {
            psimilar = similarArray;
            ppattern = pattern;
            preplacement = replacement;
            nppattern = npatternArray;
            npreplacement = nreplacementArray;
            }
        ptrdiff_t simil(
            const char * const s1,
            const char * const n1,
            const char * const s2,
            const char * const n2,
            const char * const start,
            const char * const end);
        ptrdiff_t isimil(
            const int * const s1,
            const int * const n1,
            const int * const s2,
            const int * const n2,
            const int * const start,
            const int * const end);
        void match(trainingPair * trainingpair)
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
        bool mergeTemplates(const char * predef, optionStruct * options)
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
        void print(FILE *fp)
            {
            fprintf(fp, "\npattern:\n%s\nreplacement:\n%s\n", pattern, replacement);
            fprintf(fp, "npattern:\n%s\nnreplacement:\n%s\n", npatternArray, nreplacementArray);
            }
    };


const char * find
(const char * string
, int c
, const char * end
)
    {
    while (string < end && *string != c)
        {
        ++string;
        }
    return (string < end) ? string : NULL;
    }

struct aFile
    {
    union pointers file;
    union pointers * Lines;
    int filesize;
    int lines;
    long size;
    const char * eob;
    char * fname;
    aFile(const char * fname, optionStruct * options) :Lines(NULL), filesize(0), lines(0), eob(NULL)
        {
        assert(fname);

        FILE * fp = fopenOrExit(fname, "rb", "Input file");
        if (options->verbose())
            printf("reading file %s ...", fname);
        CHECK("CglobTempDir");
        this->fname = new char[strlen(fname) + 1];
        strcpy(this->fname, fname);
        file.chars = NULL;
        size = 0;
        while (getc(fp) != EOF)
            {
            ++size;
            }
        file.chars = new char[size + 1];
        rewind(fp);
        size = 0;
        int kar;
        while ((kar = getc(fp)) != EOF)
            {
            file.chars[size++] = (char)kar;
            }
        file.chars[size] = '\0';
        fclose(fp);

        eob = file.chars + size;

        const char * buf;
        int line;
        const char * nl;
        const char * tab;
        for (buf = file.chars, line = 0; buf < eob;)
            {
            for (;;)
                {
                nl = find(buf, '\n', eob);
                if (nl)
                    {
                    tab = find(buf, '\t', nl);
                    if (tab)
                        ++line;
                    else
                        {
                        buf = nl + 1;
                        continue;
                        }
                    }
                break;
                }
            if(nl)
                {
                buf = nl + 1;
                }
            else
                {
                tab = find(buf, '\t', eob);
                if (tab)
                    ++line;
                buf = eob;
                }
            }
        if (options->verbose())
            {
            printf("line %d\n", line);
            }

        lines = line;
        Lines = new pointers[lines + 1];

        for (buf = file.cchars, line = 0; buf < eob;)
            {
            for (;;)
                {
                nl = find(buf, '\n', eob);
                if (nl)
                    {
                    tab = find(buf, '\t', nl);
                    if (tab)
                        Lines[line++].cchars = buf;
                    else
                        {
                        buf = nl + 1;
                        continue;
                        }
                    }
                break;
                }
            if(nl)
                {
                buf = nl + 1;
                }
            else
                {
                tab = find(buf, '\t', eob);
                if (tab)
                    Lines[line++].cchars = buf;
                buf = eob;
                }
            }
        if (options->verbose())
            {
            printf("line %d read\n", line);
            }
        if (options->verbose())
            {
            printf("%s:", fname);
            }
        }
    ~aFile()
        {
        delete[] file.chars;
        delete[] Lines;
        delete[] this->fname;
        }
    };

static const int unequalBit = 1 << 0;
static const int equalBit = 1 << 1;
static const int transition = unequalBit | equalBit;

static bool glob_wildcards = false;

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
        /*The strings s1 and s2 are completely different*/
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

fullRulePair::fullRulePair(shortRulePair * Rule)
    {
    CHECK("EglobTempDir");
    char * ppat = Pattern;
    char * prep = Replacement;

    size_t patlen = strlen((const char*)Rule->itsPatternArray());
    size_t replen = strlen((const char*)Rule->itsReplacementArray());
    if (*Rule->itsPatternArray() != START)
        {
        if (*Rule->itsReplacementArray() != START)
            {
            strcpy(Pattern, StartAny/*"^*"*/);
            strcpy(Replacement, StartAny/*"^*"*/);
            ppat += 2;
            prep += 2;
            }
        else
            {
            strcpy(Pattern, Start/*"^"*/);
            ++ppat;
            }
        }
    else if (*Rule->itsReplacementArray() != START)
        {
        strcpy(Replacement, Start/*"^"*/);
        ++prep;
        }
    strcpy(ppat, Rule->itsPatternArray());
    strcpy(prep, Rule->itsReplacementArray());
    ppat += patlen - 1;
    prep += replen - 1;
    if (*ppat == ANY && *prep == ANY)
        {
        strcpy(++ppat, End/*"$"*/);
        strcpy(++prep, End/*"$"*/);
        }
    else
        {
        if (*ppat != END)
            strcpy(++ppat, AnyEnd/*"*$"*/);
        if (*prep != END)
            strcpy(++prep, AnyEnd/*"*$"*/);
        }
    ++FullRulePairCount;
    }


bool rulePair::apply(trainingPair * trainingpair, size_t lemmalength, char * lemma, char * mask, optionStruct * options)
    {
    CHECK("FglobTempDir");
    char wrd[100];
    if (100 <= sprintf(wrd, "%c%.*s%c", START, (int)(trainingpair->itsWordlength()), trainingpair->itsWord(), END))
        {
        printf("rulePair::apply too small buffer");
        exit(1);
        }
    char * p = itsPattern();
    char * r = itsReplacement();
    if (!*p) // root
        {
        static char lStartAnyEnd[4] = { START, ANY, END, 0 };
        p = lStartAnyEnd/*"^*$"*/;
        r = p;
        }
    char * w = wrd;
    char * d = lemma;
    char * last = lemma + lemmalength - 7;
    char * m = mask;
    size_t inc;
    int P = UTF8char(p, UTF8);
    int W = UTF8char(w, UTF8);
    int R = UTF8char(r, UTF8);
    //    while(*p && *r)
    while (P && R)
        {
        //if(*p == *w)
        if (P == W)
            {
            *m = unequal;
            while (R && R != ANY && d < last)
                //while(*r && *r != ANY)
                {
                inc = copyUTF8char(r, d);
                r += inc;
                d += inc;
                //*d++ = *r++;
                R = UTF8char(r, UTF8);
                }
            do
                {
                *++m = unequal;
                p += skipUTF8char(p);
                w += skipUTF8char(w);
                P = UTF8char(p, UTF8);
                W = UTF8char(w, UTF8);
                //++p;
                //++w;
                } while (P && P != ANY && P == W);
                if (P != R)
                    {
                    if (options->suffixOnly())
                        {
                        suffix(mask);
                        }
                    return false;
                    }
            }
        else if (R == ANY)
            {
            p += skipUTF8char(p);
            r += skipUTF8char(r);
            P = UTF8char(p, UTF8);
            R = UTF8char(r, UTF8);
            char * ep = strchr(p, ANY);
            if (ep)
                *ep = '\0';
            char * sub = strstr(w, p);
            if (sub)
                {
                while (w < sub && d < last)
                    {
                    *m++ = equal;
                    //*d++ = *w++;
                    inc = copyUTF8char(w, d);
                    w += inc;
                    d += inc;
                    }
                W = UTF8char(w, UTF8);
                }
            else
                {
                if (ep)
                    *ep = ANY;
                d = lemma;
                break;
                }
            if (ep)
                *ep = ANY;
            }
        else
            {
            d = lemma;
            break;
            }
        }
    if (P || R)
        {
        if (options->suffixOnly())
            {
            suffix(mask);
            }
        return false;
        }
    *--m = '\0';
    for (m = mask; *m; ++m)
        m[0] = m[1];
    *d = '\0';
    d = lemma;
    char * oldd = d;
    inc = skipUTF8char(d);
    d += inc;
    while (*d)
        {
        *oldd++ = *d++;
        }
    *oldd = '\0';
    if (oldd > lemma + 1)
        oldd[-1] = '\0';
    if (options->suffixOnly())
        {
        suffix(mask);
        }
    return true;
    }

shortRulePair::shortRulePair(trainingPair * trainingpair, ruleTemplate * Template)
    {
    CHECK("GglobTempDir");
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

static void shiftleft(char * a)
    {
    CHECK("HglobTempDir");
    while (*a)
        {
        a[0] = a[1];
        ++a;
        }
    }

static void shiftright(char * A)
    {
    CHECK("IglobTempDir");
    char * a = A + strlen(A);
    a[1] = '\0';
    while (a > A)
        {
        a[0] = a[-1];
        --a;
        }
    }

static int cnt(char * s, int k)
    {
    CHECK("JglobTempDir");
    int ret = 0;
    while (*s)
        {
        if (*s == k)
            ++ret;
        ++s;
        }
    return ret;
    }

void shortRulePair::trim()
    {
    CHECK("KglobTempDir");
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
    CHECK("MglobTempDir");
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


bool shortRulePair::checkRule(/*ruleTemplate * Template,*/trainingPair * trainingpair, rulePair * parentPat, optionStruct * options)
    {
    CHECK("NglobTempDir");
    char Lemma[100] = "";
    char Mask[100] = "";
    fullRulePair FullRule(this);
    bool ret;
    /*
    if(parentPat)
    fprintf(flog,"parent:pat %s rep %s\n",parentPat->pattern(),parentPat->replacement());
    fprintf(flog,"this  :pat %s rep %s\n",FullRule.pattern(),FullRule.replacement());
    fprintf(flog,"trainingpair:");
    trainingpair->print(flog);
    fprintf(flog,"\n");
    */
    if ((!parentPat || FullRule.dif(parentPat) == dif_bigger)
        && FullRule.apply(trainingpair, sizeof(Lemma), Lemma, Mask, options)
        )
        {
        ret = trainingpair->isCorrect(Lemma);
        /*
        fprintf(flog,"ret:%s\n",ret ? "true" : "false");
        */
        }
    else
        {
        /*
        fprintf(flog,"ret is false\n");
        */
        ret = false;
        }
    return ret;
    }

static int storeRule(hash * Hash, shortRulePair * Rule, vertex *& V)
    {
    CHECK("OglobTempDir");
    fullRulePair FullRule(Rule);
    bool New;
    V = Hash->getVertex(&FullRule, New);
    /*
    V->print1(flog);
    fprintf(flog,"\n");
    */
    if (New)
        {
        return 1;
        }
    else
        {
        return 0;
        }
    }

int trainingPair::makeCorrectRules(hash * Hash, ruleTemplate * Template, const char * similar, vertex * parent, int mlow, int recurse, optionStruct * options)
    {
    CHECK("PglobTempDir");
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
        if (Rule.checkRule(/*&locTemplate,*/this, parent, options))
            {
            different = false;
            vertex * e;
            storeRule(Hash, &Rule, e);
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

int trainingPair::makeRuleEx(hash * Hash, vertex * parent, bool alreadyRight, optionStruct * options)
    {
    CHECK("QglobTempDir");
    /*
    this->print(flog);
    */
    char similarArray[1000];
    similData SimilData(similarArray);
    SimilData.match(this);
    /*
    SimilData.print(flog);
    */
    const char * predef = getMask();
    /*
    fprintf(flog,"predef %s\n",predef);
    */
    SimilData.mergeTemplates(predef, options);
    /*
    fprintf(flog,"mergeTemplates DONE:\n");
    SimilData.print(flog);
    */
    shortRulePair Rule(this, &SimilData);
    if (Rule.checkRule(/*&SimilData,*/this, parent, options))
        {
        vertex * e;
        int ret = storeRule(Hash, &Rule, e);
        /** /
        fprintf(flog,"makeRuleEx:");
        e->print1(flog);
        fprintf(flog,"\n");
        */
        return ret;
        }
    int nr = makeCorrectRules(Hash, &SimilData, similarArray, parent, 1, options->maxRecursionDepthRuleCreation(), options);
    if (nr == 0 && !alreadyRight)
        {
        if (options->verbose())
            { // This is by design. Increasing RECURSE will eventually help.
            fprintf(stderr, "Error (makeRuleEx): Cannot construct rule for trainingpair ");
            this->print(stderr);
            fprintf(stderr, " Based on parent ");
            parent->print1(stderr);
            fprintf(stderr, "\n");
            }
        }
    return nr;
    }

static bool isSpace(int a)
    {
    CHECK("RglobTempDir");
    switch (a)
        {
        case ' ':
        case '\t':
        case '\r':
        case '\n':
        return true;
        default:
        return false;
        }
    }
static tagClass * collectTags(optionStruct * options)
// "123456" means Word, Lemma, Wordfreq, Lemmafreq, Wordclass, Lemmaclass
    {
    struct aFile afile(options->wordList(), options);

    CHECK("SglobTempDir");
    tagClass * Tags = NULL;
    int line;
    pointers * Lines = afile.Lines;
    if (options->verbose())
        {
        printf("Checking %d lines in %s for tags %s\n", afile.lines, afile.fname, options->columns());
        }
    for (line = 0; line < afile.lines; ++line)
        {
        const char * cols[18];
        ptrdiff_t lengths[18];
        const char * limit = (line < afile.lines - 1) ? Lines[line + 1].cchars : afile.eob;

        const char * q = Lines[line].cchars;
        cols[0] = q;
        unsigned int ii = 0;
        while ((q = find(q, '\t', limit)) != NULL
               && (ii < sizeof(cols) / sizeof(cols[0]) - 1)
               )
            {
            lengths[ii] = q - cols[ii];
            cols[++ii] = ++q;
            }
        lengths[ii] = limit - cols[ii] - 1;
        unsigned int maxii = ++ii;
        cols[maxii] = q ? q : limit;
        const char * column;
        for (column = options->columns(), ii = 0
             ; *column && (ii < maxii)
             ; ++column, ++ii
             )
            {
            switch (*column)
                {
                case '3':
                case 'T':
                case 't':
                if (!Tags || !(Tags->has(cols[ii], lengths[ii])))
                    {
                    Tags = new tagClass(cols[ii], lengths[ii], Tags);
                    }
                break;
                default:
                ;
                }
            }
        }
    return Tags;
    }

static int compare(const void * arg1, const void * arg2)
    {
    CHECK("UglobTempDir");
    const trainingPair * A = *(const trainingPair * const *)arg1;
    const trainingPair * B = *(const trainingPair * const *)arg2;
    int ret = A->cmpWord(B);
    if (!ret)
        {
        ret = A->cmpLemma(B);
#if WORDCLASS
        if(!ret)
            {
            ret = A->cmpWordClass(B);
            }
#endif
#if LEMMACLASS
        if(!ret)
            {
            ret = A->cmpLemmaClass(B);
            }
#endif
        }
    return ret;
    }

static int markAmbiguous(int allPairs, trainingPair * TrainingPair, FILE * famb, FILE * fallFile, optionStruct * options)
    {
    if (options->verbose())
        printf("markAmbiguous\n");
    trainingPair ** pTrainingPair = new trainingPair *[allPairs];
    int j;
    for (j = 0; j < allPairs; ++j)
        pTrainingPair[j] = TrainingPair + j;
    qsort((void *)pTrainingPair, allPairs, sizeof(trainingPair *), compare);

    for (j = 0; j < allPairs; ++j)
        // Mark words containing space for skipping
        {
        size_t nn = pTrainingPair[j]->itsWordlength();
        const char * s;
        // do not use words containing space
        for (s = pTrainingPair[j]->itsWord(); nn-- > 0;)
            {
            if (s[nn] == ' ')
                {
                pTrainingPair[j]->set(b_skip);
                }
            }
        }
    int n = 0;



    for (j = 0; j < allPairs; ++j)
        {
        if (!pTrainingPair[j]->isset(b_skip))
            {
            int k;
            for (k = j + 1
                 ;    k < allPairs
                 //&& !pTrainingPair[j]->isset(b_skip)
                 && !pTrainingPair[j]->cmpWord(pTrainingPair[k])
#if WORDCLASS
                 && !pTrainingPair[j]->cmpWordClass(pTrainingPair[k])
#endif
                 ; ++k
                 )
                {
                if (!pTrainingPair[j]->cmpLemma(pTrainingPair[k]) /*|| pTrainingPair[j]->cmpLemmaClass(pTrainingPair[k])*/) //TODO Tags!
                    pTrainingPair[k]->set(b_doublet | b_skip);
                }
            }
        }



    for (j = 0; j < allPairs;// - 1;
#if AMBIGUOUS
#else
         ++j
#endif
         )
        {
        if (pTrainingPair[j]->isset(b_skip))
            ++j;
        else
            {
            int k;
#if AMBIGUOUS
            trainingPair * pAlt = pTrainingPair[j];
            if (famb)
                {
                //if(!pTrainingPair[j]->isset(b_skip))
                pTrainingPair[j]->fprint(famb);
                }
#endif
            for (k = j + 1
                 ;    k < allPairs
                 //&& !pTrainingPair[j]->isset(b_skip)
                 && !pTrainingPair[j]->cmpWord(pTrainingPair[k])
#if WORDCLASS
                 && !pTrainingPair[j]->cmpWordClass(pTrainingPair[k])
#endif
                 ; ++k
                 )
                {
                if (pTrainingPair[j]->cmpLemma(pTrainingPair[k]) /*|| pTrainingPair[j]->cmpLemmaClass(pTrainingPair[k])*/) //TODO Tags!
                    {
#if AMBIGUOUS
                    if (!pTrainingPair[k]->isset(b_skip))
                        {
                        pTrainingPair[k]->set(b_ambiguous);
                        pAlt->Alt = pTrainingPair[k];
                        pAlt = pTrainingPair[k];
                        if (famb)
                            {
                            fprintf(famb, "HOMOGRAPH:");
                            pTrainingPair[k]->fprint(famb);
                            }
                        }
#else
#if LEMMAINL
                    switch(pTrainingPair[j]->cmpFreq(pTrainingPair[k]))
                        {
                        case 1:
                        pTrainingPair[k]->set(b_skip);
                        break;
                        case -1:
                        pTrainingPair[j]->set(b_skip);
                        break;
                        default:
#else
                            {
#endif
                            /* Take the first one */
                            if(!pTrainingPair[k]->isset(b_ambiguous|b_skip))
                                {
                                pTrainingPair[k]->set(b_ambiguous);
                                if(famb)
                                    pTrainingPair[k]->fprint(famb);
                                ++n;
                                }
                            }
#endif
                    }
                else
                    pTrainingPair[k]->set(b_doublet | b_skip);
                }
#if AMBIGUOUS
            // close the ring of homographs
            if (pAlt != pTrainingPair[j])
                {
                pTrainingPair[j]->set(b_ambiguous); // ALL members in the ring are marked b_ambiguous!
                pAlt->Alt = pTrainingPair[j];
                }
            j = k;
#else
            //        if(allFile && j < allPairs)
            //          pTrainingPair[j]->fprintAll(allFile);
#endif
            }
        }

    for (j = 0; j < allPairs; ++j)
        {
        if (!pTrainingPair[j]->isset(b_skip))
            pTrainingPair[j]->fprintAll(fallFile);
        }

    delete[] pTrainingPair; // Bart 20081008
    if (options->verbose())
        printf("markAmbiguous DONE\n");
    return n;
    }

#if PESSIMISTIC
static int compare2(const void * arg1, const void * arg2)
    {
    const trainingPair * A = *(const trainingPair * const *)arg1;
    const trainingPair * B = *(const trainingPair * const *)arg2;
    int ret = A->cmpLemma(B);
    if(!ret)
        {
        ret = A->cmpLemmaClass(B);
        if(!ret)
            {
            ret = A->cmpWord(B);
            }
        }
    return ret;
    }

static int markParadigms(int allPairs,trainingPair * TrainingPair,FILE * fparadigms)
    {
    //    FILE * allFile = fopenwb("allFile.txt");
    trainingPair ** pTrainingPair = new trainingPair * [allPairs];
    int j;
    for(j = 0;j < allPairs;++j)
        pTrainingPair[j] = TrainingPair + j;
    qsort((void *)pTrainingPair, allPairs, sizeof(trainingPair *), compare2);

    int n = 0;
    for(j = 0;j < allPairs/* - 1*/;)
        {
        int k;
        trainingPair * pAlt = pTrainingPair[j];
        if(fparadigms)
            {
            if(!pTrainingPair[j]->isset(b_skip))
                pTrainingPair[j]->fprint(fparadigms);
            }
        for(k = j+1
            ;    k < allPairs
            && !pTrainingPair[j]->isset(b_skip)
            && !pTrainingPair[j]->cmpLemma(pTrainingPair[k])
            ; ++k
            )
            {
            if(!pTrainingPair[k]->isset(b_skip))
                {
                pAlt->AltLemma = pTrainingPair[k];
                pAlt = pTrainingPair[k];
                if(fparadigms)
                    {
                    fprintf(fparadigms,"PARADIGM:");
                    pTrainingPair[k]->fprint(fparadigms);
                    }
                }
            }
        if(pAlt != pTrainingPair[j])
            {
            pAlt->AltLemma = pTrainingPair[j];
            }
        j = k;
        }
    delete [] pTrainingPair; // Bart 20081008
    return n;
    }
#endif

// Mark all pairs that are going to be used for testing after the training has completed. (If there is such testing.)
static void markTest(int allPairs, trainingPair * TrainingPair, FILE * ftrain, FILE * ftestp, FILE * ftests, optionStruct * options)
    {
    if (options->verbose())
        printf("markTest\n");
    int pairs = 0;
    for (pairs = 0; pairs < allPairs; ++pairs)
        {
        if (!TrainingPair[pairs].isset(
#if !AMBIGUOUS
            b_ambiguous|
#endif
            b_doublet | b_skip))
            {
            if (rand() % 100 <= options->percentageTestPairs())
                {
                TrainingPair[pairs].set(b_test);
                TrainingPair[pairs].print(ftestp);
                fprintf(ftestp, "\n");
                fprintf(ftests, "%.*s\n", (int)(TrainingPair[pairs].itsWordlength()), TrainingPair[pairs].itsWord());
                }
            else
                {
                TrainingPair[pairs].print(ftrain);
                fprintf(ftrain, "\n");
                }
            }
        }
    if (options->verbose())
        printf("markTest DONE\n");
    }



static void writeAvailableTrainingData(int allPairs, trainingPair * TrainingPair, FILE * ftrainavailable
#if WRITEINVERTED
                                       ,FILE * fpinverted
#endif
                                       , optionStruct * options
                                       )
    {
    if (options->verbose())
        printf("writeAvailableTrainingData\n");
    int pairs = 0;
    for (pairs = 0; pairs < allPairs; ++pairs)
        {
#if AMBIGUOUS
        if (!TrainingPair[pairs].isset(b_doublet | b_skip))
#else
        if(!TrainingPair[pairs].isset(b_ambiguous|b_doublet|b_skip))
#endif
            {
            TrainingPair[pairs].fprintTraining(ftrainavailable
#if WRITEINVERTED
                                               ,fpinverted
#endif
                                               );
            }
        }
    if (options->verbose())
        printf("writeAvailableTrainingData DONE\n");
    }

int trainingPair::makeNextTrainingSet
(int allPairs
, trainingPair * TrainingPair
, FILE * train
, FILE * done
, FILE * combined
, FILE * disamb
, optionStruct * options
)
    {
    if (options->verbose())
        printf("makeNextTrainingSet\n");

    int pairs = 0;
    int donepairs = 0;
    for (pairs = 0; pairs < allPairs; ++pairs)
        {
#if AMBIGUOUS
        if (!TrainingPair[pairs].isset(b_doublet | b_skip))
#else
        if(!TrainingPair[pairs].isset(b_ambiguous|b_doublet|b_skip))
#endif
            {
            if (!TrainingPair[pairs].isset(b_ambiguous))
                {
                TrainingPair[pairs].printMore(train);
                TrainingPair[pairs].printMore(disamb);
                fprintf(combined, "\t");
                TrainingPair[pairs].printMore(combined);
                }
            else if (TrainingPair[pairs].isset(b_bench))
                {
                TrainingPair[pairs].printMore(train);
                fprintf(combined, "\t");
                TrainingPair[pairs].printMore(combined);
                }
            else
                {
                donepairs++;
                TrainingPair[pairs].printMore(done);
                TrainingPair[pairs].printMore(disamb);
                fprintf(combined, "*\t");
                TrainingPair[pairs].printMore(combined);
                }
            }
        }
    if (options->verbose())
        printf("makeNextTrainingSet DONE\n");

    return donepairs;
    }

void trainingPair::makeChains(int allPairs, trainingPair * TrainingPair, trainingPair ** train, trainingPair ** test, optionStruct * options)
    {
    if (options->verbose())
        printf("makeChains\n");

    int pairs = 0;
    trainingPair ** ptrain = train;
    trainingPair ** ptest = test;
    for (pairs = 0; pairs < allPairs; ++pairs)
        {
#if AMBIGUOUS
        if (!TrainingPair[pairs].isset(b_doublet | b_skip))
#else
        if(!TrainingPair[pairs].isset(b_ambiguous|b_doublet|b_skip))
#endif
            {
            if (TrainingPair[pairs].isset(b_test))
                {
                *ptest = TrainingPair + pairs;
                ptest = &TrainingPair[pairs].Next;
                }
            else
                {
                *ptrain = TrainingPair + pairs;
                ptrain = &TrainingPair[pairs].Next;
                }
            }
        }
    if (options->verbose())
        printf("makeChains DONE\n");
    }

#if DOTEST
// Testing. (Not called during training!)
static void showResults(trainingPair * TrainingPair,int & wrong,int & right,int & both,FILE * fr)
    {
    for(trainingPair * tp = TrainingPair;tp;tp = tp->next())
        {
        if(!tp->isset(b_ambiguous|b_doublet|b_skip)) // b_ambiguous is ok if AMBIGUOUS == 1, or what?
            {
                    {
                    if(tp->isset(b_wrong))
                        {
                        if(tp->isset(b_ok))
                            {
                            if(fr)
                                {
                                fprintf(fr,"/ ");
                                tp->fprint(fr);
                                }
                            both++;
                            }
                        else
                            {
                            if(fr)
                                {
                                fprintf(fr,"- ");
                                tp->fprint(fr);
                                }
                            wrong++;
                            }
                        }
                    else if(tp->isset(b_ok))
                        {
                        if(fr)
                            {
                            fprintf(fr,"+ ");
                            tp->fprint(fr);
                            }
                        right++;
                        }
                    tp->unset(b_ok|b_wrong);
                    }
            }
        }
    }

// Testing. (Not called during training!)
static void lemmatise(node * tree,trainingPair * TrainingPair)
    {
    /*
    printf("lemmatise\n");
    */
    for(trainingPair * tp = TrainingPair;tp;tp = tp->next())
        {
        if(!tp->isset(b_ambiguous|b_doublet|b_skip)) // Skip all homographs, even if AMBIGUOUS == 1 !
            {
            //tp->print(stdout);printf("\n");
            tree->lemmatise(tp);
            }
        }
    /*
    printf("lemmatise DONE\n");
    */
    }
static FILE * fresults = NULL;
static double maxCorrectness = 0.0;
#endif

static trainingPair * globTrainingPair;
//static int globlines;

static trainingPair * readTrainingPairs(aFile & afile, int & pairs, const char * columns, const char * tag, optionStruct * options)
    {
    if (options->verbose())
        printf("readTrainingPairs\n");

    //    globlines = afile.lines;
    trainingPair * TrainingPair = new trainingPair[afile.lines];
    globTrainingPair = TrainingPair;
    // "123456" means Word, Lemma, Wordfreq, Lemmafreq, Wordclass, Lemmaclass
    CHECK("TglobTempDir");
    pairs = 0;
    int line;
    int taglength = tag ? strlen(tag) : 0;
    pointers * Lines = afile.Lines;
    if (options->verbose())
        {
        printf("readLines with tag %s\n", (tag && *tag) ? tag : "_ (UNDEFINED)");
        }
    for (line = 0; line < afile.lines; ++line)
        {
        const char * cols[18];
        const char * Word = NULL;
        ptrdiff_t lengths[18];
        size_t wordlength = 0;
#if WORDCLASS
        const char * WordClass = NULL; // unused
        size_t wordclasslength = 0; // unused
#endif
        const char * LemmaHead = NULL;
        size_t lemmalength = 0;
#if LEMMACLASS
        const char * LemmaClass = NULL;
        size_t lemmaclasslength = 0;
#endif
#if LEMMAINL
        long Inl = 0;
        long Lemma_Inl = 0;
#endif
        const char * limit = (line < afile.lines - 1) ? Lines[line + 1].cchars : afile.eob;

        const char * q = Lines[line].cchars;
        cols[0] = q;
        unsigned int ii = 0;
        while ((q = find(q, '\t', limit)) != NULL
               && (ii < sizeof(cols) / sizeof(cols[0]) - 1)
               )
            {
            lengths[ii] = q - cols[ii];
            cols[++ii] = ++q;
            }
        lengths[ii] = limit - cols[ii] - 1;
        unsigned int maxii = ++ii;
        cols[maxii] = q ? q : limit;
        const char * column;
        bool doUse = ((tag == NULL) || !*tag);
        for (column = columns, ii = 0
             ;    *column
             && (ii < maxii)
             ; ++column, ++ii
             )
            {
            switch (*column)
                {
                case '1':
                case 'F':
                case 'f':
                case 'W':
                case 'w':
                Word = cols[ii];
                wordlength = lengths[ii];
                break;
                case '2':
                case 'B':
                case 'b':
                case 'L':
                case 'l':
                LemmaHead = cols[ii];
                lemmalength = lengths[ii];
                break;
                case '3':
                case 'T':
                case 't':
                if (tag && (taglength == lengths[ii]) && !strncmp(cols[ii], tag, taglength))
                    {
                    doUse = true;
                    }
                break;
#if LEMMAINL
                case '3':
                Inl = strtol(cols[ii],NULL,10);
                break;
                case '4':
                Lemma_Inl = strtol(cols[ii],NULL,10);
                break;
#endif
#if WORDCLASS
                case '5':
                WordClass = cols[ii];
                wordclasslength = lengths[ii];
                break;
#endif
#if LEMMACLASS
                case '6':
                LemmaClass = cols[ii];
                lemmaclasslength = lengths[ii];
                break;
#endif
                default:
                ;
                }
            }
        if (Word && doUse)
            {
            doUse = ((tag == NULL) || !*tag);
            while (wordlength > 0 && isSpace(Word[wordlength - 1]))
                --wordlength;
            while (lemmalength > 0 && isSpace(LemmaHead[lemmalength - 1]))
                --lemmalength;
#if WORDCLASS
            while(wordclasslength > 0 && isSpace(WordClass[wordclasslength- 1]))
                --wordclasslength;
#endif
#if LEMMACLASS
            while(lemmaclasslength > 0 && isSpace(LemmaClass[lemmaclasslength - 1]))
                --lemmaclasslength;
#endif
            if (wordlength > 0 && lemmalength > 0)
                {
                TrainingPair[pairs].init
                    (Word
                    , wordlength
                    , LemmaHead
                    , lemmalength
#if WORDCLASS
                    ,WordClass
                    ,wordclasslength
#endif
#if LEMMACLASS
                    ,LemmaClass
                    ,lemmaclasslength
#endif
#if LEMMAINL
                    ,Inl,Lemma_Inl
#endif
                    );
                ++pairs;
                }
            }
        }
    if (options->verbose())
        {
        if (tag && *tag)
            printf("%ld characters and %d lines, %d selected with tag %s       \n", afile.size, afile.lines, pairs, tag);
        else
            printf("%ld characters and %d lines\r", afile.size, afile.lines);
        }
    if (options->verbose())
        printf("readTrainingPairs DONE\n");
    return TrainingPair;
    }

static void markTheAmbiguousPairs(trainingPair * TrainingPair, const char * ext, int pairs, optionStruct * options)
    {
    char filename[256];
    if (256 <= sprintf(filename, "ambiguouspairs_%s.txt", ext))
        {
        printf("markTheAmbiguousPairs: buffer too small");
        exit(-1);
        }
    FILE * famb = fopenOrExit(tempDir(filename, options), "wb", "famb");

    if (256 <= sprintf(filename, "allFile_%s.txt", ext))
        {
        printf("markTheAmbiguousPairs: buffer too small");
        exit(-1);
        }
    FILE * fallFile = fopenOrExit(tempDir(filename, options), "wb", "fallFile");

    /*int ambi =*/ markAmbiguous(pairs, TrainingPair, famb, fallFile, options);
    fclose(famb);
    fclose(fallFile);

#if PESSIMISTIC
    sprintf(filename,"paradigms_%s.txt",ext);
    FILE * fparadigms = fopenOrExit(filename,"wb","fparadigms");
    markParadigms(pairs,TrainingPair,fparadigms);
    if(fparadigms)
        fclose(fparadigms);
#endif
    }

static void writeAllAvailablePairs(trainingPair * TrainingPair, const char * ext, int pairs, optionStruct * options)
    {
    CHECK("cglobTempDir");
    char train[256];
    if (256 <= sprintf(train, "availabletrainingpairs_%s.txt", ext))
        {
        printf("writeAllAvailablePairs: buffer 2 small");
        exit(0);
        }
    FILE * ftrain = fopenOrExit(tempDir(train, options), "wb", "ftrain");

#if WRITEINVERTED
    sprintf(train,"availabletrainingpairs_inverted_%s.txt",ext);
    FILE * fpinverted = fopenOrExit(train,"wb","inverted");
#endif
    writeAvailableTrainingData(pairs, TrainingPair, ftrain
#if WRITEINVERTED
                               ,fpinverted
#endif
                               , options
                               );
    if (ftrain)
        {
        fclose(ftrain);
        ftrain = NULL;
        }
#if WRITEINVERTED
    if(fpinverted)
        {
        fclose(fpinverted);
        ftrain = NULL;
        }
#endif
    }

static void writeAllTestPairs(trainingPair * TrainingPair, const char * ext, int pairs, optionStruct * options)
    {
    CHECK("dglobTempDir");
    char bigbuf[10000];

    sprintf(bigbuf, "trainingpairs_%s.txt", ext);
    FILE * ftrain = fopenOrExit(tempDir(bigbuf, options), "wb", "ftrain");

    sprintf(bigbuf, "testpairs_%s.txt", ext);
    FILE * ftestp = fopenOrExit(tempDir(bigbuf, options), "wb", "ftestp");

    sprintf(bigbuf, "tests_%s.txt", ext);
    FILE * ftests = fopenOrExit(tempDir(bigbuf, options), "wb", "ftests");

    markTest(pairs, TrainingPair, ftrain, ftestp, ftests, options);

    if (options->verbose())
        {
        printf("trainingpairs_%s.txt testpairs_%s.txt tests_%s.txt written\n", ext, ext, ext);
        }
    fclose(ftrain);
    fclose(ftestp);
    fclose(ftests);
    }
/*
static void checkTrainingPairIntegrity()
{
int ii;
for(ii = 0;ii < globlines;++ii)
if(!globTrainingPair[ii].isset(b_ambiguous|b_doublet|b_test|b_skip))
globTrainingPair[ii].checkIntegrity();
}
*/
static void doTheRules(hash * Hash, trainingPair * TrainingPair, node ** top, optionStruct * options)
    {
    CHECK("eglobTempDir");
    //    fullRulePair ROOT("^*$","^*$");
    fullRulePair ROOT(StartAnyEnd, StartAnyEnd/*"^*$","^*$"*/);
    bool New;
    vertex * best = Hash->getVertex(&ROOT, New);

    *top = new node(best);
    trainingPair * Right = NULL;
    (*top)->init(&Right, &TrainingPair, 0/*,0,0*/, options);
    (*top) = (*top)->cleanup(NULL, options);
    }


/* Create binary output */
static void rearrange
(const char * filename
, FILE * folel
#if RULESASTEXT
, FILE * foleltxt
#endif
)
    {
    CHECK("fglobTempDir");
    FILE * fo = fopenOrExit(filename, "wb", "rearrange output");
    long end;
    end = ftell(folel);
    char * buf = new char[end + 1]; // contents of folel, a textual file.
    rewind(folel);
    if (fread(buf, end, 1, folel) != 1)
        return;
    buf[end] = '\0'; // 20140224 new
    unsigned int n = 0;
    int i;
    for (i = 0; i < end; ++i)
        if (buf[i] == '\n')
            ++n;
    char ** pbuf = new char *[n + 1]; // lines
    unsigned int * size = new unsigned int[n + 1]; // binary size taken up by each rule
    unsigned int * cumsize = new unsigned int[n + 1]; // idem, cumulative
    int * ind = new int[n + 1]; // nesting levels, also used to store cumulative sizes
    pbuf[0] = buf + 0;
    n = 0;
    bool doind = true;  // first field on line is an indentation (or nesting
    // level) number
    ind[0] = 0;
    cumsize[0] = 0;
    size[0] = 0;
    for (i = 0; i < end; ++i)
        {
        if (buf[i] == '\n')
            {
            //printf("\n");
            size[n] += sizeof(int) + sizeof(int); // room for index and '\0' byte
            size[n] >>= 2;
            size[n] <<= 2; // rounded to nearest word boundary
            pbuf[++n] = buf + i + 1; // start of next line
            cumsize[n] = size[n - 1];
            cumsize[n] += cumsize[n - 1];
            size[n] = 0;
            doind = true; // read first data on line as nesting level
            ind[n] = 0; // initialize nesting level to 0
            }
        else if (buf[i] == '\t')
            {
            if (doind)
                doind = false;
            else
                size[n]++;
            }
        else if (doind)
            {// read nesting level
            ind[n] *= 10;
            //printf("%c",buf[i]);
            ind[n] += buf[i] - '0';
            }
        else
            size[n]++;
        }
    pbuf[n] = NULL;
    int lev[50];
    unsigned int j;
    for (j = 0; j < sizeof(lev) / sizeof(lev[0]); ++j)
        lev[j] = 0;
    for (j = 0; j < n; ++j)
        {
        int oldj = lev[ind[j]]; // find previous sibling
        if (oldj)
            ind[oldj] = cumsize[j]; // tell previous sibling where its next 
        // sibling is
        for (int k = ind[j] + 1; k < 50 && lev[k]; ++k)
            {
            lev[k] = 0; // forget about previous sibling'c children
            }
        lev[ind[j]] = j; // update who's the current node at this level.
        ind[j] = 0; // the level information is not needed anymore. Space can 
        // be reused to keep the position to the next sibling 
        // (if any).
        }
    char * p = buf;
    for (j = 0; j < n; ++j)
        {
        long pos = ftell(fo);
#ifdef _DEBUG
        long ppos = cumsize[j];
        assert(pos == ppos);
#endif
        if (ind[j] >= pos)
            ind[j] -= pos; // We only need to know how far to jump from here.
        fwrite(ind + j, sizeof(int), 1, fo);

#if RULESASTEXT
        fprintf(foleltxt,"%d",ind[j]);
#endif
        unsigned int written = sizeof(int);

        while (*p && *p != '\t')
            ++p;
        assert(*p);
        ++p;
        while (*p != '\n')
            { // write the pattern and the replacement (which are intertwined)
#if RULESASTEXT
            fputc(*p,foleltxt);
#endif
            fputc(*p++, fo);
            written++;
            }
#if RULESASTEXT
        fputc(*p,foleltxt); 
#endif
        fputc(*p++, fo); // write the end-of-line marker
        written++;
        assert(written <= size[j]);
        assert(written > size[j] - sizeof(int));
        while (written < size[j]) // write until word boundary hit
            {
#if RULESASTEXT
            fputc('\n',foleltxt);
#endif
            fputc('\n', fo);
            written++;
            }
        }
    fclose(fo);
    delete[] buf;
    delete[] size;
    delete[] cumsize;
    delete[] ind;
    }

//int Nnodes = 0;

static bool writeAndTest(node * tree, const char * ext, int threshold, const char * nflexrules, int & Nnodes, double & weight, optionStruct * options)
    {
    CHECK("gglobTempDir");
    char filename[1000];
#if RULESASTEXTINDENTED
    sprintf(filename,"tree_%d%s.txt",threshold,ext);
    FILE * foo = fopenOrExit(filename,"wb","indented rules");

    if(foo)
#endif
        {
#if RULESASTEXTINDENTED
        fprintf(foo,"threshold %d\n",threshold);
        int NnodesR = 0;
        int N = tree->print(foo,0,Nnodes,NnodesR);
        fprintf(foo,"threshold %d:%d words %d nodes %d weight %f nodes with words\n\n",threshold,N,Nnodes,NnodesR,weight);
        fclose(foo);
        sprintf(filename,"rules_%d%s.txt",threshold,ext);
        foo = fopenOrExit(filename,"wb","rules");
#else
        Nnodes = 0;
        Nnodes = tree->count();
        weight = 0.0;
        weight = tree->weightedcount();
#endif
        sprintf(filename, "numberOfRules_%s_%d.txt", ext, threshold);
        FILE * fono = fopenOrExit(tempDir(filename, options), "wb", "writeAndTest");
#if BRACMATOUTPUT
        sprintf(filename, "rules_%d%s.bra", threshold, ext);
        FILE * fobra = fopenOrExit(tempDir(filename, options), "wb", "Bracmat output");
#endif
        if (256 <= sprintf(filename, "rules_%d%s.lel", threshold, ext))
            {
            printf("writeAndTest: filename 2 small");
            exit(-1);
            }
        FILE * folel = fopenOrExit(tempDir(filename, options), "wb+", "writeAndTest");
#if RULESASTEXT
        filename[strlen(filename)-1] += 2; // change ".lel" to ".len"
        FILE * foleltxt = fopenOrExit(filename,"wb","Text version");
        filename[strlen(filename)-1] -= 2; // change ".len" back to ".lel"
#endif
        if (fono
#if RULESASTEXTINDENTED
            && foo
#endif
#if BRACMATOUTPUT
            && fobra
#endif
            && folel
#if RULESASTEXT
            && foleltxt
#endif
            )
            {
#if RULESASTEXTINDENTED
            fprintf(foo,"tree={%d %f}\n",Nnodes,weight); // "rules_%d%s.txt"
#endif
            fprintf(fono, "%d\n%f", Nnodes, weight); // "numberOfRules_%d.txt"
            fclose(fono);
            int nr = 0;
            strng L("");
            strng R("");
            printRules
                (tree
#if RULESASTEXTINDENTED
                , foo // "rules_%d%s.txt"
#endif
#if BRACMATOUTPUT
                , fobra // "rules_%d%s.bra"
#endif
                , folel // "rules_%d%s.lel"
                , 0
                , &L
                , &R
                , nr
                , options
                );
#if BRACMATOUTPUT
            fclose(fobra);
#endif
#if RULESASTEXTINDENTED
            fclose(foo);
#endif
            if (nflexrules)
                {
                rearrange
                    (nflexrules// the binary output of the training, 
                    // third command line argument
                    , folel     // .lel file, textual output with relative 
                    // positions
#if RULESASTEXT
                    , foleltxt  // .len file, textual output with absolute 
                    // positions, like in the binary output.
#endif
                    );
                }
            else
                {
                filename[strlen(filename) - 1]++; // change ".lel" to ".lem"
                rearrange
                    (filename  // .lem file, the binary output of the training
                    , folel     // .lel file, textual output with relative 
                    // positions
#if RULESASTEXT
                    , foleltxt  // .len file, textual output with absolute 
                    // positions, like in the binary output.
#endif
                    );
                filename[strlen(filename) - 1]--; // change ".lem" back to ".lel"
                }
            fclose(folel);
            if (remove(tempDir(filename, options))) // del ".lel"
                {
                if (options->verbose())
                    {
                    printf("cannot remove %s\n", filename);
                    getchar();
                    }
                }
#if RULESASTEXT
            fclose(foleltxt);
#endif
            return true;
            }
        else
            {
#if RULESASTEXTINDENTED
            if(foo)
                fclose(foo);

#endif
#if BRACMATOUTPUT
            if (fobra)
                fclose(fobra);
#endif
            if (folel)
                fclose(folel);
#if RULESASTEXT
            if(foleltxt)
                fclose(foleltxt);
#endif
            }
        }
    return false;
    }

#if DOTEST
static void testf(node * tree,trainingPair * TestPair,const char * ext,int threshold,char * nflexrules)
    {
    char filename[256];
    int wrong = 0;
    int right = 0;
    int both = 0;
    if(TestPair) // requires PERC > 0
        {
        lemmatise(tree,TestPair);
        sprintf(filename,"test_%d%s.txt",threshold,ext);
        FILE * ftest = fopen(filename,"ab+");
        if(!ftest)
            fprintf(stderr,"Error (testf): Cannot open \"%s\" for appending\n",filename);
        showResults(TestPair,wrong,right,both,ftest);
        int tot = right+wrong+both;
        if(ftest)
            fprintf(ftest,"test pairs %d threshold %d vertices %d right %d (%f) wrong %d (%f) both %d (%f)\n\n",tot,threshold,tree->count(),right,(right*100.0)/tot,wrong,(wrong*100.0)/tot,both,(both*100.0)/tot);
        fresults = fopen("results.txt","ab+");
        if(!fresults)
            fprintf(stderr,"Error (testf): Cannot open \"%s\" for appending\n","results.txt");
        if(fresults)
            {
            fprintf(fresults,"tot %d\tthreshold %d\ttree->count() %d\tright %d\t%% %f\twrong %d\t%% %f\n",tot,threshold,tree->count(),right,(right*100.0)/tot,wrong,(wrong*100.0)/tot);
            }
        if(maxCorrectness < (right*100.0)/tot)
            {
            maxCorrectness = (right*100.0)/tot;
            if(fresults)
                {
                fprintf(fresults,"\nhighest correctness %f\n",maxCorrectness);
                }
            }
        if(fresults)
            fclose(fresults);
        if(ftest)
            fclose(ftest);
        }
    else // requires PERC <= 0
        {
        if(options->verbose())
            {
            sprintf(filename,"rules_%d%s.lem",threshold,ext);
            printf("readRules %d\n",threshold);
            if(readRules(filename) || readRules(nflexrules))
                {
                printf("readRules done %d\n",threshold);
                char word[100];
                printf("\nType first word:\n");
                while(gets(word)[0])
                    {
                    const char * result = applyRules(word);
                    printf("%s\n",result);
                    printf("\nType word:\n");
                    }
                deleteRules();
                }
            }
        }
    }
#endif

static bool doTraining
(const char * fname
, const char * ext
, int cutoff
, const char * nflexrulesFormat
, const char * columns
, char * pairsToTrainInNextPassName
, char * ingestedFractionOfAmbiguousPairsName
, char * allPairsName
, char * allIngestedPairsName
, char * wordsGroupedByRuleName
, char * numbersName
, int & Nnodes
, double & weight
, const char * tag
, int * filelines
, optionStruct * options
)
    {
    bool moreToDo = false;
    Nnodes = 0;

    VertexPointerCount = 0;

    aFile afile(fname, options);

    if (filelines)
        *filelines = afile.lines;

    int pairs;
    trainingPair * TrainingPair = readTrainingPairs(afile, pairs, columns, tag, options);
    markTheAmbiguousPairs(TrainingPair, ext, pairs, options);
    writeAllAvailablePairs(TrainingPair, ext, pairs, options);
    if (options->percentageTestPairs() > 0)
        writeAllTestPairs(TrainingPair, ext, pairs, options);
    hash Hash(10);
    trainingPair * train = NULL;
    trainingPair * test = NULL;
    // Split list of pairs in those that are to be used for training and those
    // that are to be used for testing.
    // Pairs that are doublets are not added to either list.
    // Nor are pairs that are not well-formed (e.g. contain a ' ').
    trainingPair::makeChains(pairs, TrainingPair, &train, &test, options);
    node * top;
    doTheRules(&Hash, train, &top, options);

    FILE * nexttrain = pairsToTrainInNextPassName ? fopenOrExit(tempDir(pairsToTrainInNextPassName, options), "wb", "nexttrain") : NULL;
    FILE * done = ingestedFractionOfAmbiguousPairsName ? fopenOrExit(tempDir(ingestedFractionOfAmbiguousPairsName, options), "wb", "done") : NULL;
    FILE * combined = allPairsName ? fopenOrExit(tempDir(allPairsName, options), "wb", "combined") : NULL;
    FILE * disamb = allIngestedPairsName ? fopenOrExit(tempDir(allIngestedPairsName, options), "wb", "disamb") : NULL;
    if (nexttrain && done && combined && disamb)
        {
        int donepairs = trainingPair::makeNextTrainingSet(pairs, TrainingPair, nexttrain, done, combined, disamb, options);
        if (donepairs)
            {
            moreToDo = true;
            }
        }

    if (nexttrain)
        fclose(nexttrain); /* The training that still has to be done, containing all unambiguous pairs and all that remains of the ambiguous pairs. */
    if (done)
        fclose(done); /* Those parts of ambiguous pairs that are done in this pass. These are typically the most regular parts of ambiguous pairs. */
    if (combined)
        fclose(combined); /* The sum of the above two. */
    if (disamb)
        fclose(disamb);  /* The training that has been done, containing all unambiguous pairs and all parts of ambiguous pairs that were done in this pass. */

    FILE * wordsFile = fopenOrExit(tempDir(wordsGroupedByRuleName, options), "wb", "words file");
    ambivalentWords = 0;
    alternatives = 0;
    allwords = 0;
    top->printSep(wordsFile, 0);
    fprintf(wordsFile
            , "\nAll words: %d    Words with alternative lemmas: %d    Words that are not lemmatized correctly: %d\n"
            , allwords
            , ambivalentWords
            , alternatives);
    fclose(wordsFile); /* Lists all words, grouped by the rule that creates each word's lemma. */

    FILE * fcounting = fopenOrExit(tempDir(numbersName, options), "wb", "counting");
    fprintf(fcounting, "Bottom up left to right traversal of tree.\n");
    fprintf(fcounting, "Nodes\tPairs\tlog(Nodes)\tlog(Pairs)\n");
    int locnodes = 0, locpairs = 0;
    top->Counting(locnodes, locpairs, fcounting);
    fclose(fcounting);

    if (nflexrulesFormat)
        {
        if (cutoff >= 0)
            {
            char naam[500];
            if (sizeof(naam) <= (size_t)sprintf(naam, nflexrulesFormat, 0))
                {
                printf("doTraining: naam 2 small");
                exit(-1);
                }
            writeAndTest(top, ext, 0, naam, Nnodes, weight, options);
#if DOTEST
            testf(top,test,ext,0,naam);
#endif
            for (int thresh = 1; thresh <= cutoff; thresh++)
                {
                top->pruneAll(thresh);
                top = top->cleanup(NULL, options);
                if (sizeof(naam) <= (size_t)sprintf(naam, nflexrulesFormat, thresh))
                    {
                    printf("doTraining: naam 2 small");
                    exit(-1);
                    }
                writeAndTest(top, ext, thresh, naam, Nnodes, weight, options);
#if DOTEST
                testf(top,test,ext,0,naam);
#endif
                }
            }
        else
            {
            char naam[500];
            if (sizeof(naam) <= (size_t)sprintf(naam, nflexrulesFormat, cutoff))
                {
                printf("doTraining: naam 2 small");
                exit(-1);
                }
            writeAndTest(top, ext, 0, naam, Nnodes, weight, options);
#if DOTEST
            testf(top,test,ext,0,naam);
#endif
            int max = 3;
            if (cutoff >= 0)
                max = cutoff;
            for (int thresh = 1; thresh <= max; thresh++)
                {
                top->pruneAll(thresh);
                writeAndTest(top, ext, thresh, naam, Nnodes, weight, options);
#if DOTEST
                testf(top,test,ext,thresh,naam);
#endif
                }
            }
        }
    else
        {
        writeAndTest(top, ext, 0, 0, Nnodes, weight, options);
#if DOTEST
        testf(top,test,ext,0,0);
#endif
        int max = 3;
        if (cutoff >= 0)
            max = cutoff;
        for (int thresh = 1; thresh <= max; thresh++)
            {
            top->pruneAll(thresh);
            writeAndTest(top, ext, thresh, 0, Nnodes, weight, options);
#if DOTEST
            testf(top,test,ext,thresh,0);
#endif
            }
        }
    delete top;
    building = false; // Signal to ~vertexPointer() to not access nodes.
    delete[] TrainingPair;
    building = true; // Signal to ~vertexPointer() to access nodes.
    return moreToDo;
    }

const int partOfFile(const char * fbuf, const double fraction, optionStruct * options)
    {
    if (options->currentParms() && !flog)
        {
        flog = fopenOrExit(options->currentParms(), "a", "log file");
        fprintf(flog, "%s: blobs=%d lines=%d fraction=%f most penalized=%d\n", options->wordList(), options->blobs(), options->lines(), fraction, node::mostPenalized);
        fclose(flog);
        flog = 0;
        }
    FILE * f2 = fopenOrExit(fbuf, "w", "computeParms");
    double bucket = fraction;
    int kar;
    FILE * f = fopen(options->wordList(), "r");
    if ((double)options->blobs() * fraction > 1.0)
        {
        int bl = 1;
        int li = 0;
        int blbs = 0; // If there are no non-empty lines, there are no blobs
                     // either.
        int prevkar = 0;
        unsigned int lineProfile = 0;
        while ((kar = fgetc(f)) != EOF)
            {
            if (kar == '\n')
                {
                if(lineProfile > 1)
                    { // Just finished non-empty line
                    lineProfile >>= 1; // push line before previous line out of
                                       // memory
                    // if lineProfile is 1, we may have a bunch of empty lines
                    // before a new non-empty line appears
                    }
                if (prevkar == '\n')
                    {
                    if (bucket >= 1.0)
                        {
                        ++bl;
                        bucket -= 1.0;
                        }
                    bucket += fraction;
                    }
                else if (prevkar != 0)
                    {
                    if (bucket >= 1.0)
                        {
                        ++li;
                        fputc(kar, f2);
                        }
                    }
                }
            else
                {
                if(kar != '\r')
                    {
                    if(lineProfile == 1)
                        {
                        // previous line was empty. Before that, there was a 
                        // non-empty line. Blob boundary detected!
                        ++blbs;
                        }
                    lineProfile |= 4; // set third bit when making non-empty line
                    }
                if (bucket >= 1.0)
                    fputc(kar, f2);
                }
            prevkar = kar;
            }
        if(lineProfile & 4)
            { // last line was not finished off with newline
            ++li;
            }
        if(li > 0)
            ++blbs;
        printf("li %d blbs %d\n",li,blbs);
        if (options->currentParms() && !flog)
            {
            flog = fopenOrExit(options->currentParms(), "a", "log file");
            fprintf(flog, "Read %d blobs\n", bl);
            fclose(flog);
            flog = 0;
            }
        options->setReadBlobs(bl);
        }
    else if ((double)options->lines() * fraction > 1.0)
        {
        while ((kar = fgetc(f)) != EOF)
            {
            if (bucket >= 1.0)
                fputc(kar, f2);
            if (kar == '\n')
                {
                if (bucket >= 1.0)
                    bucket -= 1.0;
                bucket += fraction;
                }
            }
        }
    else
        {
        while ((kar = fgetc(f)) != EOF)
            {
            fputc(kar, f2);
            }
        }
    fclose(f2);

    f2 = fopenOrExit(fbuf, "r", "computeParms");
    int fraclines = 0;
    while ((kar = fgetc(f2)) != EOF)
        {
        if (kar == '\n')
            ++fraclines;
        }
    fclose(f2);
    options->setReadLines(fraclines);
    if (options->currentParms() && !flog)
        {
        flog = fopenOrExit(options->currentParms(), "a", "log file");
        fprintf(flog, "Use %d lines of %d\n", fraclines, options->lines());
        fclose(flog);
        flog = 0;
        }
    if (options->verbose())
        {
        printf("reading %d lines\n", fraclines);
        printf("tmpnam %s\n", fbuf);
        }
    return fraclines;
    }

void computeParms(optionStruct * options)
    {
    CHECK("iglobTempDir");
    int maxswath = options->swaths();//MAXSWATH;
    int currentNo = 0;
    int brownNo = 0;
    double currentweight = 0.0;
    double brownweight = 0.0;
    double fraction = 0.0; // 0.0 <= fraction <= 1.0
    double factor = 0.0;
    double iterationsfactor = 1;
    double miniterations = options->minIterations();//MINITERATIONS;
    double maxiterations = options->maxIterations();//MAXITERATIONS;
    const char * tag = "";
    node::mostPenalized = options->expectedCutoff() + 1; // parameter to weight function
    if (options->minfraction() > 0.0)
        {
        factor = pow(options->maxfraction() / options->minfraction(), 1.0 / (double)maxswath);
        // minfraction * factor^0 == minfraction
        // minfraction * factor^maxswath == maxfraction
        }

    if (miniterations > 0)
        {
        iterationsfactor = pow(maxiterations / miniterations, 1.0 / (double)maxswath);
        // maxiterations * iterationsfactor^0 == maxiterations
        // maxiterations * iterationsfactor^-maxswath == miniterations
        }
    else
        iterationsfactor = 1;
    const char * filename = options->wordList();
    const char * fbuf;
    fbuf = dup(tempDir("trainFraction", options));

    char ext[100];
    ext[0] = '\0';
    if (sizeof(ext) <= (size_t)sprintf(ext, "%s", options->extra() ? options->extra() : ""))
        {
        printf("computeParms: ext 2 small");
        exit(-1);
        }
    int Nnodes = 0;
    double weight = 0.0;
    int br1 = 0, br2 = 0;
    for (int swath = 0; swath <= maxswath; ++swath)
        {
        int blobs = 1;
        int lines = 0;
        int fraclines = 0;
        if (options->minfraction() > 0.0)
            {
            fraction = (swath == maxswath) ? options->maxfraction() : options->minfraction() * pow(factor, (double)swath);
            if (fraction == 1.0)
                filename = options->wordList();
            else
                {
                filename = fbuf;
                fraclines = partOfFile(fbuf, fraction, options);
                }
            CHECK("D1globTempDir");
            brown(); // until not all parms are zero
            ++br1;
            if (swath == 0)
                init(options);
            else
                copybest(); // go on with best result so far.
            char wordsGroupedByRuleName[1024];
            if (sizeof(wordsGroupedByRuleName) <= (size_t)sprintf(wordsGroupedByRuleName, "words_%s%s.txt", ext, tag))
                {
                printf("computeParms: wordsGroupedByRuleName 2 small");
                exit(-1);
                }
            char numbersName[1024];
            if (sizeof(numbersName) <= (size_t)sprintf(numbersName, "numbers_%s%s.tab", ext, tag))
                {
                printf("computeParms: numbersName 2 small");
                exit(-1);
                }
            int filelines;
            doTraining
                (/* const char *                                */  filename
                ,/* const char *                                */  ext
                ,/* int cutoff                                  */  0
                ,/* const char * nflexrulesFormat               */  options->flexrules()
                ,/* const char *                                */  options->columns()
                ,/* char * pairsToTrainInNextPassName           */  NULL
                ,/* char * ingestedFractionOfAmbiguousPairsName */  NULL
                ,/* char * allPairsName                         */  NULL
                ,/* char * allIngestedPairsName                 */  NULL
                ,/* char *                                      */  wordsGroupedByRuleName
                ,/* char *                                      */  numbersName
                ,/* int &                                       */  Nnodes
                ,/* double &                                    */  weight
                ,/* const char *                                */  tag
                ,/* int * filelines                             */  &filelines
                , options
                ); // sets Nnodes
            if (lines == 0)
                fraclines = lines = filelines;

            brownNo = Nnodes;
            currentNo = brownNo;
            brownweight = weight;
            currentweight = brownweight;
            betterfound(currentNo, currentweight, swath, -1, blobs, lines, fraction, fraclines, false, options);
            printparms(Nnodes, weight, options);
            }
        int looplimit = (int)(maxiterations*pow(iterationsfactor, -swath));
#if FLOATINGPOINTPARMS
        for (int iterations = 0; iterations < looplimit; ++iterations)
#else
        for(int iterations = 0;iterations < 64;++iterations)
#endif
            {
            CHECK("D2aglobTempDir");
            brown();
            ++br2;
            CHECK("D2bglobTempDir");
            if (options->currentParms() && !flog)
                {
                flog = fopenOrExit(options->currentParms(), "a", "log file");
                CHECK("D2dglobTempDir");
                fprintf(flog, "//iteration:%d.%d %s\n", swath, iterations, options->doweights() ? "weights" : "count");
                CHECK("D2eglobTempDir");
                fclose(flog);
                flog = 0;
                }
            CHECK("D2fglobTempDir");
            if (options->verbose())
                {
                printf("%d.%d ", swath, iterations);
                }
            CHECK("D2globTempDir");

            char wordsGroupedByRuleName[1024];
            if (sizeof(wordsGroupedByRuleName) <= (size_t)sprintf(wordsGroupedByRuleName, "words_%s%s.txt", ext, tag))
                {
                printf("computeParms: wordsGroupedByRuleName 2 small");
                exit(-1);
                }
            char numbersName[1024];
            if (sizeof(numbersName) <= (size_t)sprintf(numbersName, "numbers_%s%s.tab", ext, tag))
                {
                printf("computeParms: numbersName 2 small");
                exit(-1);
                }

            int filelines;
            doTraining
                (/* const char *                                */  filename
                ,/* const char *                                */  ext
                ,/* int cutoff                                  */  0
                ,/* const char * nflexrulesFormat               */  options->flexrules()
                ,/* const char *                                */  options->columns()
                ,/* char * pairsToTrainInNextPassName           */  NULL
                ,/* char * ingestedFractionOfAmbiguousPairsName */  NULL
                ,/* char * allPairsName                         */  NULL
                ,/* char * allIngestedPairsName                 */  NULL
                ,/* char *                                      */  wordsGroupedByRuleName
                ,/* char *                                      */  numbersName
                ,/* int &                                       */  Nnodes
                ,/* double &                                    */  weight
                ,/* const char *                                */  tag
                ,/* int * filelines                             */  &filelines
                , options
                ); // sets Nnodes
            if (lines == 0)
                fraclines = lines = filelines;

            printparms(Nnodes, weight, options);
            if (options->verbose())
                {
                printf("\r%d %d %f %d %f           \n", iterations, currentNo, currentweight, brownNo, brownweight);
                }

            if (currentNo == 0)
                {
                currentNo = Nnodes;
                currentweight = weight;
                }
            else
                {
                brownNo = Nnodes;
                brownweight = weight;
                if (options->verbose())
                    printf("swath %d brownNo %d currentNo %d\n", swath, brownNo, currentNo);
                if ((!options->doweights() && brownNo <= currentNo) || (options->doweights() && brownweight <= currentweight))
                    {
                    bool improvement = (!options->doweights() && (brownNo < currentNo)) || (options->doweights() && (brownweight < currentweight));
                    if (options->verbose())
                        printf("%s\n", improvement ? "IMPROVEMENT" : "same");
                    currentNo = brownNo;
                    currentweight = brownweight;
                    betterfound(currentNo, currentweight, swath, iterations, blobs, lines, fraction, fraclines, improvement, options);
                    }
                else
                    worsefound();
                }
            }
        if (options->verbose())
            printf("br1 %d br2 %d\n", br1, br2);
        }
    remove(fbuf);
    delete[] fbuf;
    }

static void trainRules(const char * tag, optionStruct * options)
    {
    CHECK("jglobTempDir");
    assert(options->flexrules() != NULL); // 20130125
    const char * nflexrules = options->flexrules();
    const char * fname = options->wordList();
    bool moreToDo = true;
    int passes = 0;
    char pairsToTrainInNextPassName[1024];
    char ingestedFractionOfAmbiguousPairsName[1024];
    char allPairsName[1024];
    char allIngestedPairsName[1024];
    char wordsGroupedByRuleName[1024];
    char numbersName[1024];


    char pairsToTrainInNextPassFormat[1024];
    char ingestedFractionOfAmbiguousPairsFormat[1024];
    char allPairsFormat[1024];
    char allIngestedPairsFormat[1024];
    char wordsGroupedByRuleFormat[1024];
    char bestRulesFormat[1024];



    char numbersFormat[1024];
    char command[1024];
    const char * FlexrulePassFormat = "%s.pass%d.cutoff%%d"; // flexrules.passN.cutoffM
    const char * AccumulatedFlexrulePassFormat = "%s.pass%d.cutoff%%d.accumulated"; // flexrules.passN.cutoffM
    clock_t start = clock();
    int Nnodes = 0;
    double weight = 0.0;
    sprintf(pairsToTrainInNextPassFormat, "pairsToTrainInNextPass.%s%s.pass%%d", options->extra(), tag);
    sprintf(ingestedFractionOfAmbiguousPairsFormat, "ingestedFractionOfAmbiguousPairs.%s%s.pass%%d", options->extra(), tag);
    sprintf(allPairsFormat, "allPairs.%s%s.pass%%d", options->extra(), tag);
    sprintf(allIngestedPairsFormat, "allIngestedPairs.%s%s.pass%%d", options->extra(), tag);
    sprintf(wordsGroupedByRuleFormat, "wordsGroupedByRuleName.%s%s.pass%%d", options->extra(), tag);
    sprintf(numbersFormat, "numbers.%s%s.pass%%d", options->extra(), tag);



    char nflexrulesTag[1256];
    if (nflexrules)
        {
        if (tag && *tag)
            sprintf(nflexrulesTag, "%s.%s", nflexrules, tag);
        else
            strcpy(nflexrulesTag, nflexrules);
        nflexrules = nflexrulesTag;
        }

    const char * accumulatedFormat = FlexrulePassFormat;
    const char * accumulatedFormatPrev = NULL;

    do
        {
        ++passes;
        char flexrulesPass[1256];
        sprintf(flexrulesPass, FlexrulePassFormat, nflexrules, passes);
        char ext[100];
        ext[0] = '\0';
        if (options->extra())
            strcpy(ext, options->extra());

        sprintf(pairsToTrainInNextPassName, pairsToTrainInNextPassFormat, passes);
        sprintf(ingestedFractionOfAmbiguousPairsName, ingestedFractionOfAmbiguousPairsFormat, passes);
        sprintf(allPairsName, allPairsFormat, passes);
        sprintf(allIngestedPairsName, allIngestedPairsFormat, passes);
        sprintf(wordsGroupedByRuleName, wordsGroupedByRuleFormat, passes);
        sprintf(numbersName, numbersFormat, passes);
        char spass[10];
        if (sizeof(spass) <= (size_t)sprintf(spass, ".pass%d", passes))
            {
            printf("trainRules: spass 2 small");
            exit(-1);
            }
        if (sizeof(ext) <= (size_t)sprintf(ext, "%s%s", options->extra() ? options->extra() : "", spass))
            {
            printf("trainRules: ext 2 small");
            exit(-1);
            }


        moreToDo = doTraining
            (/* const char *                                */  fname
            ,/* const char *                                */  ext
            ,/* int cutoff                                  */  options->cutoff()
            ,/* const char * nflexrulesFormat               */  flexrulesPass
            ,/* const char *                                */  passes > 1 ? "12" : options->columns()
            ,/* char * pairsToTrainInNextPassName           */  pairsToTrainInNextPassName
            ,/* char * ingestedFractionOfAmbiguousPairsName */  ingestedFractionOfAmbiguousPairsName
            ,/* char * allPairsName                         */  allPairsName
            ,/* char * allIngestedPairsName                 */  allIngestedPairsName
            ,/* char *                                      */  wordsGroupedByRuleName
            ,/* char *                                      */  numbersName
            ,/* int &                                       */  Nnodes
            ,/* double &                                    */  weight
            ,/* const char *                                */  passes > 1 ? NULL : tag
            ,/* int * filelines                             */  NULL
            , options
            );
        /*
        Re-do the training, but only with those pairs that made it into
        the set of ingested pairs. The idea is to avoid the "noise" caused
        by the siblings of ambiguous pairs that didn't make it.
        The result of this excercise can be a better rule tree.
        */
        if (options->redo() && moreToDo)
            {
            if (options->verbose())
                {
                printf("More training to do with file \"%s\"\n", pairsToTrainInNextPassName);
                }
            sprintf(wordsGroupedByRuleName, "words_%s%s.txt", ext, tag);
            sprintf(numbersName, "numbers_%s%s.tab", ext, tag);
            if (doTraining
                (/* const char *                                */  tempDir(allIngestedPairsName, options)
                ,/* const char *                                */  ext
                ,/* int                                         */  options->cutoff()
                ,/* const char * nflexrulesFormat               */  flexrulesPass
                ,/* const char *                                */  options->columns()
                ,/* char * pairsToTrainInNextPassName           */  NULL
                ,/* char * ingestedFractionOfAmbiguousPairsName */  NULL
                ,/* char * allPairsName                         */  NULL
                ,/* char * allIngestedPairsName                 */  NULL
                ,/* char *                                      */  wordsGroupedByRuleName
                ,/* char *                                      */  numbersName
                ,/* int &                                       */  Nnodes
                ,/* double &                                    */  weight
                ,/* const char *                                */  tag
                ,/* int * filelines                             */  NULL
                , options
                )
                ) // sets Nnodes
                {
                if (options->verbose())
                    {
                    printf("This should return 0\n");
                    getchar();
                    }
                }
            }
        else if (options->verbose())
            {
            if (moreToDo)
                printf("No retraining done on ingested pairs, although ambiguous pairs were found and may have caused noise. (Faster) \n");
            }
        char filename[256];
        if (sizeof(filename) <= (size_t)sprintf(filename, "statistics_%s.txt", ext))
            {
            printf("trainRules: filename 2 small");
            exit(-1);
            }
        FILE * fo = fopenOrExit(tempDir(filename, options), "ab", "trainRules");
        fprintf(fo, "VertexPointerCount          %d\n", VertexPointerCount);
        fprintf(fo, "VertexCount                 %d\n", VertexCount);
        fprintf(fo, "TrainingPairCount           %d\n", TrainingPairCount);
        fprintf(fo, "HashCount                   %d\n", HashCount);
        fprintf(fo, "RulePairCount               %d\n", RulePairCount);
        fprintf(fo, "StrngCount                  %d\n", StrngCount);
        fprintf(fo, "RuleTemplateCount           %d\n", RuleTemplateCount);
        fprintf(fo, "ShortRulePairCount          %d\n", ShortRulePairCount);
        fprintf(fo, "FullRulePairCount           %d\n", FullRulePairCount);

        clock_t finish = clock();
        double  duration = 0;
        duration = (double)(finish - start) / CLOCKS_PER_SEC;

        fprintf(fo, "%2.1f seconds\n", duration);
        fclose(fo);
        if (options->verbose())
            {
            printf("%2.1f seconds\n", duration);
            }

        char AccumulatedFlexrulesPass[256];
        char NextAccumulatedFlexrulesPassFormat[256];
        if (sizeof(bestRulesFormat) <= (size_t)sprintf(bestRulesFormat, accumulatedFormat, nflexrules, passes))
            {
            printf("trainRules: bestRulesFormat 2 small");
            exit(-1);
            }
        if (passes > 1)
            {
            assert(accumulatedFormatPrev != NULL);
            if (sizeof(AccumulatedFlexrulesPass) <= (size_t)sprintf(AccumulatedFlexrulesPass, accumulatedFormatPrev, nflexrules, passes - 1))
                {
                printf("trainRules: AccumulatedFlexrulesPass 2 small");
                exit(-1);
                }
            if (sizeof(NextAccumulatedFlexrulesPassFormat) <= (size_t)sprintf(NextAccumulatedFlexrulesPassFormat, AccumulatedFlexrulePassFormat, nflexrules, passes))
                {
                printf("trainRules: NextAccumulatedFlexrulesPassFormat 2 small");
                exit(-1);
                }
            }

        for (int cut = 0; cut <= options->cutoff(); ++cut)
            {
            char nextbestflexrules[1150];
            if (sizeof(nextbestflexrules) <= (size_t)sprintf(nextbestflexrules, flexrulesPass, cut))
                {
                printf("trainRules: nextbestflexrules 2 small");
                exit(-1);
                }
            prettyPrint(nextbestflexrules);
            if (passes > 1)
                {
                char bestflexrules[1150], newbestflexrules[1150];
                if (sizeof(bestflexrules) <= (size_t)sprintf(bestflexrules, AccumulatedFlexrulesPass, cut))
                    {
                    printf("trainRules: bestflexrules 2 small");
                    exit(-1);
                    }
                if (sizeof(newbestflexrules) <= (size_t)sprintf(newbestflexrules, NextAccumulatedFlexrulesPassFormat, cut))
                    {
                    printf("trainRules: newbestflexrules 2 small");
                    exit(-1);
                    }
                if (options->verbose())
                    {
                    printf("flexcombi best %s + next best %s -> combined %s\n", bestflexrules, nextbestflexrules, bestflexrules);
                    }
                if (!flexcombi(bestflexrules, nextbestflexrules, newbestflexrules))
                    break;
                prettyPrint(newbestflexrules);
                }
            }
        fname = tempDir(pairsToTrainInNextPassName, options);

        accumulatedFormatPrev = accumulatedFormat;
        accumulatedFormat = AccumulatedFlexrulePassFormat;
        } while (moreToDo && passes < 30);



        if (options->verbose())
            {
            putchar('\n'); // to compensate for missing newline in readTrainingPairs
            }
        if (options->cutoff() >= 0)
            {
            char dest[1000];
            dest[0] = '\0';
            for (int cut = 0; cut <= options->cutoff(); ++cut)
                {
                sprintf(dest, bestRulesFormat, cut);

                char dirname[500];
                const char * lastslash = strrchr(nflexrules, DIRSEP);
                const char * filename;
                if (lastslash)
                    {
                    filename = lastslash + 1;
                    sprintf(dirname, "%.*s%c%d", (int)(lastslash - nflexrules), nflexrules, DIRSEP, cut);
                    }
                else
                    {
                    filename = nflexrules;
                    sprintf(dirname, "%d", cut);
                    }
                if (options->verbose())
                    printf("dirname %s\n", dirname);

                char testfile[500];
                sprintf(testfile, "%s%cTESTFILE", dirname, DIRSEP);
                if (options->verbose())
                    printf("testfile %s\n", testfile);
                FILE * fptest = fopen(testfile, "w");
                bool hasDir = false;
                if (fptest)
                    {
                    if (options->verbose())
                        printf("testfile created\n");
                    fclose(fptest);
                    remove(testfile);
                    if (options->verbose())
                        printf("testfile deleted\n");
                    hasDir = true;
                    }
                else
                    {
                    if (options->verbose())
                        printf("testfile not created\n");
                    sprintf(command, "%s %s", MKDIR, dirname);
                    if (options->verbose())
                        printf("command %s\n", command);
                    system(command);
                    fptest = fopen(testfile, "w");
                    if (fptest)
                        {
                        if (options->verbose())
                            printf("testfile created\n");
                        fclose(fptest);
                        remove(testfile);
                        if (options->verbose())
                            printf("testfile deleted\n");
                        hasDir = true;
                        }
                    }
                if (hasDir)
                    {
                    sprintf(command, "%s%c%s", dirname, DIRSEP, filename);
                    remove(command);
                    rename(dest, command);
                    }
                }
            }
        for (int cut = 0; cut <= options->cutoff(); ++cut)
            {
            char scut[20];
            sprintf(scut, "%s%d", SCUT, cut);
            for (int passes = 1;; ++passes)
                {
                char dest[1000];
                char spass[10];
                sprintf(spass, ".pass%d", passes);
                sprintf(dest, "%s%s%s", nflexrules, scut, spass);
                if (options->verbose())
                    printf("Remove %s \n", dest);
                if (remove(dest))
                    break;
                }
            }
    }


FILE * flog = NULL;

static void initOutput(const char * path)
    {
    if (!path)
        return;
    FILE * fp = fopenOrExit(path, "w", "initOutput");
    fclose(fp);
    }

int main(int argc, char **argv)
    {
    if (argc < 2)
        {
        printf("affixtrain - supervised learning of affix rules for AFFIXTRAIN, version " VERSION "\n");
        printf("%ssing floating point parameters.\n", FLOATINGPOINTPARMS ? "U" : "Not u");
        printf("%s -h for usage\n", argv[0]);
        exit(0);
        }
    optionStruct options;
    switch (options.readArgs(argc, argv))
        {
        case GoOn:
        break;
        case Leave:
        exit(0);
        case Error:
        fprintf(stderr, "Error: Error in options. Exiting\n");
        exit(1);
        }

    if (options.verbose())
        {
        options.print(stdout);
        options.printArgFile();
        }


    if (options.rawRules())
        {
        prettyPrint(options.rawRules());
        prettyPrintBracmat(options.rawRules());
        }
    else
        {
        FILE * fptest = fopen(tempDir("testFile", &options), "wb");
        if (fptest)
            {
            fclose(fptest);
            remove(tempDir("testFile", &options));
            }
        else
            {
            if (options.tempDir())
                printf("Cannot create file %s. Did you specify an existing and writable temp directory? (option -j)\n", tempDir("testFile", &options));
            else
                printf("Cannot create file %s. Is the working directory writable?\n", "testFile");
            exit(-1);
            }

        if (options.computeParms())
            {
            if (options.currentParms())
                initOutput(options.currentParms());

            if (options.bestParms())
                initOutput(options.bestParms());
            computeParms(&options);
            }
        else
            {
            setCompetitionFunction(&options);
            tagClass * Tags = NULL;
            options.setReadLines(options.lines());
            Tags = collectTags(&options);

            tagClass * theTag = Tags;
            if (theTag)
                {
                if (options.verbose())
                    printf("Doing Tags\n");
                while (theTag)
                    {
                    if (options.verbose())
                        {
                        printf("Doing tag %s\n", theTag->name);
                        }
                    trainRules(theTag->name, &options);
                    theTag = theTag->next;
                    }
                }
            else
                {
                if (options.verbose())
                    printf("NOT doing Tags\n");
                trainRules("", &options);
                }
            }

        if (options.verbose())
            {
            if (argc < 3)
                getchar();
            printf("\nAffixTrain OK\n");
            }
        }
    }
