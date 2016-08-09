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
#ifndef RULETEMPLATE_H
#define RULETEMPLATE_H

#include <stddef.h>
#include <stdio.h>

extern const int unequal;
extern const int equal;
extern const char Start[2];
extern const char End[2];
extern const char AnyEnd[3];
extern const char StartAnyEnd[4];
extern const char StartAny[3];

void suffix(char * mask);

class optionStruct;
class trainingPair;
class ruleTemplate
    {
    protected:
        char npatternArray[1000];
        char nreplacementArray[1000];
    public:
        const char * itsNPatternArray(){ return npatternArray; }
        const char * itsNReplacementArray(){ return nreplacementArray; }
        void copy(ruleTemplate * Template);
        bool makebigger(int countdown, int & anihilatedGuards, optionStruct * options);
        ruleTemplate();
        ~ruleTemplate();
    };

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
        void match(trainingPair * trainingpair);
        bool mergeTemplates(const char * predef, optionStruct * options);
        void print(FILE *fp)
            {
            fprintf(fp, "\npattern:\n%s\nreplacement:\n%s\n", pattern, replacement);
            fprintf(fp, "npattern:\n%s\nnreplacement:\n%s\n", npatternArray, nreplacementArray);
            }
    };


#endif
