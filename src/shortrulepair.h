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
#ifndef SHORTRULEPAIR_H
#define SHORTRULEPAIR_H

#include <string.h>

class trainingPair;
class ruleTemplate;
class vertex;

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
        bool checkRule(trainingPair * trainingpair, vertex * parentPat);
        ~shortRulePair();
    };

#endif