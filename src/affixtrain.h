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
#ifndef AFFIXTRAIN_H
#define AFFIXTRAIN_H

#include "settingsaffixtrain.h"

class optionStruct;

class countAndWeight
    {
    double Weight;
    LONG CountByDepth;
    LONG CountBySize;
    int Nnodes;
    public:
        void setNnodes(int N){Nnodes = N;};
        void setWeight(double W){Weight = W;};
        void setCountByDepth(LONG C){CountByDepth = C;};
        void setCountBySize(LONG C){CountBySize = C;};
        int getNnodes()const{return Nnodes;}
        double getWeight()const{return Weight;}
        LONG getCountByDepth()const{return CountByDepth;}
        LONG getCountBySize()const{return CountBySize;}
    };


//void trainRules(const char * tag, optionStruct * options);
void trainRules(const char * tag, optionStruct * options,countAndWeight * Counts);

extern int openfiles;

#endif
