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

typedef enum {GoOn = 0,Leave = 1,Error = 2} OptReturnTp;

struct optionStruct
    {
/* 
-c <cutoff> 
-e <extra> 
-f <compfunc>
-i <word list> 
-n <columns> 
-o <flexrules> 
*/
    const char * c; // cutoff
    const char * e; // extra
    const char * f; // compfunc
    const char * i; // word list
    const char * j; // temp dir
    const char * n; // columns
    const char * o; // flexrules
    const char * B; // Best parms
    const char * P; // Current parms
    const char * b; // raw rules
    const char * t; // pretty printed rules
    bool computeParms;// compute parms
    bool suffixOnly;// suffix only
    bool verbose;         // verbose
    int minperc;
    int maxperc;
    bool doweights;
    optionStruct();
    ~optionStruct();
    OptReturnTp doSwitch(int c,char * locoptarg,char * progname);
    OptReturnTp readOptsFromFile(char * locoptarg,char * progname);
    OptReturnTp readArgs(int argc, char * argv[]);
    };
