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

class optionStruct
    {
    private:
        /* 
        -c <cutoff> 
        -e <extra> 
        -f <compfunc>
        -i <word list> 
        -n <columns> 
        -o <flexrules> 
        */
        int c; // cutoff
        int C; // expected cutoff
        const char * e; // extra
        const char * f; // compfunc
        const char * i; // word list
        const char * j; // temp dir
        const char * n; // columns
        const char * o; // flexrules
        const char * B; // Best parms
        const char * P; // Current parms
        const char * b; // raw rules
        bool ComputeParms;// compute parms
        bool SuffixOnly;// suffix only
        bool SuffixOnlyParmSeen;
        bool Verbose;         // verbose
        double Minfraction;
        double Maxfraction;
        bool Doweights;
        bool Redo; // option R. Set to true if training has to be done once more 
        // after the removal of the homographs that will be handled in
        // the next iteration.
        int Q; // Max recursion depth when attempting to create candidate rule
        int q; // Percentage of training pairs to set aside for testing
        // set to positive value if you want to set PERC percent
        // of the avaliable data aside for testing.
        OptReturnTp doSwitch(int c,char * locoptarg,char * progname);
        OptReturnTp readOptsFromFile(char * locoptarg,char * progname);
    public:
        const int cutoff() const{return c;}
        const int expectedCutoff() const{return C;}
        const char * extra() const{return e;}
        const char * compfunc() const{return f;}
        const char * wordList() const{return i;}
        const char * tempDir() const{return j;}
        const char * columns() const{return n;}
        const char * flexrules() const{return o;}
        const char * bestParms() const{return B;}
        const char * currentParms() const{return P;}
        const char * rawRules() const{return b;}
        const int maxRecursionDepthRuleCreation()const{return Q;}
        const int percentageTestPairs()const{return q;}
        const bool computeParms()const{return ComputeParms;}
        const bool suffixOnly()const{return SuffixOnly;}
        const bool verbose()const{return Verbose;}
        const double minfraction()const{return Minfraction;}
        const double maxfraction()const{return Maxfraction;}
        const bool doweights()const{return Doweights;}
        const bool redo()const{return Redo;}
        
        optionStruct();
        ~optionStruct();
        OptReturnTp readArgs(int argc, char * argv[]);
    };
