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

#include "optionaff.h"
#include "applyaffrules.h"
#include "affixtrain.h"
#include "argopt.h"
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

static char opts[] = "?@:A:B:b:C:c:D:d:E:e:F:f:G:hH:I:i:j:K:k:L:M:N:n:O:o:P:p:Q:"/*q:*/"R:s:T:t:V:v:X:x:W:" /* GNU: */ "wr";
static char *** Ppoptions = NULL;
static char ** Poptions = NULL;
static int optionSets = 0;

char * dupl(const char * s)
    {
    if(s)
        {
        char * d = new char[strlen(s) + 1];
        strcpy(d, s);
        return d;
        }
    return 0;
    }

optionStruct::optionStruct(optionStruct & O)
    {
    c = O.c;
    C = O.C;
    e = dupl(O.e);
    f = dupl(O.f);
    I = dupl(O.i);
    this->O = dupl(O.O);
    i = dupl(O.i);
    n = dupl(O.n);
    k = dupl(O.k);
//    o = dupl(O.o);
    o = 0;
    B = dupl(O.B);
    P = dupl(O.P);
    X = dupl(O.X);
    j = dupl(O.j);
    b = dupl(O.b);
    E = dupl(O.E);
    G = dupl(O.G);
    nD = O.nD;
    D = new double[nD];
    for(int ii = 0;ii < nD;++ii)
        D[ii] = O.D[ii];
    ComputeParms = O.ComputeParms;
    SuffixOnly = O.SuffixOnly;
    ExpensiveInfix = O.ExpensiveInfix;
    SuffixOnlyParmSeen = O.SuffixOnlyParmSeen;
    Verbose = O.Verbose;
    Remove = O.Remove;
    VX = O.VX;
    Minfraction = O.Minfraction;
    Maxfraction = O.Maxfraction;
    Redo = O.Redo;
    Test = O.Test;
    F = O.F;
    TrainTest = O.TrainTest;
    Q = O.Q;
    //q = O.q;
    K = O.K;
    M = O.M;
    N = O.N;
    Argstring = 0;
    WeightFunction = econstant;
    Blobs = O.Blobs;
    Lines = O.Lines;
    FracBlobs = O.FracBlobs;
    FracLines = O.FracLines;
    }

optionStruct::optionStruct()
    {
    c = 5; // pruning threshold, cutoff 
    C = -1; // expected pruning threshold when determining the parameters using weightedcount() (-XW)
    // Rules with C+1 supporting word/lemma pairs have the highest penalty.
    // Rules with C or fewer supporting word/lemma pairs are probably best cut off.
    // To decide whether that is the case, the rules must be tested with OOV word/lemma pairs.
    e = NULL; // extra
    f = NULL; // compfunc
    I = NULL; // word list (with -b option)
    O = NULL; // lemmas of word list (with -b option)
    i = NULL; // word/lemma list
    n = NULL; // columns
    k = NULL; // specific POS tag
    o = NULL; // flexrules
    B = NULL;
    P = NULL;
    X = NULL;
    j = NULL; // temp dir
    b = NULL; // raw file (see t)
    E = NULL;
    G = NULL;
    D = NULL;
    nD = 0;
    ComputeParms = true;// compute parms
    SuffixOnly = false;// suffix rules only
    ExpensiveInfix = false;
    SuffixOnlyParmSeen = false;
    Argstring = 0;
    Verbose = 0;// verbose
    Remove = true;
    VX = false;
    Minfraction = 0.01; // L
    Maxfraction = 1.0; // H
    Redo = false;
    Test = true;
    F = true;
    TrainTest = true;
    Q = 1;
    //q = 1;
    K = 20;   // Number of differently sized fractions of trainingdata 
    M = 10.0; // # Iterations when training with Maxfraction of input
    N = 100.0;// # Iterations when training with Minfraction of input

    WeightFunction = econstant;
    Blobs = 0; // Number of blobs found in word/lemma list
    Lines = 0; // Number of lines found in word/lemma list
    FracBlobs = 0; // Number of blobs used for training
    FracLines = 0; // Number of lines used for training
    }

void cleanUpOptions()
    {
	if (optionSets > 0)
		{
		for (int I = 0; I < optionSets; ++I)
			{
			delete[] Poptions[I];
			delete[] Ppoptions[I];
			}
		delete[] Poptions;
		delete[] Ppoptions;
		}
    }

optionStruct::~optionStruct()
    {
    delete[] Argstring;
    delete[] e;
    delete[] f;
    delete[] i;
    delete[] j;
    delete[] n;
    delete[] k;
    delete[] o;
    delete[] B;
    delete[] P;
    delete[] b;
    delete[] G;
    delete[] E;
    delete[] X;
    delete[] I;
    delete[] O;
    }

void optionStruct::detectFloatingPointNumbers(const char * S)
    {
    int n = 0;
    const char * t = 0;
    double Sum2 = 0.0;
    double Sum = 0.0;
    double d;
    const char * s;
    const char * endptr = S;
    char * newendptr = 0;
    for (s = S;;s = t + 1,endptr = newendptr)
        {
        d = strtod(s, &newendptr);
        if  (   newendptr && newendptr > s 
            &&  newendptr > endptr 
            &&  (   newendptr[-1] == '.' 
                ||  (   newendptr[-1] >= '0' 
                    &&  newendptr[-1] <= '9'
                    )
                )
            )
            {
            Sum2 += d*d;
            ++n;
            }
        t = strpbrk(s, ",;:");
        if(t == 0)
            break;
        }
    if (Sum2 <= 0)
        {
        fprintf(stderr, "Sum of squared penalty parameters shall not be zero\n");
        exit(-3);
        }
    Sum = sqrt(Sum2);
    if (n == 4 || n == 6)
        {
        nD = n;
        D = new double[n];
        n = 0;
        endptr = S;
        newendptr = 0;
        for (s = S;;s = t + 1,endptr = newendptr)
            {
            d = strtod(s, &newendptr);
            if  (   newendptr && newendptr > s 
                &&  newendptr > endptr 
                &&  (   newendptr[-1] == '.' 
                    ||  (   newendptr[-1] >= '0' 
                        &&  newendptr[-1] <= '9'
                        )
                    )
                )
                {
                D[n] = d / Sum;
                ++n;
                }
            t = strchr(s, ';');
            if(t == 0)
                break;
            }
        
        if(D[1] > 0 && D[2] < 0)
            { // The penalty for W=>R and R=>W must be negative and positive, respectively. (W=>R is "good" and R=>W is "bad")
            for(n = 0;n < nD;++n)
                D[n] = -D[n];
            }
           
        }
    else
        {
        fprintf(stderr, "There must be 4 or 6 penalty parameters, separated by semicolon and no spaces\n");
        exit(-4);
        }
    }

OptReturnTp optionStruct::doSwitch(int optchar, char * locoptarg, char * progname)
    {
    switch (optchar)
        {
        case 'W':
        case 'd':
            printf("Obsolete option %c. Use option -X\n",optchar);
            exit(1);
            break;
        case '@':
            readOptsFromFile(locoptarg, progname);
            break;
        case 'B': // best parms
            B = dupl(locoptarg);
            break;
        case 'P': // current parms
            P = dupl(locoptarg);
            break;
        case 'X':
            X = dupl(locoptarg);
            break;
        case 'c': // pruning threshold, cutoff
            if (locoptarg && *locoptarg)
                c = *locoptarg - '0';

            if (c < 0 || c > 9)
                c = 0;
            break;
        case 'C': // pruning threshold (cutoff) for weightedcount
            if (locoptarg && *locoptarg)
                C = *locoptarg - '0';

            if (C > 9 || C < 0)
                C = -1;
            break;
        case 'e': // extra
            if (locoptarg)
                {
                e = dupl(locoptarg);
                if (strstr(e, "_suffix"))
                    {
                    if (this->SuffixOnlyParmSeen && SuffixOnly == false)
                        {
                        fprintf(stderr, "Option -e [%s] and option -s %s are contradictory\n", e, locoptarg);
                        exit(-1);
                        }
                    SuffixOnly = true;
                    }
                }
            break;
        case 'f': // compfunc
            f = dupl(locoptarg);
            break;
        case 'L':
            Minfraction = strtod(locoptarg, (char**)0);
            if (Minfraction <= 0.0 || 1.0 < Minfraction)
                {
                printf("%s", "Option -L: value must be greater than 0.0 and not greater than 1.0");
                exit(-1);
                }
            break;
        case 'H':
            Maxfraction = strtod(locoptarg, (char**)0);
            if (Maxfraction > 1.0)
                {
                printf("%s", "Option -H: value must be greater than 0.0 and not greater than 1.0");
                exit(-1);
                }
            break;
        case 'K':
            K = strtol(locoptarg,0,10);
            if(K <= 0)
                {
                printf("%s","Option -K: value must be 1 or more");
                exit(-1);
                }
        case 'M':
            M = strtod(locoptarg, (char**)0);
            if (M < 1.0)
                {
                printf("%s", "Option -M: value must be 1.0 or greater");
                exit(-1);
                }
            break;
        case 'N':
            N = strtod(locoptarg, (char**)0);
            if (N < 1.0)
                {
                printf("%s", "Option -N: value must be 1.0 or greater");
                exit(-1);
                }
            break;
        case 'h':
        case '?':
            printf("usage:\n"
                "affixtrain [-@ <option file>] -i <word/lemma list> [-c <pruning threshold>] [-C <expected pruning threshold>] [-o <flexrules>] [-e <extra>] [-n <columns>] [-f <compfunc>] [-p[-]] [-s[-]] [-v[-]] [-j <tempdir>] [-L<n>] [-H<n>]"
                "\nor\n"
                "affixtrain [<word/lemma list> [<pruning threshold> [<flexrules> [<extra> [<columns> [<compfunc>]]]]]]"
                "\nor\n"
                "affixtrain -b <rules>"
                "\n");
            printf("-@: Options are read from file with lines formatted as: -<option letter> <value>\n"
                "    A semicolon comments out the rest of the line.\n"
                );
            printf("-I: (input) full form list. (Together with -b option.)\n");
            printf("-i: (input) full form / lemma list\n");
            printf("-c: discard rules with little support. 0 <= pruning threshold <= 9\n");
            printf("-C: (together with -p and -XW) expected pruning threshold (parameter for tree penalty function). 0 <= expected pruning threshold <= 9\n");
            printf("-D: penalty parameters. Four or six: R__R;W__R;R__W;W__W[;R__NA;W__NA]\n");
            printf("-O: (output) default is stdout (Together with -b and -I options.)\n");
            printf("-o: (output) default is to automatically generate file name\n");
            printf("-e: language code (da,en,is,nl,ru, etc)\n");
            printf("-p: compute parameters (overrules -f)\n");
            printf("-s: create suffix-only rules\n");
            printf("-v: verbose\n");
            printf("-j: directory to store the bulk of intermediate results, also textual presentations of rules. Must not include final (back)slash.\n");
            printf("-L: minimum fraction (with option -f0 or -p)\n");
            printf("-H: maximum fraction (with option -f0 or -p)\n");
            printf("-K: number of differently sized fractions of trainingdata\n");
            printf("-N: number of iterations of training with same fraction of training data when fraction is minimal\n");
            printf("-M: number of iterations of training with same fraction of training data when fraction is maximal\n");
            printf("-X: Tree penalty function C D E or W:\n");
            printf("  C All rules contribute with the same penalty\n");
            printf("  D Rules deeper in the tree contribute more than rules near the root. Or: a rule's contribution is proportional with number of non-wildcard characters in the rule pattern.\n");
            printf("  E Trees evenly distributing examples over the nodes have lowest penalty. (Strive for max entropy)\n");
            printf("  W Rules with one less supporting examples than the expected pruning threshold. Set expected pruning threshold with -C parameter.\n");
            printf("  S Rules with longer patterns contribute more. Or: a rule's contribution is proportional with number of characters in the rule pattern, including wildcards.\n");
            printf("-P: write parameters to file (default parms.txt if -p or -f0, otherwise no parameter file)\n");
            printf("-t: test the rules with data not used for training\n");
            printf("-T: test the rules with the training data\n");
            printf("-F; create flexrules. Can be combined with computation (-p) and testing (-t, -T)\n");
            printf("-R; Redo training with unambiguous training set after each training with\n    ambiguous training set, not using the candidates that would not be\n    lemmatised correctly anyway by the produced rules.\n");
            printf("-n: columns (default FBT if there are three or more columns in the input, otherwise FB):\n");
            printf("  1 or F or f or W or w:Word\n");
            printf("  2 or B or b or L or l:Lemma\n");
            printf("  3 or T or t:POS tag (create rules for each of them)\n");
#if LEMMAINL
            printf("  3:Word Frequency\n");
            printf("  4:Lemma Frequency\n");
#endif
#if WORDCLASS
            printf("  5:Word Class\n");
#endif
#if LEMMACLASS
            printf("  6:Lemma Class\n");
#endif
            printf("  0 or O or o:Other (don't care). (Use this to skip and ignore e.g. column with POS tags.)\n");
            printf("-k: Can be set to a tag name. In that case only lines with that tag name are processed.\n");
            printf("-f: Name or index of comparison function (which look at number of right, wrong and not applicable cases):\n");
            printf("  0:parms (Compute good parameter settings. Use -L and -H if training set must grow.)\n");
#if _NA
            printf("  1:fairly_good\n");
            printf("  2:even_better (Icelandic)\n");
            printf("  3:affiksFEW3 (Slovene)\n");
            printf("  4:affiksFEW\n");
            printf("  5:affiksFEW2 (Dutch, Norwegian)\n");
            printf("  6:fixNA\n");
            printf("  7:fruit\n");
            printf("  8:ice\n");
            printf("  9:pisang\n");
            printf("  10:kiwi\n");
            printf("  11:carrot\n");
            printf("  12:peen (Danish, Polish)\n");
            printf("  13:beet\n");
            printf("  14:sugar (English, Greek, German, Swedish, Russian)\n");
            printf("  15:affiksFEW2org ()\n");
            printf("  16:honey ()\n");
#endif
            printf("  17:koud ()\n");
            printf("  18:parms0 (Use computed parameter settings.)\n");
            printf("  19:parmsoff (obsolete, same as -f18)\n");
            printf("-b: Name of binary rule file. Pretty print rule file and create Bracmat\n    version of rule file. Optionally (-I) lemmatise file.\n");
            printf("-Q: Max recursion depth when attempting to create candidate rule\n");
            printf("-G: External training program\n");
            printf("-E: External lemmatizer program (arg1=input, arg2=rules, arg3=output)\n");
            printf("-VX: 10-fold cross validation\n");
            printf("-x: Keep all intermediary files. (default: -x- delete them)\n");
//            printf("-q: Percentage of training pairs to set aside for testing\n");
            return Leave;
        case 'D':
            detectFloatingPointNumbers(locoptarg);
            break;
        case 'I': // word list
            I = dupl(locoptarg);
            break;
        case 'i': // full form/lemma list
            i = dupl(locoptarg);
            break;
        case 'j': // temp dir
            {
            size_t jlen = strlen(locoptarg);
            if(locoptarg[jlen-1] == '\\' || locoptarg[jlen-1] == '/')
                locoptarg[jlen-1] = '\0'; // Remove separator at end
            j = dupl(locoptarg);
            break;
            }
        case 'n': // columns
            n = dupl(locoptarg);
            break;
        case 'k': // specific POS tag
            k = dupl(locoptarg);
            break;
        case 'O': // lemmas
            O = dupl(locoptarg);
            break;
        case 'o': // flexrules
            o = dupl(locoptarg);
            break;
        case 'p': // compute parms
            ComputeParms = locoptarg && *locoptarg == '-' ? false : true;
            break;
            // GNU >>
        case 'w':
            printf("11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY\n");
            printf("FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN\n");
            printf("OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES\n");
            printf("PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED\n");
            printf("OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF\n");
            printf("MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS\n");
            printf("TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE\n");
            printf("PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,\n");
            printf("REPAIR OR CORRECTION.\n");
            return Leave;
        case 'r':
            printf("12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING\n");
            printf("WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR\n");
            printf("REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,\n");
            printf("INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING\n");
            printf("OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED\n");
            printf("TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY\n");
            printf("YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER\n");
            printf("PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE\n");
            printf("POSSIBILITY OF SUCH DAMAGES.\n");
            return Leave;
            // << GNU
        case 's': // create suffix-only rules
            if (SuffixOnly) // If earlier -e option has "_suffix" suffix
                {
                if (locoptarg && *locoptarg == '-')
                    {
                    fprintf(stderr, "Option -e [%s] and option -s %s are contradictory\n", e, locoptarg);
                    exit(-1);
                    }
                }
            SuffixOnly = locoptarg && *locoptarg == '-' ? false : true;
            SuffixOnlyParmSeen = true;
            break;
        case 'A':
            ExpensiveInfix = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'v': // verbose
            if (locoptarg && *locoptarg)
                if (*locoptarg == '-')
                    Verbose = 0;
                else
                    {
                    Verbose = strtol(locoptarg, (char**)0, 10);
                    if (Verbose < 0)
                        {
                        fprintf(stderr, "Option v:Invalid value [%d]. Verbosity can be - or a natural number.\n", Verbose);
                        exit(-1);
                        }
                    }
            else
                Verbose = INT_MAX;
            break;
        case 'x':
            Remove = locoptarg && *locoptarg == '-' ? true : false;
            break;
        case 'T': // test with the training data
            TrainTest = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 't': // test with data that is not use for training
            Test = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'F': // create flexrules
            F = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'R':
            Redo = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'b': // raw rules
            b = dupl(locoptarg);
            break;
        case 'Q': // max recursion depth when attempting to create candidate rule
            if (locoptarg && *locoptarg)
                {
                Q = strtol(locoptarg, (char**)0, 10);
                if (Q < 1)
                    {
                    fprintf(stderr, "Option Q:Invalid value [%d]. Max recursion depth when attempting to create candidate rule must be >= 1\n", Q);
                    exit(-1);
                    }
                }
            else
                {
                fprintf(stderr, "Option Q:No parameter value found.\n");
                exit(-1);
                }
            break;
        case 'E': // External training program
            E = dupl(locoptarg);
            break;
        case 'G': // External lemmatizer
            G = dupl(locoptarg);
            break;
        case 'V': // Tenfold cross-validation
            VX = locoptarg && *locoptarg == 'X';
            break;
            /*
        case 'q':
            if (locoptarg && *locoptarg)
                {
                q = strtol(locoptarg, (char**)0, 10);
                if (q < 0 || 100 < q)
                    {
                    fprintf(stderr, "Option q:Invalid value [%d]. Percentage of training pairs to set aside for testing must be 0 <= q <= 100\n", q);
                    exit(-1);
                    }
                }
            else
                {
                fprintf(stderr, "Option q:No parameter value found.\n");
                exit(-1);
                }
                */
        }
    return GoOn;
    }



OptReturnTp optionStruct::readOptsFromFile(char * locoptarg, char * progname)
    {
    char ** poptions;
    char * options;
    FILE * fpopt = fopen(locoptarg, "r");
    ++openfiles;
    OptReturnTp result = GoOn;
    if (fpopt)
        {
        char * p;
        char line[1000];
        int lineno = 0;
        size_t bufsize = 0;
        while (fgets(line, sizeof(line) - 1, fpopt))
            {
            lineno++;
            size_t off = strspn(line, " \t");
            if (line[off] == ';')
                continue; // comment line
            if (line[off] == '-')
                {
                off++;
                if (line[off])
                    {
                    char * optarg2 = line + off + 1;
                    size_t off2 = strspn(optarg2, " \t");
                    if (!optarg2[off2])
                        optarg2 = NULL;
                    else
                        optarg2 += off2;
                    if (optarg2)
                        {
                        for (p = optarg2 + strlen(optarg2) - 1; p >= optarg2; --p)
                            {
                            if (!isspace((unsigned char)*p))
                                break;
                            *p = '\0';
                            }
                        bool string = false;
                        if (*optarg2 == '\'' || *optarg2 == '"')
                            {
                            for (p = optarg2 + strlen(optarg2) - 1; p > optarg2; --p)
                                {
                                if (*p == *optarg2)
                                    {
                                    string = true;
                                    for (char * q = p + 1; *q; ++q)
                                        {
                                        if (*q == ';')
                                            break;
                                        if (!isspace((unsigned char)*q))
                                            {
                                            string = false;
                                            }
                                        }
                                    if (string)
                                        {
                                        *p = '\0';
                                        ++optarg2;
                                        }
                                    break;
                                    }
                                }
                            }
                        if (!*optarg2 && !string)
                            optarg2 = NULL;
                        }
                    if (optarg2)
                        {
                        bufsize += strlen(optarg2) + 1;
                        }
                    else
                        ++bufsize;
                    char * optpos = strchr(opts, line[off]);
                    if (optpos)
                        {
                        if (optpos[1] != ':')
                            {
                            if (optarg2)
                                {
                                printf("Option argument %s provided for option letter %c that doesn't use it on line %d in option file \"%s\"\n", optarg2, line[off], lineno, locoptarg);
                                exit(1);
                                }
                            }
                        }
                    }
                else
                    {
                    printf("Missing option letter on line %d in option file \"%s\"\n", lineno, locoptarg);
                    exit(1);
                    }
                }
            }
        rewind(fpopt);

        poptions = new char *[lineno];
        options = new char[bufsize];
        // update stacks that keep pointers to the allocated arrays.
        optionSets++;
        char *** tmpPpoptions = new char **[optionSets];
        char ** tmpPoptions = new char *[optionSets];
        int g;
        for (g = 0; g < optionSets - 1; ++g)
            {
            tmpPpoptions[g] = Ppoptions[g];
            tmpPoptions[g] = Poptions[g];
            }
        tmpPpoptions[g] = poptions;
        tmpPoptions[g] = options;
        delete[] Ppoptions;
        Ppoptions = tmpPpoptions;
        delete[] Poptions;
        Poptions = tmpPoptions;

        lineno = 0;
        bufsize = 0;
        while (fgets(line, sizeof(line) - 1, fpopt))
            {
            poptions[lineno] = options + bufsize;
            size_t off = strspn(line, " \t");
            if (line[off] == ';')
                continue; // comment line
            if (line[off] == '-')
                {
                off++;
                if (line[off])
                    {
                    char * optarg2 = line + off + 1;
                    size_t off2 = strspn(optarg2, " \t");
                    if (!optarg2[off2])
                        optarg2 = NULL;
                    else
                        optarg2 += off2;
                    if (optarg2)
                        {
                        for (p = optarg2 + strlen(optarg2) - 1; p >= optarg2; --p)
                            {
                            if (!isspace((unsigned char)*p))
                                break;
                            *p = '\0';
                            }
                        bool string = false;
                        if (*optarg2 == '\'' || *optarg2 == '"')
                            {
                            for (p = optarg2 + strlen(optarg2) - 1; p > optarg2; --p)
                                {
                                if (*p == *optarg2)
                                    {
                                    string = true;
                                    for (char * q = p + 1; *q; ++q)
                                        {
                                        if (*q == ';')
                                            break;
                                        if (!isspace((unsigned char)*q))
                                            {
                                            string = false;
                                            }
                                        }
                                    if (string)
                                        {
                                        *p = '\0';
                                        ++optarg2;
                                        }
                                    break;
                                    }
                                }
                            }
                        if (!*optarg2 && !string)
                            optarg2 = NULL;
                        }
                    if (optarg2)
                        {
                        strcpy(poptions[lineno], optarg2);
                        bufsize += strlen(optarg2) + 1;
                        }
                    else
                        {
                        poptions[lineno][0] = 0;
                        ++bufsize;
                        }
                    OptReturnTp res = doSwitch(line[off], poptions[lineno], progname);
                    if (res > result)
                        result = res;
                    }
                }
            lineno++;
            }
        --openfiles;
        fclose(fpopt);
        }
    else
        {
        printf("Cannot open option file %s\n", locoptarg);
        }
    return result;
    }

void optionStruct::completeArgs()
    {
    if(b != 0)
        return;

    if(!j)
        {
        j = dupl("tmp");
        }

    if (SuffixOnly)
        {
        if (e)
            {
            if (!strstr(e, "_suffix"))
                {
                static char * ext = new char[strlen(e) + strlen("_suffix") + 1];
                sprintf(ext, "%s_suffix", e);
                delete[] e;
                e = ext;
                }
            }
        else
            {
            e = dupl("xx_suffix");
            }
        }
    else if (ExpensiveInfix)
        {
        if (e)
            {
            if (!strstr(e, "_fewinfix"))
                {
                static char * ext = new char[strlen(e) + strlen("_fewinfix") + 1];
                sprintf(ext, "%s_fewinfix", e);
                delete[] e;
                e = ext;
                }
            }
        else
            {
            e = dupl("xx_fewinfix");
            }
        }

    if (P == 0)
        {
        if (ComputeParms)
            {
            if (f == NULL)
                {
                if (Verbose)
                    printf("computeParms == true\n");
                P = dupl("parms.txt");
                }
            else
                {
                printf("Option -p should not be set if option -f is set.\n");
                exit(1);
                }
            }
        }

    if (P)
        {
        if (!B)
            {
            char * tmp;
            const char * best = "best";  
            const char * exte = ".txt";  
            if(e)
                {
                tmp = new char[strlen(best) + 1 + strlen(e) + strlen(exte) + 1];
                sprintf(tmp, "%s_%s%s", best, e, exte);
                }
            else
                {
                tmp = new char[strlen(best) + strlen(exte) + 1];
                sprintf(tmp, "%s%s", best, exte);
                }
            B = tmp;
            }
        }

    if(X)
        {
        if(!strcmp(X,"W"))
            {
            WeightFunction = esupport;
            if(C < 0)
                C = 0;
            }
        else if(!strcmp(X,"D"))
            {
            WeightFunction = edepth;
            }
        else if(!strcmp(X,"S"))
            {
            WeightFunction = esize;
            }
        else if(!strcmp(X,"C"))
            {
            WeightFunction = econstant;
            }
        else if (!strcmp(X, "E"))
            {
            WeightFunction = eentropy;
            }
        else
            {
            printf("Option -X can only have values C, D, E or W\n");
            exit(1);
            }
        }

    if(WeightFunction != esupport)
        C = -1; // This option is only used in combination with -XW.

    if (!i)
        {
        fprintf(stderr, "Error: No training data (-i option)\n");
        exit(2);
        }

    //const char * columns = "12634"; // Word, Lemma, LemmaClass, WordFreq, LemmaFreq
    if (!n)
        {
        if (k && k[0])
            n = dupl("FBT");// Word, Lemma, TAG
//        else
//            n = dupl("FBO");// Word, Lemma, Other
        }

    if (!f)
        {
        if (!ComputeParms && !D)
            {
            fprintf(stderr, "Error: No parameters provided. (-D option).\n"
                "Add -p if you want to compute the parameters.\n"
                "If the parameters are already computed, edit the generated parameter file\n"
                "(e.g., parms.txt) and then run the program with\n\n\taffixtrain -@ parms.txt\n"
                );
            exit(2);
            }
        }

    if (ComputeParms)
        {
        if(M > N)
            {
            fprintf(stderr, "M (%f) must be <= N (%f)\n",M,N);
            exit(2);
            }

        if (Minfraction > Maxfraction)
            {
            fprintf(stderr, "L (%f) must be <= H (%f)\n",Minfraction ,Maxfraction);
            exit(2);
            }
        if(nD == 0)
            {
            detectFloatingPointNumbers("0;-0.7;0.7;-0.1;0.1;0");
            //detectFloatingPointNumbers("0.088468;-0.685320;0.711931;-0.083827;0.090956;0.019069");
            }
        }


    FILE * fpWrdLem = fopen(wordLemmaList(), "rb");
    if(!fpWrdLem)
        {
        fprintf(stderr,"Input %s could not be opened for reading (-i parameter)\n",wordLemmaList());
        exit(-1);
        }
    ++openfiles;
    Blobs = 0; // If there are no non-empty lines, there are no blobs either.
    Lines = 0;
    int kar = 0;
    unsigned int lineProfile = 0;
    while ((kar = fgetc(fpWrdLem)) != EOF)
        {
        if (kar == '\n')
            {
            if(lineProfile & 4)
                ++Lines;

            if(lineProfile > 1)
                { // Just finished non-empty line
                lineProfile >>= 1; // push line before previous line out of
                                   // memory
                // if lineProfile is 1, we may have a bunch of empty lines
                // before a new non-empty line appears
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
                    ++Blobs;
                    }
                lineProfile |= 4; // set third bit when making non-empty line
                }
            else
                ;
            }
        }
    if(lineProfile & 4)
        { // last line was not finished off with newline
        ++Lines;
        }
    if(Lines > 0)
        ++Blobs;

    if (verbose() > 5)
        printf("blobs:%d lines %d\n", Blobs, Lines);
    --openfiles;
    fclose(fpWrdLem);
    if(Minfraction < MINLINES/(double)Lines)
        Minfraction = MINLINES/(double)Lines;
    if(Minfraction > 1.0)
        Minfraction = 1.0;
    if(Maxfraction < Minfraction)
        Maxfraction = Minfraction;
    }


OptReturnTp optionStruct::readArgs(int argc, char * argv[])
    {
    int optchar;

    OptReturnTp result = GoOn;
    while ((optchar = getopt(argc, argv, opts)) != -1 && result != Leave)
        {
        OptReturnTp res = doSwitch(optchar, optarg, argv[0]);
        if (res > result)
            result = res;
        }
    if(result != Leave)
        {
        if (!i && !c && !o && !e && !n)
            {
            while (optind < argc)
                {
                if (!i)
                    i = dupl(argv[optind++]);
                else if (c < -1)
                    {
                    if (argv[optind] && *argv[optind])
                        c = *argv[optind] - '0';

                    if (c > 9 || c < 0)
                        c = -1;
                    }
                else if (!o)
                    o = dupl(argv[optind++]);
                else if (!e)
                    e = dupl(argv[optind++]);
                else if (!n)
                    n = dupl(argv[optind++]);
                else
                    printf("Too many arguments:%s\n", argv[optind]);
                }
            }
        else if (optind < argc)
            {
            if (i && c && o && e && n)
                printf("Too many arguments:%s\n", argv[optind]);
            else
                printf("You cannot have a command line with both option-style arguments and option-less-fixed-position arguments:%s\n", argv[optind]);
            }

        completeArgs();
        setArgstring();
        }

    return result;
    }

void optionStruct::print(FILE * fp) const
    {
    fprintf(fp, "               ; verbose (-v <n>: yes (n=1 high priority n > 1 lower priority) -v-: no)\n-v %s\n", Verbose ? "" : "-");
    fprintf(fp, "               ; keep intermediary files (-x: yes -x-: no)\n-x %s\n", Remove ? "-" : "");
    fprintf(fp, "               ; 10-fold cross validation (-VX: yes, overrules T and t options ;-VX: no)\n%s-VX\n", VX ? "" : ";");
    fprintf(fp, "               ; External training program\n%s-G%s\n", G ? "" : ";", G ? G : "<program name>");
    fprintf(fp, "               ; External lemmatizer\n%s-E%s\n", E ? "" : ";", E ? E : "<program name>");
    if (b)
        {
        fprintf(fp, "               ; flex rules (input file, binary format)\n-b %s\n-b %s\n", b, b);
        fprintf(fp, "               ; Word list (Optional if -b is specified. Otherwise N/A) (-I filename)\n");
        if(I)
            {
            fprintf(fp, "-I %s\n", I);
            }
        else
            {
            fprintf(fp, ";-I (Not specified, no lemmatisation takes place)\n");
            }
        fprintf(fp, "               ; Output, lemmas of words in input (-I option)\n");
        if(I)
            {
            if(O)
                fprintf(fp, "-O %s\n", O);
            else
                fprintf(fp, ";-O (Not specified, defaults to stdout)\n");
            }
        else
            fprintf(fp, ";-O (N/A)\n");
        fprintf(fp, "               ; word/lemma list\n;-i %s (N/A)\n", i ? i : "");
        fprintf(fp, "               ; extra file name affix\n;-e %s (N/A)\n", e ? e : "");
        fprintf(fp, "               ; suffix only (-s: yes -s-: no)\n;-s %s (N/A)\n", SuffixOnly ? "" : "-");
        fprintf(fp, "               ; make rules with infixes less prevalent(-A: yes -A-: no)\n;-A %s (N/A)\n", ExpensiveInfix ? "" : "-");
        fprintf(fp, "               ; columns (1=word,2=lemma,3=tags,0=other)\n;-n %s (N/A)\n", n ? n : "");
        fprintf(fp, "               ; specific POS tag (default: empty string)\n;-k %s (N/A)\n", k ? k : "");
        fprintf(fp, "               ; max recursion depth when attempting to create candidate rule\n;-Q %d (N/A)\n", Q);
        fprintf(fp, "               ; flex rules (output, binary format)\n;-o %s (N/A)\n", o ? o : "");
        fprintf(fp, "               ; temp dir (including separator at end!)\n;-j %s (N/A)\n", j ? j : "");
//        fprintf(fp, "               ; percentage of training pairs to set aside for testing\n;-q %d (N/A)\n", q);
        fprintf(fp, "               ; penalties to decide which rule survives (4 or 6 floating point numbers: R=>R;W=>R;R=>W;W=>W[;R=>N/A;W=>NA], where R=#right cases, W=#wrong cases, N/A=#not applicable cases, previous success state=>success state after rule application)\n"); if (nD > 0){ fprintf(fp, ";-D "); for (int i = 0; i < nD; ++i)fprintf(fp, "%.10f;", D[i]); fprintf(fp, " (N/A)\n"); } else fprintf(fp, ";-D (N/A)\n");
        fprintf(fp, "               ; compute parms (-p: yes -p-: no)\n;-p %s (N/A)\n", ComputeParms ? "" : "-");
        fprintf(fp, "               ; expected optimal pruning threshold (only effective in combination with -XW)\n;-C %d (N/A)\n", C);
        fprintf(fp, "               ; tree penalty (-XC: constant -XD: more support is better -XE: higher entropy is better -XW: Fewer pattern characters other than wildcards is better)\n;-X %s (N/A)\n", X);
        fprintf(fp, "               ; current parameters (-P filename)\n;-P %s (N/A)\n", P ? P : "");
        fprintf(fp, "               ; best parameters (-B filename)\n;-B %s (N/A)\n", B ? B : "");
        fprintf(fp, "               ; start training with minimal fraction of training pairs (-Ln: 0.0 < n <= 1.0)\n;-L %f (N/A)\n", Minfraction);
        fprintf(fp, "               ; end training with maximal fraction of training pairs (-Hn: 0.0 < n <= 1.0)\n;-H %f (N/A)\n", Maxfraction);
        fprintf(fp, "               ; number of differently sized fractions of trainingdata (natural number)\n;-K %d (N/A)\n",K);
        fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is minimal (positive number)\n;-N %f (N/A)\n", N);
        fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is maximal (positive number)\n;-M %f (N/A)\n", M);
        fprintf(fp, "               ; competition function (deprecated)\n;-f %s (N/A)\n", f ? f : "");
        fprintf(fp, "               ; redo training after homographs for next round are removed (-R: yes -R-: no)\n;-R %s (N/A)\n", Redo ? "" : "-");
        fprintf(fp, "               ; max. pruning threshold to evaluate\n;-c %d (N/A)\n", c);
        fprintf(fp, "               ; test with the training data (-T: yes -T-: no)\n;-T %s (N/A)\n", TrainTest ? "" : "-");
        fprintf(fp, "               ; test with data not used for training (-t: yes -t-: no)\n;-t %s (N/A)\n", Test ? "" : "-");
        fprintf(fp, "               ; create flexrules using full training set (-F: yes -F-: no)\n;-F %s (N/A)\n", F ? "" : "-");
        }
    else
        {
        fprintf(fp, "               ; flex rules (input file, binary format)\n;-b (not specified, not doing actions with already created flex rules.)\n");
        fprintf(fp, "               ; Word list (Optional if -b is specified. Otherwise N/A) (-I filename)\n;-I %s (N/A)\n",I ? I : "");
        fprintf(fp, "               ; Output, lemmas of words in input (-I option)\n;-O %s (N/A)\n",O ? O : "");
        fprintf(fp, "               ; word/lemma list\n-i %s\n", i);
        fprintf(fp, "               ; extra file name affix\n-e %s\n", e ? e : "");
        fprintf(fp, "               ; suffix only (-s: yes -s-: no)\n-s %s\n", SuffixOnly ? "" : "-");
        fprintf(fp, "               ; make rules with infixes less prevalent(-A: yes -A-: no)\n-A %s\n", ExpensiveInfix ? "" : "-");
        fprintf(fp, "               ; columns (1 or F or W=word,2 or B or L=lemma,3 or T=tags,0 or O=other)\n-n %s\n", n);
        fprintf(fp, "               ; specific POS tag (default: empty string)\n-k %s\n", k);
        fprintf(fp, "               ; max recursion depth when attempting to create candidate rule\n-Q %d\n", Q);
        fprintf(fp, "               ; flex rules (output, binary format, can be left unspecified)\n%s-o %s\n",o ? "" : ";", o ? o : "(Not specified, autogenerated)");
        fprintf(fp, "               ; temp dir\n-j %s\n", j ? j : "");
//        fprintf(fp, "               ; percentage of training pairs to set aside for testing\n-q %d\n", q);
        fprintf(fp, "               ; penalties to decide which rule survives (4 or 6 floating point numbers: R=>R;W=>R;R=>W;W=>W[;R=>N/A;W=>NA], where R=#right cases, W=#wrong cases, N/A=#not applicable cases, previous success state=>success state after rule application)\n"); if (nD > 0){ fprintf(fp, "-D "); for (int i = 0; i < nD; ++i)fprintf(fp, "%.10f;", D[i]); fprintf(fp, "\n"); } else fprintf(fp, ";-D not specified\n");
        fprintf(fp, "               ; compute parms (-p: yes -p-: no)\n-p %s\n", ComputeParms ? "" : "-");
        if (ComputeParms)
            {
            assert(P);
            assert(B);
            fprintf(fp, "               ; expected optimal pruning threshold (only effective in combination with -XW)\n-C %d\n", C);
            fprintf(fp, "               ; tree penalty (-XC: constant -XD: more support is better -XE: higher entropy is better -XW: Fewer pattern characters other than wildcards is better)\n-X %s\n", X ? X : "C");
            fprintf(fp, "               ; current parameters (-P filename)\n-P %s\n", P);
            fprintf(fp, "               ; best parameters (-B filename)\n-B %s\n", B);
            fprintf(fp, "               ; start training with minimal fraction of training pairs (-Ln: 0.0 < n <= 1.0)\n-L %f\n", Minfraction);
            fprintf(fp, "               ; end training with maximal fraction of training pairs (-Hn: 0.0 < n <= 1.0)\n-H %f\n", Maxfraction);
            fprintf(fp, "               ; number of differently sized fractions of trainingdata (natural number)\n-K %d\n",K);
            fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is minimal (positive number)\n-N %f\n", N);
            fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is maximal (positive number)\n-M %f\n", M);
            fprintf(fp, "               ; competition function (deprecated)\n;-f %s (N/A)\n", f ? f : "");
            fprintf(fp, "               ; redo training after homographs for next round are removed (-R: yes -R-: no)\n;-R %s (N/A)\n", Redo ? "" : "-");
            }
        else
            {
//            assert(f);
            fprintf(fp, "               ; expected optimal pruning threshold (only effective in combination with -XW)\n;-C %d (N/A)\n", C);
            fprintf(fp, "               ; tree penalty (-XC: constant -XD: more support is better -XE: higher entropy is better -XW: Fewer pattern characters other than wildcards is better)\n;-X %s (N/A)\n", X ? X : "C");
            fprintf(fp, "               ; current parameters (-P filename)\n;-P %s (N/A)\n", P ? P : "");
            fprintf(fp, "               ; best parameters (-B filename)\n;-B %s (N/A)\n", B ? B : "");
            fprintf(fp, "               ; start training with minimal fraction of training pairs (-Ln: 0.0 < n <= 1.0)\n;-L %f (N/A)\n", Minfraction);
            fprintf(fp, "               ; end training with maximal fraction of training pairs (-Hn: 0.0 < n <= 1.0)\n;-H %f (N/A)\n", Maxfraction);
            fprintf(fp, "               ; number of differently sized fractions of trainingdata (natural number)\n;-K %d (N/A)\n",K);
            fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is minimal (positive number)\n;-N %f (N/A)\n", N);
            fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is maximal (positive number)\n;-M %f (N/A)\n", M);
            fprintf(fp, "               ; competition function (deprecated)\n%s-f %s\n",f ? "" : ";", f ? f : "");
            fprintf(fp, "               ; redo training after homographs for next round are removed (-R: yes -R-: no)\n-R %s\n", Redo ? "" : "-");
            }
        if(TrainTest || Test)
            fprintf(fp, "               ; max. pruning threshold to evaluate\n-c %d\n", c);
        else
            fprintf(fp, "               ; max. pruning threshold to evaluate\n;-c %d (N/A)\n", c);
        fprintf(fp, "               ; test with the training data (-T: yes -T-: no)\n-T %s\n", TrainTest ? "" : "-");
        fprintf(fp, "               ; test with data not used for training (-t: yes -t-: no)\n-t %s\n", Test ? "" : "-");
        fprintf(fp, "               ; create flexrules using full training set (-F: yes -F-: no)\n-F %s\n", F ? "" : "-");
        fprintf(fp, "               ; Number of clusters found in word/lemma list: %d\n", Blobs);
        fprintf(fp, "               ; Number of lines found in word/lemma list:    %d\n", Lines);
        }
    }

void optionStruct::seti(const char * WordLemmaList)
    {
    delete i;
    i = dupl(WordLemmaList);
    }

void optionStruct::setI(const char * WordList)
    {
    delete I;
    I = dupl(WordList);
    }

void optionStruct::seto(const char * Result)
    {
    delete o;
    o = dupl(Result);
    }

void optionStruct::setO(const char * Result)
    {
    delete O;
    O = dupl(Result);
    }

void optionStruct::sete(const char * Extra)
    {
    delete e;
    e = dupl(Extra);
    }

void optionStruct::setn(const char * Columns)
    {
    delete n;
    n = dupl(Columns);
    }

void optionStruct::setk(const char * PoS)
    {
    delete k;
    k = dupl(PoS);
    }

void optionStruct::setf(const char * Compfunc)
    {
    delete f;
    f = dupl(Compfunc);
    }

void optionStruct::setP(const char * ParamFile)
    {
    delete P;
    P = dupl(ParamFile);
    }

void optionStruct::setArgstring()
    {
    if(!i)
        return;

    if(Argstring)
        delete[] Argstring;

    size_t nameLength;
    if(G)
        {
        nameLength = strlen(i) + (e ? 1 + strlen(e) : 0) + (1+strlen(G)) + strlen("_externalTrainer")+1+(VX ? strlen("_VX"):0);
        
        Argstring = new char[nameLength];
        strcpy(Argstring, i);
        if (e)
            {
            strcat(Argstring, "_");
            strcat(Argstring, e);
            }
        strcat(Argstring, "_");
        strcat(Argstring, G);
        strcat(Argstring, "_externalTrainer");
        if(VX)
            strcat(Argstring, "_VX");
        }
    else
        {
        nameLength = strlen(i) + (e ? 1 + strlen(e) : 0) + (SuffixOnly ? strlen("_suf") : 0) + (ExpensiveInfix ? strlen("_inf") : 0) + (C < 0 ? 0 : strlen("_C") + 1) + (WeightFunction == econstant ? strlen("_XE") : 0) + (WeightFunction == esupport ? strlen("_XW") : 0) + (WeightFunction == eentropy ? strlen("_XE") : 0) + (WeightFunction == edepth ? strlen("_XD") : 0) + (WeightFunction == esize ? strlen("_XS") : 0) + (Redo ? strlen("_R") : 0) /*+ (k && k[0] ? 1 + strlen(k) : 0)*/ + 1 + (VX ? strlen("_VX") : 0);
        
        Argstring = new char[nameLength];
        strcpy(Argstring, i);
        if (e)
            {
            strcat(Argstring, "_");
            strcat(Argstring, e);
            }
        if (SuffixOnly)
            {
            strcat(Argstring, "_suf");
            }
/*        if (k && k[0])
            {
            strcat(Argstring, "_");
            strcat(Argstring, k);
            }*/
        if (ExpensiveInfix)
            {
            strcat(Argstring, "_inf");
            }
        if (C >= 0)
            {
            strcat(Argstring, "_C");
            size_t L = strlen(Argstring);
            Argstring[L] = (char)(C + '0');
            Argstring[L + 1] = 0;
            }
        switch(WeightFunction)
            {
            case econstant:
                strcat(Argstring, "_XC");
                break;
            case esupport:
                strcat(Argstring, "_XW");
                break;
            case eentropy:
                strcat(Argstring, "_XE");
                break;
            case edepth:
                strcat(Argstring, "_XD");
                break;
            case esize:
                strcat(Argstring, "_XS");
                break;
            default:
                break;
            }
        if (Redo)
            {
            strcat(Argstring, "_R");
            }
        if(VX)
            strcat(Argstring, "_VX");
        }
    }

const char * optionStruct::argstring_no_path() const
    {
    assert(Argstring);
    const char * fw_slash = strrchr(Argstring,'/');
    const char * back_slash = strrchr(Argstring,'\\');
    if(fw_slash != 0)
        {
        if(back_slash != 0)
            {
            if(back_slash < fw_slash)
                return fw_slash + 1;
            else
                return back_slash + 1;
            }
        else
            return fw_slash + 1;
        }
    else if(back_slash != 0)
        return back_slash + 1;
    else
        return Argstring;
    }


void optionStruct::printArgFile() const
    {
    const char * Args = argstring_no_path();
    char * name = new char[strlen("parms.")+strlen(Args)+1];
    strcpy(name, "parms.");
    strcat(name, Args);
    FILE * fp = fopen(name, "wb");
    ++openfiles;
    print(fp);
    --openfiles;
    fclose(fp);
    }


void optionStruct::printEvaluation(const char * introduction,char * evaluation,char * postScriptum) const
    {
    if(evaluation)
        {
        const char * Args = argstring_no_path();
        char * name = new char[strlen("parms.")+strlen(Args)+1];
        strcpy(name, "parms.");
        strcat(name, Args);
        FILE * fp = fopen(name, "ab");
        ++openfiles;
        fprintf(fp,"\n; Evaluation:\n; -----------\n%s\n%s",introduction,evaluation);
        if(postScriptum && *postScriptum)
            {
            fprintf(fp,"\n; Postscriptum\n\n%s",postScriptum);
            }
        --openfiles;
        fclose(fp);
        }
    }


const char * optionStruct::flexrules()
    {
    if (o == NULL)
        {
        const char * Args = argstring_no_path();
        char * name = new char[strlen("flexrules.")+strlen(Args)+1];
        strcpy(name, "flexrules.");
        strcat(name, Args);
        o = name;
        }
    return o;
    }
