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


//        printf("usage: makeaffixrules -w <word list> -c <cutoff> -C <expected cutoff> -o <flexrules> -e <extra> -n <columns> -f <compfunc> [<word list> [<cutoff> [<flexrules> [<extra> [<columns> [<compfunc>]]]]]]\n");

static char opts[] = "?@:B:b:c:C:D:e:f:hH:i:j:K:L:M:N:n:o:O:p:P:Q:q:R:s:T:t:v:W:" /* GNU: */ "wr";
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
    i = dupl(O.i);
    n = dupl(O.n);
    o = dupl(O.o);
    B = dupl(O.B);
    P = dupl(O.P);
    j = dupl(O.j);
    b = dupl(O.b);
    nD = O.nD;
    D = new double[nD];
    for(int ii = 0;ii < nD;++ii)
        D[ii] = O.D[ii];
    ComputeParms = O.ComputeParms;
    SuffixOnly = O.SuffixOnly;
    SuffixOnlyParmSeen = O.SuffixOnlyParmSeen;
    Verbose = O.Verbose;
    Minfraction = O.Minfraction;
    Maxfraction = O.Maxfraction;
    Doweights = O.Doweights;
    Redo = O.Redo;
    Test = O.Test;
    TrainTest = O.TrainTest;
    Q = O.Q;
    q = O.q;
    K = O.K;
    M = O.M;
    N = O.N;

    Blobs = O.Blobs;
    Lines = O.Lines;
    FracBlobs = O.FracBlobs;
    FracLines = O.FracLines;
    }

optionStruct::optionStruct()
    {
    c = -1; // cutoff 
    C = -1; // expected cutoff when determining the parameters using weightedcount()
    // Rules with C+1 supporting word/lemma pairs have the highest penalty.
    // Rules with C or fewer supporting word/lemma pairs are probably best cut off.
    // To decide whether that is the case, the rules must be tested with OOV word/lemma pairs.
    e = NULL; // extra
    f = NULL; // compfunc
    i = NULL; // word list
    n = NULL; // columns
    o = NULL; // flexrules
    B = NULL;
    P = NULL;
    j = NULL; // temp dir
    b = NULL; // raw file (see t)
    D = NULL;
    nD = 0;
    ComputeParms = false;// compute parms
    SuffixOnly = false;// suffix rules only
    SuffixOnlyParmSeen = false;
    Verbose = false;// verbose
    Minfraction = 1.0; // L
    Maxfraction = 1.0; // H
    Doweights = false;
    Redo = false;
    Test = false;
    TrainTest = false;
    Q = 1;
    q = 1;
    K = 20;   // Number of differently sized fractions of trainingdata 
    M = 10.0; // # Iterations when training with Maxfraction of input
    N = 100.0;// # Iterations when training with Minfraction of input

    Blobs = 0; // Number of blobs found in word list
    Lines = 0; // Number of lines found in word list
    FracBlobs = 0; // Number of blobs used for training
    FracLines = 0; // Number of lines used for training
    }

void cleanUpOptions()
    {
    for (int I = 0; I < optionSets; ++I)
        {
        delete[] Poptions[I];
        delete[] Ppoptions[I];
        }
    delete[] Poptions;
    delete[] Ppoptions;
    }

optionStruct::~optionStruct()
    {
    delete[] e;
    delete[] f;
    delete[] i;
    delete[] j;
    delete[] n;
    delete[] o;
    delete[] B;
    delete[] P;
    delete[] b;
    delete[] D;
    }

void optionStruct::detectFloatingPointNumbers(char * S)
    {
    int n = 0;
    char * t = 0;
    double Sum = 0.0;
    double d;
    char * s;
    char * endptr = 0;
    for (s = S; (t = strchr(s, ';')) != 0; s = t + 1)
        {
        d = strtod(s, &endptr);
        if (endptr && endptr > s && (endptr[-1] == '.' || (endptr[-1] >= '0' && endptr[-1] <= '9')))
            {
            if (d >= 0.0)
                Sum += d;
            else
                Sum -= d;
            ++n;
            }
        }
    d = strtod(s, &endptr);
    if (endptr && endptr > s && (endptr[-1] == '.' || (endptr[-1] >= '0' && endptr[-1] <= '9')))
        {
        if (d >= 0.0)
            Sum += d;
        else
            Sum -= d;
        ++n;
        }
    if (Sum <= 0)
        {
        fprintf(stderr, "Sum of penalty parameters shall not be zero\n");
        exit(-3);
        }
    if (n == 4 || n == 6)
        {
        nD = n;
        D = new double[n];
        n = 0;
        for (s = S; (t = strchr(s, ';')) != 0; s = t + 1)
            {
            d = strtod(s, &endptr);
            if (endptr && endptr > s && (endptr[-1] == '.' || (endptr[-1] >= '0' && endptr[-1] <= '9')))
                D[n++] = d / Sum;
            }
        d = strtod(s, &endptr);
        if (endptr && endptr > s && (endptr[-1] == '.' || (endptr[-1] >= '0' && endptr[-1] <= '9')))
            D[n] = d / Sum;
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
        case '@':
            readOptsFromFile(locoptarg, progname);
            break;
        case 'B': // best parms
            B = dupl(locoptarg);
            break;
        case 'P': // current parms
            P = dupl(locoptarg);
            break;
        case 'c': // cutoff
            if (locoptarg && *locoptarg)
                c = *locoptarg - '0';

            if (c > 9 || c < 0)
                c = -1;
            break;
        case 'C': // cutoff for weightedcount
            if (locoptarg && *locoptarg)
                C = *locoptarg - '0';

            if (C > 9 || C < 0)
                C = -1;
            break;
        case 'e': // extra
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
                "affixtrain [-@ <option file>] -i <word list> [-c <cutoff>] [-C <expected cutoff>] [-o <flexrules>] [-e <extra>] [-n <columns>] [-f <compfunc>] [-p[-]] [-s[-]] [-v[-]] [-j <tempdir>] [-L<n>] [-H<n>]"
                "\nor\n"
                "affixtrain [<word list> [<cutoff> [<flexrules> [<extra> [<columns> [<compfunc>]]]]]]"
                "\nor\n"
                "affixtrain -b <rules>"
                "\n");
            printf("-@: Options are read from file with lines formatted as: -<option letter> <value>\n"
                "    A semicolon comments out the rest of the line.\n"
                );
            printf("-i: (input) full form / lemma list\n");
            printf("-c: discard rules with little support. 0 <= cutoff <= 9\n");
            printf("-C: (together with -p) expected cutoff (parameter for rule weight function). 0 <= cutoff <= 9\n");
            printf("-D: penalty parameters. Four or six: R__R;W__R;R__W;W__W[;R__NA;W__NA]\n");
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
            printf("-W: minimise weight, not count (sum of rules) (with -p or -f0)\n");
            printf("-P: write parameters to file (default parms.txt if -p or -f0, otherwise no parameter file)\n");
            printf("-t: test the rules\n");
            printf("-n: columns (default 120):\n");
            printf("  1:Word\n");
            printf("  F:Word\n");
            printf("  f:Word\n");
            printf("  W:Word\n");
            printf("  w:Word\n");
            printf("  2:Lemma\n");
            printf("  B:Lemma\n");
            printf("  b:Lemma\n");
            printf("  L:Lemma\n");
            printf("  l:Lemma\n");
            printf("  3:POS tag (create rules for each of them)\n");
            printf("  T:POS tag (create rules for each of them)\n");
            printf("  t:POS tag (create rules for each of them)\n");
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
            printf("  0:Other (don't care). (Use this to skip and ignore e.g. column with POS tags.)\n");
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
            printf("-b: compiled, raw rule file. Exits after pretty printing.\n");
            printf("-Q: Max recursion depth when attempting to create candidate rule\n");
            printf("-q: Percentage of training pairs to set aside for testing\n");
            return Leave;
        case 'D':
            detectFloatingPointNumbers(locoptarg);
            break;
        case 'i': // word list
            i = dupl(locoptarg);
            break;
        case 'j': // temp dir
            j = dupl(locoptarg);
            break;
        case 'n': // columns
            n = dupl(locoptarg);
            break;
        case 'o': // flexrules
            o = dupl(locoptarg);
            break;
        case 'p': // compute parms
            ComputeParms = locoptarg && *locoptarg == '-' ? false : true;
            break;
            // GNU >>
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
        case 'v': // verbose
            Verbose = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'T': // test with training data
            TrainTest = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 't': // test (not with training data)
            Test = locoptarg && *locoptarg == '-' ? false : true;
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
            // << GNU
        case 'W':
            Doweights = locoptarg && *locoptarg == '-' ? false : true;
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
                            if (!isspace(*p))
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
                                        if (!isspace(*q))
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
                            if (!isspace(*p))
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
                                        if (!isspace(*q))
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
            char * tmp = new char[strlen(e) + 10];
            B = tmp;
            sprintf(tmp, "best_%s.txt", e);
            }
        }

    if (o == NULL)
        o = dupl("rules"); //20130125

    if (!i)
        {
        fprintf(stderr, "Error: No training data (-i option)\n");
        exit(2);
        }

    //const char * columns = "12634"; // Word, Lemma, LemmaClass, WordFreq, LemmaFreq
    if (!n)
        n = dupl("FB");// Word, Lemma

    if (!f)
        {
        if (!ComputeParms)
            {
            fprintf(stderr, "Error: No competition function defined (-f option)\n");
            getchar();
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
        }


    FILE * f = fopen(wordList(), "r");
    ++openfiles;
    Blobs = 0; // If there are no non-empty lines, there are no blobs either.
    Lines = 0;
    int kar = 0;
    int prevkar = 0;
    unsigned int lineProfile = 0;
    while ((kar = fgetc(f)) != EOF)
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
            }
        prevkar = kar;
        }
    if(lineProfile & 4)
        { // last line was not finished off with newline
        ++Lines;
        }
    if(Lines > 0)
        ++Blobs;

    if (verbose())
        printf("blobs:%d lines %d\n", Blobs, Lines);
    --openfiles;
    fclose(f);
    }


OptReturnTp optionStruct::readArgs(int argc, char * argv[])
    {
    int optchar;

    OptReturnTp result = GoOn;
    while ((optchar = getopt(argc, argv, opts)) != -1)
        {
        OptReturnTp res = doSwitch(optchar, optarg, argv[0]);
        if (res > result)
            result = res;
        }
    if (!i && !c && !o && !e && !n && !f)
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
            else if (!f)
                f = dupl(argv[optind++]);
            else
                printf("Too many arguments:%s\n", argv[optind]);
            }
        }
    else if (optind < argc)
        {
        if (i && c && o && e && n && f)
            printf("Too many arguments:%s\n", argv[optind]);
        else
            printf("You cannot have a command line with both option-style arguments and option-less-fixed-position arguments:%s\n", argv[optind]);
        }

    completeArgs();
    return result;
    }

void optionStruct::print(FILE * fp) const
    {
    fprintf(fp,"               ; verbose\n-v %s\n", Verbose ? "" : "-");
    if (b)
        {
        fprintf(fp, "               ; raw rules\n-b %s\n", b); if (b) fprintf(fp, "-b %s\n", b); else fprintf(fp, ";-b not specified\n");
        fprintf(fp, "               ; (N/A) raw rules\n;-b not specified\n");
        fprintf(fp, "               ; (N/A) word list\n-i %s\n", i ? i : "?");
        fprintf(fp, "               ; (N/A) extra name suffix\n"); if (e) fprintf(fp, "-e %s\n", e); else fprintf(fp, ";-e not specified\n");
        fprintf(fp, "               ; (N/A) suffix only (%s)\n-s %s\n", SuffixOnly ? "yes" : "no", SuffixOnly ? "" : "-");
        fprintf(fp, "               ; (N/A) columns (1=word,2=lemma,3=tags,0=other)\n-n %s\n", n ? n : "?");
        fprintf(fp, "               ; (N/A) max recursion depth when attempting to create candidate rule\n-Q %d\n", Q);
        fprintf(fp, "               ; (N/A) flex rules\n-o %s\n", o ? o : "?");
        fprintf(fp, "               ; (N/A) temp dir (including separator at end!)\n"); if (j) fprintf(fp, "-j %s\n", j); else fprintf(fp, ";-j not specified\n");
        fprintf(fp, "               ; (N/A) percentage of training pairs to set aside for testing\n-q %d\n", q);
        fprintf(fp, "               ; (N/A) penalties to decide which rule survives\n"); if (nD > 0){ fprintf(fp, "-D "); for (int i = 0; i < nD; ++i)fprintf(fp, "%f;", D[i]); fprintf(fp, "\n"); } else fprintf(fp, ";-D not specified\n");
        fprintf(fp, "               ; (N/A) compute parms (%s)\n-p %s\n", ComputeParms ? "yes" : "no", ComputeParms ? "" : "-");
        fprintf(fp, "               ; (N/A) expected cutoff\n-C %d\n", C);
        fprintf(fp, "               ; (N/A) do weights (%s)\n-W %s\n", Doweights ? "yes" : "no", Doweights ? "" : "-");
        fprintf(fp, "               ; (N/A) current parameters\n-P %s\n", P ? P : "?");
        fprintf(fp, "               ; (N/A) best parameters\n-B %s\n", B ? B : "?");
        fprintf(fp, "               ; (N/A) start training with minimal fraction of training pairs\n-L %f\n", Minfraction);
        fprintf(fp, "               ; (N/A) end training with maximal fraction of training pairs\n-H %f\n", Maxfraction);
        fprintf(fp, "               ; (N/A) number of differently sized fractions of trainingdata\n-K %d\n",K);
        fprintf(fp, "               ; (N/A) number of iterations of training with same fraction of training data when fraction is minimal\n-N %f\n", N);
        fprintf(fp, "               ; (N/A) number of iterations of training with same fraction of training data when fraction is maximal\n-M %f\n", M);
        fprintf(fp, "               ; (N/A) competition function\n-f %s\n", f ? f : "?");
        fprintf(fp, "               ; (N/A) redo training after homographs for next round are removed (%s)\n-R %s\n", Redo ? "yew" : "no", Redo ? "" : "-");
        fprintf(fp, "               ; (N/A) cutoff\n-c %d\n", c);
        }
    else
        {
        fprintf(fp, "               ; raw rules\n;-b not specified\n");
        fprintf(fp, "               ; word list\n-i %s\n", i);
        fprintf(fp, "               ; extra name suffix\n"); if (e) fprintf(fp, "-e %s\n", e); else fprintf(fp, ";-e not specified\n");
        fprintf(fp, "               ; suffix only (%s)\n-s %s\n", SuffixOnly ? "yes" : "no", SuffixOnly ? "" : "-");
        fprintf(fp, "               ; columns (1=word,2=lemma,3=tags,0=other)\n-n %s\n", n);
        fprintf(fp, "               ; max recursion depth when attempting to create candidate rule\n-Q %d\n", Q);
        fprintf(fp, "               ; flex rules\n-o %s\n", o);
        fprintf(fp, "               ; temp dir (including separator at end!)\n"); if (j) fprintf(fp, "-j %s\n", j); else fprintf(fp, ";-j not specified\n");
        fprintf(fp, "               ; percentage of training pairs to set aside for testing\n-q %d\n", q);
        fprintf(fp, "               ; penalties to decide which rule survives\n"); if (nD > 0){ fprintf(fp, "-D "); for (int i = 0; i < nD; ++i)fprintf(fp, "%f;", D[i]); fprintf(fp, "\n"); } else fprintf(fp, ";-D not specified\n");
        fprintf(fp, "               ; compute parms (%s)\n-p %s\n", ComputeParms ? "yes" : "no", ComputeParms ? "" : "-");
        if (ComputeParms)
            {
            assert(P);
            assert(B);
            fprintf(fp, "               ; expected cutoff\n-C %d\n", C);
            fprintf(fp, "               ; do weights (%s)\n-W %s\n", Doweights ? "yes" : "no", Doweights ? "" : "-");
            fprintf(fp, "               ; current parameters\n-P %s\n", P);
            fprintf(fp, "               ; best parameters\n-B %s\n", B);
            fprintf(fp, "               ; start training with minimal fraction of training pairs\n-L %f\n", Minfraction);
            fprintf(fp, "               ; end training with maximal fraction of training pairs\n-H %f\n", Maxfraction);
            fprintf(fp, "               ; number of differently sized fractions of trainingdata\n-K %d\n",K);
            fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is minimal\n-N %f\n", N);
            fprintf(fp, "               ; number of iterations of training with same fraction of training data when fraction is maximal\n-M %f\n", M);
            fprintf(fp, "               ; (N/A) competition function\n-f %s\n", f ? f : "?");
            fprintf(fp, "               ; (N/A) redo training after homographs for next round are removed (%s)\n-R %s\n", Redo ? "yew" : "no", Redo ? "" : "-");
            fprintf(fp, "               ; (N/A) cutoff\n-c %d\n", c);
            }
        else
            {
            assert(f);
            fprintf(fp, "               ; (N/A) expected cutoff\n-C %d\n", C);
            fprintf(fp, "               ; (N/A) do weights (%s)\n-W %s\n", Doweights ? "yes" : "no", Doweights ? "" : "-");
            fprintf(fp, "               ; (N/A) current parameters\n-P %s\n", P ? P : "?");
            fprintf(fp, "               ; (N/A) best parameters\n-B %s\n", B ? B : "?");
            fprintf(fp, "               ; (N/A) start training with minimal fraction of training pairs\n-L %f\n", Minfraction);
            fprintf(fp, "               ; (N/A) end training with maximal fraction of training pairs\n-H %f\n", Maxfraction);
            fprintf(fp, "               ; (N/A) number of differently sized fractions of trainingdata\n-K %d\n",K);
            fprintf(fp, "               ; (N/A) number of iterations of training with same fraction of training data when fraction is minimal\n-N %f\n", N);
            fprintf(fp, "               ; (N/A) number of iterations of training with same fraction of training data when fraction is maximal\n-M %f\n", M);
            fprintf(fp, "               ; competition function\n-f %s\n", f);
            fprintf(fp, "               ; redo training after homographs for next round are removed (%s)\n-R %s\n", Redo ? "yew" : "no", Redo ? "" : "-");
            fprintf(fp, "               ; cutoff\n-c %d\n", c);
            }
        if(this->TrainTest)
            {
            fprintf(fp, "               ; test without training data (%s)\n-T %s\n", TrainTest ? "yes" : "no", TrainTest ? "" : "-");
            }
        if(this->Test)
            {
            fprintf(fp, "               ; test (without training data) (%s)\n-t %s\n", Test ? "yes" : "no", Test ? "" : "-");
            }
        fprintf(fp, "               ; Number of blobs found in word list: %d whereof used for training %d\n", Blobs, FracBlobs == 0 ? Blobs : FracBlobs);
        fprintf(fp, "               ; Number of lines found in word list: %d whereof used for training %d\n", Lines, FracLines == 0 ? Lines : FracLines);
        fprintf(fp, "               ; Current training size step: %d\n", Swath);
        fprintf(fp, "               ; Current iteration in current training size step: %d\n", SwathIteration);
        fprintf(fp, "               ; Current number of nodes: %d\n", NumberOfNodes);
        fprintf(fp, "               ; Current number of lines: %d\n", TrainingPairsLines);
        fprintf(fp, "               ; Nodes/line: %.*e\n", DBL_DIG+2,(double)NumberOfNodes/(double)TrainingPairsLines);
        if(Doweights)
            fprintf(fp, "               ; Current weight: %.*e\n", DBL_DIG+2,Weight);
        else
            fprintf(fp, "               ; Current weight: N/A\n");
        }
    }

void optionStruct::seti(const char * WordList)
    {
    delete i;
    i = dupl(WordList);
    }

void optionStruct::seto(const char * Result)
    {
    delete o;
    o = dupl(Result);
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


void optionStruct::printArgFile() const
    {
    int nameLength = strlen(i) + 1 + (e ? strlen(e) + 1 : 0) + (SuffixOnly ? strlen("suf_") : 0) + (C < 0 ? 0 : 2) + (Doweights ? 3 : 0) + 1;
    char * name = new char[nameLength];
    strcpy(name, i);
    strcat(name, "_");
    if (e)
        {
        strcat(name, e);
        strcat(name, "_");
        }
    if (SuffixOnly)
        {
        strcat(name, "suf_");
        }
    if (C >= 0)
        {
        strcat(name, "C");
        int L = strlen(name);
        name[L] = (char)(C + '0');
        name[L + 1] = 0;
        }
    if (Doweights)
        {
        strcat(name, "_W_");
        }
    int nl = strlen(name) + 1;
    assert(nameLength == nl);
    FILE * fp = fopen(name, "wb");
    ++openfiles;
    print(fp);
    --openfiles;
    fclose(fp);
    }
