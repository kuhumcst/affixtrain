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

#include "applyaffrules.h"
#include "optionaff.h"
#include "argopt.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>


//        printf("usage: makeaffixrules -w <word list> -c <cutoff> -C <expected cutoff> -o <flexrules> -e <extra> -n <columns> -f <compfunc> [<word list> [<cutoff> [<flexrules> [<extra> [<columns> [<compfunc>]]]]]]\n");

bool VERBOSE = false;
static char opts[] = "?@:B:b:c:C:e:f:hH:i:j:L:n:o:O:p:P:s:v:W:t:" /* GNU: */ "wr";
static char *** Ppoptions = NULL;
static char ** Poptions = NULL;
static int optionSets = 0;

char * dupl(const char * s)
    {
    char * d = new char[strlen(s) + 1];
    strcpy(d,s);
    return d;
    }

optionStruct::optionStruct()
    {
    c = NULL; // cutoff 
    C = NULL; // expected cutoff when determining the parameters using weightedcount()
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
	t = NULL; // pretty printed file (see b)
    computeParms = false;// compute parms
    suffixOnly = false;// suffix rules only
    verbose = false;// verbose
    minfraction = -1; // L
    maxfraction = -1; // H
    doweights = false;
    }

optionStruct::~optionStruct()
    {
    for(int I = 0;I < optionSets;++I)
        {
        delete [] Poptions[I];
        delete [] Ppoptions[I];
        }
    delete [] Poptions;
    delete [] Ppoptions;
    delete [] c;
    delete [] C;
    delete [] e;
    delete [] f;
    delete [] i;
    delete [] j;
    delete [] n;
    delete [] o;
    delete [] B;
    delete [] P;
	delete [] b;
	delete [] t;
    }

OptReturnTp optionStruct::doSwitch(int optchar,char * locoptarg,char * progname)
    {
    switch (optchar)
        {
        case '@':
            readOptsFromFile(locoptarg,progname);
            break;
        case 'B': // best parms
            B = dupl(locoptarg);
            break;
        case 'P': // current parms
            P = dupl(locoptarg);
            break;
        case 'c': // cutoff
            c = dupl(locoptarg);
            break;
        case 'C': // cutoff for weightedcount
            C = dupl(locoptarg);
            break;
        case 'e': // extra
            e = dupl(locoptarg);
            break;
        case 'f': // compfunc
            f = dupl(locoptarg);
            break;
        case 'L':
            minfraction = strtod(locoptarg,(char**)0);
            if(minfraction > 1.0)
                {
                printf("%s","Option -L: value must be less than 1.0");
                exit(-1);
                }
            if(minfraction > 0 && maxfraction == -1)
                maxfraction = minfraction;
            break;
        case 'H':
            maxfraction = strtod(locoptarg,(char**)0);
            if(maxfraction > 1.0)
                {
                printf("%s","Option -H: value must be less than 1.0");
                exit(-1);
                }
            if(maxfraction > 0 && minfraction == -1)
                minfraction = maxfraction;
            break;
        case 'h':
        case '?':
            printf("usage:\n"
                "affixtrain [-@ <option file>] -i <word list> [-c <cutoff>] [-C <expected cutoff>] [-o <flexrules>] [-e <extra>] [-n <columns>] [-f <compfunc>] [-p[-]] [-s[-]] [-v[-]] [-j <tempdir>] [-L<n>] [-H<n>]"
                "\nor\n"
                "affixtrain [<word list> [<cutoff> [<flexrules> [<extra> [<columns> [<compfunc>]]]]]]"
                "\nor\n"
				"affixtrain -b <rules> -t <prettyprint>"
                "\n");
            printf("-@: Options are read from file with lines formatted as: -<option letter> <value>\n"
                   "    A semicolon comments out the rest of the line.\n"
                );
            printf("-i: (input) full form / lemma list\n");
            printf("-c: discard rules with little support. 0 <= cutoff <= 9\n");
            printf("-C: (together with -p) expected cutoff (parameter for rule weight function). 0 <= cutoff <= 9\n");
            printf("-o: (output) default is to automatically generate file name\n");
            printf("-e: language code (da,en,is,nl,ru, etc)\n");
            printf("-p: compute parameters (overrules -f)\n");
            printf("-s: create suffix-only rules\n");
            printf("-v: verbose\n");
            printf("-j: directory to store the bulk of intermediate results, also textual presentations of rules. Must not include final (back)slash.\n");
            printf("-L: minimum fraction (with option -f0 or -p)\n");
            printf("-H: maximum fraction (with option -f0 or -p)\n");
            printf("-W: minimise weight, not count (sum of rules) (with -p or -f0)\n");
            printf("-P: write parameters to file (default parms.txt if -p or -f0, otherwise no parameter file)\n");
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
            printf("-b: compiled, raw rule file. Requires -t option. Exits after pretty printing.\n");
            printf("-t: pretty printed output. Requires -b option. Exits after pretty printing.\n");
            return Leave;
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
            computeParms = locoptarg && *locoptarg == '-' ? false : true;
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
            suffixOnly = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'v': // verbose
            verbose = locoptarg && *locoptarg == '-' ? false : true;
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
            doweights = locoptarg && *locoptarg == '-' ? false : true;
            break;
        case 'b': // raw rules
            b = dupl(locoptarg);
            break;
        case 't': // pretty printed rules
            t = dupl(locoptarg);
            break;
        }
    return GoOn;
    }



OptReturnTp optionStruct::readOptsFromFile(char * locoptarg,char * progname)
    {
    char ** poptions;
    char * options;
    FILE * fpopt = fopen(locoptarg,"r");
    OptReturnTp result = GoOn;
    if(fpopt)
        {
        char * p;
        char line[1000];
        int lineno = 0;
        size_t bufsize = 0;
        while(fgets(line,sizeof(line) - 1,fpopt))
            {
            lineno++;
            size_t off = strspn(line," \t");
            if(line[off] == ';')
                continue; // comment line
            if(line[off] == '-')
                {
                off++;
                if(line[off])
                    {
                    char * optarg2 = line + off + 1;
                    size_t off2 = strspn(optarg2," \t");
                    if(!optarg2[off2])
                        optarg2 = NULL;
                    else
                        optarg2 += off2;
                    if(optarg2)
                        {
                        for(p = optarg2 + strlen(optarg2) - 1;p >= optarg2;--p)
                            {
                            if(!isspace(*p))
                                break;
                            *p = '\0';
                            }
                        bool string = false;
                        if(*optarg2 == '\'' || *optarg2 == '"')
                            {
                            for(p = optarg2 + strlen(optarg2) - 1;p > optarg2;--p)
                                {
                                if(*p == *optarg2)
                                    {
                                    string = true;
                                    for(char * q = p + 1;*q;++q)
                                        {
                                        if(*q == ';')
                                            break;
                                        if(!isspace(*q))
                                            {
                                            string = false;
                                            }
                                        }
                                    if(string)
                                        {
                                        *p = '\0';
                                        ++optarg2;
                                        }
                                    break;
                                    }
                                }
                            }
                        if(!*optarg2 && !string)
                            optarg2 = NULL;
                        }
                    if(optarg2)
                        {
                        bufsize += strlen(optarg2) + 1;
                        }
                    char * optpos = strchr(opts,line[off]);
                    if(optpos)
                        {
                        if(optpos[1] != ':')
                            {
                            if(optarg2)
                                {
                                printf("Option argument %s provided for option letter %c that doesn't use it on line %d in option file \"%s\"\n",optarg2,line[off],lineno,locoptarg);
                                exit(1);
                                }
                            }
                        }
                    }
                else
                    {
                    printf("Missing option letter on line %d in option file \"%s\"\n",lineno,locoptarg);
                    exit(1);
                    }
                }
            }
        rewind(fpopt);

        poptions = new char * [lineno];
        options = new char[bufsize];
        // update stacks that keep pointers to the allocated arrays.
        optionSets++;
        char *** tmpPpoptions = new char **[optionSets];
        char ** tmpPoptions = new char *[optionSets];
        int g;
        for(g = 0;g < optionSets - 1;++g)
            {
            tmpPpoptions[g] = Ppoptions[g];
            tmpPoptions[g] = Poptions[g];
            }
        tmpPpoptions[g] = poptions;
        tmpPoptions[g] = options;
        delete [] Ppoptions;
        Ppoptions = tmpPpoptions;
        delete [] Poptions;
        Poptions = tmpPoptions;

        lineno = 0;
        bufsize = 0;
        while(fgets(line,sizeof(line) - 1,fpopt))
            {
            poptions[lineno] = options+bufsize;
            size_t off = strspn(line," \t");
            if(line[off] == ';')
                continue; // comment line
            if(line[off] == '-')
                {
                off++;
                if(line[off])
                    {
                    char * optarg2 = line + off + 1;
                    size_t off2 = strspn(optarg2," \t");
                    if(!optarg2[off2])
                        optarg2 = NULL;
                    else
                        optarg2 += off2;
                    if(optarg2)
                        {
                        for(p = optarg2 + strlen(optarg2) - 1;p >= optarg2;--p)
                            {
                            if(!isspace(*p))
                                break;
                            *p = '\0';
                            }
                        bool string = false;
                        if(*optarg2 == '\'' || *optarg2 == '"')
                            {
                            for(p = optarg2 + strlen(optarg2) - 1;p > optarg2;--p)
                                {
                                if(*p == *optarg2)
                                    {
                                    string = true;
                                    for(char * q = p + 1;*q;++q)
                                        {
                                        if(*q == ';')
                                            break;
                                        if(!isspace(*q))
                                            {
                                            string = false;
                                            }
                                        }
                                    if(string)
                                        {
                                        *p = '\0';
                                        ++optarg2;
                                        }
                                    break;
                                    }
                                }
                            }
                        if(!*optarg2 && !string)
                            optarg2 = NULL;
                        }
                    if(optarg2)
                        {
                        strcpy(poptions[lineno],optarg2);
                        bufsize += strlen(optarg2) + 1;
                        }
                    OptReturnTp res = doSwitch(line[off],poptions[lineno],progname);
                    if(res > result)
                        result = res;
                    }
                }
            lineno++;
            }
        fclose(fpopt);
        }
    else
        {
        printf("Cannot open option file %s\n",locoptarg);
        }
    return result;
    }

OptReturnTp optionStruct::readArgs(int argc, char * argv[])
    {
    int optchar;
    OptReturnTp result = GoOn;
    while((optchar = getopt(argc,argv, opts)) != -1)
        {
        OptReturnTp res = doSwitch(optchar,optarg,argv[0]);
        if(res > result)
            result = res;
        }
    if(!i && !c && !o && !e && !n && !f)
        {
        while(optind < argc)
            {
            if(!i)
                i = argv[optind++];
            else if(!c)
                c = argv[optind++];
            else if(!o)
                o = argv[optind++];
            else if(!e)
                e = argv[optind++];
            else if(!n)
                n = argv[optind++];
            else if(!f)
                f = argv[optind++];
            else
                printf("Too many arguments:%s\n",argv[optind]);
            }
        }
    else if(optind < argc)
        {
        if(i && c && o && e && n && f)
            printf("Too many arguments:%s\n",argv[optind]);
        else
            printf("You cannot have a command line with both option-style arguments and option-less-fixed-position arguments:%s\n",argv[optind]);
        }
    return result;
    }

