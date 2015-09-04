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

#include "comp.h"
#include "affixtrain.h"
#include "graph.h"
#include "optionaff.h"
#include <float.h>

#define ZIGGURAT 1

#if ZIGGURAT
#include "rnorrexp.c"
#endif


#define NPARMS parms.ROWPARMS

/*
ACL 2009 paper:
Icelandic   71.3    1.5 even_better (71,30 1,51 iflg. D:\dokumenter\tvärsök\even_better\icelandic.xls) peen 71,51 1,65 sugar 70,93 1,86 affiksFEW3 71,02 2,16 no pruning
Danish      92.8    0.2 peen sugar: 92,72 0,19 no pruning
Norwegian   87.6    0.3 affiksFEW2 sugar: 86,67 0,68
Greek       90.4    0.4 sugar no pruning
Slovene     86.7    0.3 affiksFEW3 affiksFEW2: 86,23 0,58 sugar: 86,27 0,41 peen:86,13 0,55  0,4
Swedish     92.3    0.3 sugar pruning 1
German      91.46   0.17 sugar no pruning
English     89.0    1.3 sugar pruning 2
Dutch       90.4    0.5 affiksFEW2 sugar: 90,17 0,31    0,3 no pruning
Polish      93.88   0.08 peen sugar: 93,88 0,08 (?) no pruning
*/

#if _NA
// IMPORTANT (20090511) R__NA and W__NA are not updated as sibling rules are 
// added and eat up the training pairs that earlier siblings did not handle.
// This error was detected after having used the weight functions for 
// the ACL-paper.

static int comp_fairly_good(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //fairly good, Icelandic 71.270883
//AMBI:
    // French ok 85.767516 ambi1 1.156051 ambi2 0.955414 diff 12.121019 rules 7337.500000 2.731849% cutoff 2
    int A1 = a->W__R           + a->R__R;
    int B1 = b->W__R           + b->R__R;
    int A2 = a->W__R + a->W__W           + a->R__NA;
    int B2 = b->W__R + b->W__W           + b->R__NA;
    int A3 = a->W__R           + a->R__R + a->R__NA;
    int B3 = b->W__R           + b->R__R + b->R__NA;
/*  int A2 = a->R__NA - a->W__NA;
    int B2 = b->R__NA - b->W__NA;
    int A3 = a->W__R - a->R__W;
    int B3 = b->W__R - b->R__W;
*/
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_even_better(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //even better, Icelandic 71.300716
    // BEST Icelandic 71.535870 +/- 1.919590 at 0.9856 of dataset, 17 iterations, 23209.882353 = 40.477646% rules, cutoff = 0
    // Icelandic 71.283167 +/- 1.714260 at 0.9856 of dataset, 17 iterations, 22719.470588 = 39.622376% rules, cutoff = 0, RECURSE == 4
//AMBI:
    // French ok 85.487261 ambi1 1.283439 ambi2 1.050955 diff 12.178344 rules 7360.125000 2.740283% cutoff 2

    int A1 = a->W__R           + a->R__R;// wr + rr
    int B1 = b->W__R           + b->R__R;
    int A2 = a->W__R           + a->R__R + a->R__NA;// wr + rr + rn - r = wr - rw
    int B2 = b->W__R           + b->R__R + b->R__NA;
    int A3 = a->W__R + a->W__W           + a->R__NA;// wr + ww + rn - w = -wn + rn
    int B3 = b->W__R + b->W__W           + b->R__NA;
//    int A2 = a->W__R - a->R__W;
//    int B2 = b->W__R - b->R__W;
//    int A3 = a->R__NA - a->W__NA;
//    int B3 = b->R__NA - b->W__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_affiksFEW3(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    // Icelandic 65.781623, cutoff 1 (old lemmatizer 73.329356, cutoff 0)
    // Icelandic 66.544995 +/- 1.943469 at 0.9856 of dataset, 17 iterations, 11134.176471 = 19.417817% rules, cutoff = 1
    // English 87.863636, cutoff 2 (old 87.954545, cutoff 1)
    // English 87.806061 +/- 1.009323 at 0.9856 of dataset, 15 iterations, 1619.133333 = 2.152101% rules (cutoff = 2)
    // BEST Slovene 86.669776 +/- 0.331106 at 0.9856 of dataset, 9 iterations, 5650.777778 = 2.888237% rules (cutoff = 2)
    // Slovene-ambi (4.23%) 83.165661, cutoff 3 (3550 rules!) (old 82.017103, cutoff 1, 9377 rules) Better than _affiksFEW2, 82.780013, 6656 rules.
    // Danish 90.942165    +/- 0.589437 at 0.9856 of dataset, 5 iterations, 32327.400000 = 5.925881% rules, cutoff = 1
    // German 90.266461 +/-    0.509202 at 0.9856 of dataset, 7 iterations, 21539.428571 = 6.930653% rules, cutoff = 1
    // Greek 89.640779 +/- 0.402079 at 0.9856 of dataset, 5 iterations, 13377.200000 = 2.472132% rules, cutoff = 2
    // Dutch 87.817059 +/- 0.366236 at 0.9856 of dataset, 7 iterations, 23493.571429 = 7.895486% rules, cutoff = 1
    // Norwegian 85.788507 +/- 0.484921 at 0.9856 of dataset, 6 iterations, 14904.000000 = 3.157580% rules, cutoff = 2
    // Polish 93.203365 +/- 0.175436 at 0.9856 of dataset, 2 iterations, 50597.500000 = 1.491153% rules, cutoff = 2
    // Swedish 91.709042 +/- 0.170094 at 0.9856 of dataset, 6 iterations, 4407.666667 = 0.935737% rules, cutoff = 3
//AMBI:
    // French ok 82.754777 ambi1 2.353503 ambi2 1.805732 diff 13.085987 rules 7360.125000 2.740283% cutoff 2

    /* Interesting because it generates far less rules than the above 
       variables, only 20 % more than the old lemmatizer.
       Also interesting is that there are not many leaves with only one
       supporting training pair.
       Yet, the leaves with only one supporter are detrimentous to the overall
       result (cutoff has to be 1 or even 2).
    */
#if 1
    int A1 = a->W__R           + a->R__R + a->R__NA; // Good: previously wrong words got it right. Bad: previously right words got it wrong.
    int B1 = b->W__R           + b->R__R + b->R__NA;
    int A2 = a->W__R           + a->R__R; // Good: any rightly lemmatized word
    int B2 = b->W__R           + b->R__R;
    int A3 = a->W__R + a->W__W           + a->R__NA; // Good: previously right words that didn't match. They may return to the parent.
    int B3 = b->W__R + b->W__W           + b->R__NA; // Bad: previously wrong words that didn't match. They must be handled by siblings.
#else
    int A1 = a->W__R - a->R__W; // Good: previously wrong words got it right. Bad: previously right words got it wrong.
    int B1 = b->W__R - b->R__W;
    int A2 = a->W__R + a->R__R; // Good: any rightly lemmatized word
    int B2 = b->W__R + b->R__R;
    int A3 = a->W__R + a->W__W - a->R__R - a->R__W; // Good: previously right words that didn't match. They may return to the parent.
    int B3 = b->W__R + b->W__W - b->R__R - b->R__W; // Bad: previously wrong words that didn't match. They must be handled by siblings.
#endif
/*
    int A1 = a->W__R - a->R__W; // Good: previously wrong words got it right. Bad: previously right words got it wrong.
    int B1 = b->W__R - b->R__W;
    int A2 = a->W__R + a->R__R; // Good: any rightly lemmatized word
    int B2 = b->W__R + b->R__R;
    int A3 = a->R__NA - a->W__NA; // Good: previously right words that didn't match. They may return to the parent.
    int B3 = b->R__NA - b->W__NA; // Bad: previously wrong words that didn't match. They must be handled by siblings.
*/
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_affiksFEW(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
/**/
    //_affiksFEW
    // Dutch 88.138224, 39943.5 flexrules cutoff 1 (old 89.656164, 47277.75 flexrules, cutoff 1)
    // German 90.266461    +/-    0.509202 at 0.9856 of dataset, 7 iterations, 21539.428571 = 6.930653% rules, cutoff = 1
//AMBI:
    // French ok 82.617834 ambi1 2.455414 ambi2 1.872611 diff 13.054140 rules 7360.125000 2.740283% cutoff 2
    int N = a->W__W + a->W__R + a->W__NA + a->R__W + a->R__R + a->R__NA;

    int A1;
    int B1;
    int A2;
    int B2;
    int A3;
    int B3;

    // good for small numbers:
    if(N < 3)
        {
        A1 = a->W__R           + a->R__R;
        B1 = b->W__R           + b->R__R;
        A2 = a->W__R           + a->R__R + a->R__NA;
        B2 = b->W__R           + b->R__R + b->R__NA;
        A3 = a->W__R + a->W__W           + a->R__NA;
        B3 = b->W__R + b->W__W           + b->R__NA;
/*      A1 = a->W__R + a->R__R;
        B1 = b->W__R + b->R__R;
        A2 = a->W__R - a->R__W;
        B2 = b->W__R - b->R__W;
        A3 = a->R__NA - a->W__NA;
        B3 = b->R__NA - b->W__NA;
*/
        }
    // good for big numbers:
    else
        {
        A1 = a->W__R + a->R__R + a->R__NA; // Good: previously wrong words got it right. Bad: previously right words got it wrong.
        B1 = b->W__R + b->R__R + b->R__NA;
        A2 = a->W__R + a->R__R; // Good: any rightly lemmatized word
        B2 = b->W__R + b->R__R;
        A3 = a->R__NA + a->W__R + a->W__W; // Good: previously right words that didn't match. They may return to the parent.
        B3 = b->R__NA + b->W__R + a->W__W; // Bad: previously wrong words that didn't match. They must be handled by siblings.
/*      A1 = a->W__R - a->R__W; // Good: previously wrong words got it right. Bad: previously right words got it wrong.
        B1 = b->W__R - b->R__W;
        A2 = a->W__R + a->R__R; // Good: any rightly lemmatized word
        B2 = b->W__R + b->R__R;
        A3 = a->R__NA - a->W__NA; // Good: previously right words that didn't match. They may return to the parent.
        B3 = b->R__NA - b->W__NA; // Bad: previously wrong words that didn't match. They must be handled by siblings.
*/
        }
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_affiksFEW2(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //_affiksFEW2
    // (OK) BEST Dutch 90.452096 +/- 0.655431 at 0.9856 of dataset, 7 iterations, 53607.714286 = 18.015948% rules, cutoff = 0
    // (OK) Norwegian 86.776860 +/- 0.642621 at 0.9856 of dataset, 6 iterations, 112374.000000 = 23.807698% rules, cutoff = 0
    // (OK) English 88.424242 +/- 1.191106 at 0.9856 of dataset, 15 iterations, 1383.000000 = 1.838240% rules, cutoff = 2
    // (OK) Icelandic 71.304226 +/- 1.453643 at 0.9856 of dataset, 17 iterations, 25635.000000 = 44.707011% rules, cutoff = 0
    // (OK) German 91.156762 +/- 0.348391 at 0.9856 of dataset, 7 iterations, 48816.571429 = 15.707506% rules, cutoff = 0
    // (OK) Slovene 86.537639 +/- 0.559484 at 0.9856 of dataset, 9 iterations, 40643.444444 = 20.773759% rules, cutoff = 0
    // (OK) Swedish 91.907598 +/- 0.224888 at 0.9856 of dataset, 6 iterations, 27958.000000 = 5.935415% rules, cutoff = 1
    // (OK) Greek 90.741209 +/- 0.312526 at 0.9856 of dataset, 5 iterations, 125306.400000 = 23.156860% rules, cutoff = 0
    // (OK) Danish 92.994605 +/- 0.210674 at 0.9856 of dataset, 5 iterations, 67278.800000 = 12.332763% rules, cutoff = 0
    
    // (?)ALMOST BEST Polish 93.398015 +/- 0.045642 at 0.9856 of dataset, 2 iterations, 165511.500000 = 4.877770% rules, cutoff = 1
//AMBI:
    // French ok 84.194268 ambi1 2.277070 ambi2 1.576433 diff 11.952229 rules 6453.250000 2.402640% cutoff 2
/* 20140922, 

nohup nice /home/zgk261/bin/testrules -I -D /home/zgk261/sandkasse/nl/dict_nl_non_ambiguous -L nl -C affiksFEW2 -A >/home/zgk261/sandkasse/nl/testout 2>/home/zgk261/sandkasse/nl/testerr &

cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          7  297719.000000   41163.857143   50127.285714      13.826413      16.837113    3830.428571       0.000000       0.000000       0.000000     519.571429    3922.428571       0.000000       0.000000       0.000000     427.571429       0.539105       0.000000       0.000000       0.000000       0.539105       0.440368       0.000000       0.000000       0.000000       0.440368      88.055829       0.000000       0.000000       0.000000      11.944171       0.000000      90.170772       0.000000       0.000000       0.000000       9.829228       0.000000       0.000000       0.000000       0.000000     100.000000       0.000000           -nan       0.000000       0.000000       0.000000     100.000000       0.000000           -nan
1       0.985600          7  297719.000000   11178.571429   22834.000000       3.754739       7.669648    3803.285714       0.000000       0.000000       0.000000     546.714286    3914.285714       0.000000       0.000000       0.000000     435.714286       0.587252       0.000000       0.000000       0.000000       0.587252       0.570207       0.000000       0.000000       0.000000       0.570207      87.431856       0.000000       0.000000       0.000000      12.568144       0.000000      89.983580       0.000000       0.000000       0.000000      10.016420       0.000000       0.000000       0.000000       0.000000     100.000000       0.000000           -nan       0.000000       0.000000       0.000000     100.000000       0.000000           -nan
2       0.985600          7  297719.000000    6855.000000    8014.285714       2.302507       2.691896    3757.857143       0.000000       0.000000       0.000000     592.142857    3895.857143       0.000000       0.000000       0.000000     454.142857       0.549072       0.000000       0.000000       0.000000       0.549072       0.530304       0.000000       0.000000       0.000000       0.530304      86.387521       0.000000       0.000000       0.000000      13.612479       0.000000      89.559934       0.000000       0.000000       0.000000      10.440066       0.000000       0.000000       0.000000       0.000000     100.000000       0.000000           -nan       0.000000       0.000000       0.000000     100.000000       0.000000           -nan
3       0.985600          7  297719.000000    5109.142857    5404.142857       1.716096       1.815182    3726.285714       0.000000       0.000000       0.000000     623.714286    3868.571429       0.000000       0.000000       0.000000     481.428571       0.586802       0.000000       0.000000       0.000000       0.586802       0.407759       0.000000       0.000000       0.000000       0.407759      85.661741       0.000000       0.000000       0.000000      14.338259       0.000000      88.932677       0.000000       0.000000       0.000000      11.067323       0.000000       0.000000       0.000000       0.000000     100.000000       0.000000           -nan       0.000000       0.000000       0.000000     100.000000       0.000000           -nan
4       0.985600          7  297719.000000    4082.571429    4193.857143       1.371283       1.408663    3694.285714       0.000000       0.000000       0.000000     655.714286    3844.714286       0.000000       0.000000       0.000000     505.285714       0.695342       0.000000       0.000000       0.000000       0.695342       0.391416       0.000000       0.000000       0.000000       0.391416      84.926108       0.000000       0.000000       0.000000      15.073892       0.000000      88.384236       0.000000       0.000000       0.000000      11.615764       0.000000       0.000000       0.000000       0.000000     100.000000       0.000000           -nan       0.000000       0.000000       0.000000     100.000000       0.000000           -nan
5       0.985600          7  297719.000000    3462.714286    3482.571429       1.163081       1.169751    3667.000000       0.000000       0.000000       0.000000     683.000000    3827.000000       0.000000       0.000000       0.000000     523.000000       0.788006       0.000000       0.000000       0.000000       0.788006       0.407141       0.000000       0.000000       0.000000       0.407141      84.298851       0.000000       0.000000       0.000000      15.701149       0.000000      87.977011       0.000000       0.000000       0.000000      12.022989       0.000000       0.000000       0.000000       0.000000     100.000000       0.000000           -nan       0.000000       0.000000       0.000000     100.000000       0.000000           -nan
cutoff 0 Affix  a       0.676984 b       0.804745: N(rules)=      1.967934*N(trainpairs)^0.804745
         Suffix a       0.181309 b       0.828157: N(rules)=      1.198786*N(trainpairs)^0.828157
cutoff 1 Affix  a       0.020304 b       0.794314: N(rules)=      1.020512*N(trainpairs)^0.794314
         Suffix a      -1.577537 b       0.865510: N(rules)=      0.206483*N(trainpairs)^0.865510
cutoff 2 Affix  a      -1.037201 b       0.799804: N(rules)=      0.354446*N(trainpairs)^0.799804
         Suffix a      -2.496919 b       0.900882: N(rules)=      0.082338*N(trainpairs)^0.900882
cutoff 3 Affix  a      -1.918028 b       0.840841: N(rules)=      0.146896*N(trainpairs)^0.840841
         Suffix a      -2.960479 b       0.913754: N(rules)=      0.051794*N(trainpairs)^0.913754
cutoff 4 Affix  a      -2.344658 b       0.853928: N(rules)=      0.095880*N(trainpairs)^0.853928
         Suffix a      -3.219691 b       0.916661: N(rules)=      0.039967*N(trainpairs)^0.916661
cutoff 5 Affix  a      -2.716219 b       0.868745: N(rules)=      0.066124*N(trainpairs)^0.868745
         Suffix a      -3.525562 b       0.927406: N(rules)=      0.029435*N(trainpairs)^0.927406

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
cutoff                      0
fraction          9856.000000
iterations                  7
trainlines      297719.000000
rules            50127.285714 (  41163.857143)
rules%              16.837113 (     13.826413)
same%stdev           0.440368
ambi1%stdev          0.000000
ambi2%stdev          0.000000
ambi3%stdev          0.000000
diff%stdev           0.440368
same%               90.170772 (     88.055829)
ambi1%               0.000000 (      0.000000)
ambi2%               0.000000 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%                9.829228 (     11.944171)
amb.rules%           0.000000 (      0.000000)
false_amb%           0.000000 (      0.000000)
false_not_amb%       0.000000 (      0.000000)
true_amb%            0.000000 (      0.000000)
true_not_amb%      100.000000 (    100.000000)
precision            0.000000 (      0.000000)
recall                   -nan (          -nan)
*/
/* 20140922, same as above, except that data are ambiguous.

nohup nice /home/zgk261/bin/testrules -I -D /home/zgk261/sandkasse/nl/dict_nl_without_doubles_UTF8 -L nl -C affiksFEW2 -A >/home/zgk261/sandkasse/nl/testout 2>/home/zgk261/sandkasse/nl/testerr &

cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          7  306848.000000   73132.000000   50311.000000      23.833299      16.396066    3635.857143      78.000000     117.000000       4.714286     648.428571    3843.428571       0.000000       0.000000       0.000000     640.571429       0.368578       0.192276       0.216604       0.027956       0.556177       0.484660       0.000000       0.000000       0.000000       0.484660      81.085128       1.739518       2.609277       0.105136      14.460940       4.982796      85.714286       0.000000       0.000000       0.000000      14.285714       0.000000       4.552695       6.203008       0.430101      88.814197       0.045105       0.064841       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000
1       0.985600          7  306848.000000   18391.714286   22871.857143       5.993754       7.453807    3743.714286       0.000000       0.000000       0.000000     740.285714    3897.142857       0.000000       0.000000       0.000000     586.857143       0.341564       0.000000       0.000000       0.000000       0.341564       0.373271       0.000000       0.000000       0.000000       0.373271      83.490506       0.000000       0.000000       0.000000      16.509494       0.000000      86.912196       0.000000       0.000000       0.000000      13.087804       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000
2       0.985600          7  306848.000000   11355.142857    8029.142857       3.700576       2.616652    3696.857143       0.000000       0.000000       0.000000     787.142857    3909.714286       0.000000       0.000000       0.000000     574.285714       0.396739       0.000000       0.000000       0.000000       0.396739       0.301886       0.000000       0.000000       0.000000       0.301886      82.445521       0.000000       0.000000       0.000000      17.554479       0.000000      87.192558       0.000000       0.000000       0.000000      12.807442       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000
3       0.985600          7  306848.000000    8112.285714    5438.857143       2.643747       1.772492    3652.571429       0.000000       0.000000       0.000000     831.428571    3888.428571       0.000000       0.000000       0.000000     595.571429       0.255206       0.000000       0.000000       0.000000       0.255206       0.205149       0.000000       0.000000       0.000000       0.205149      81.457882       0.000000       0.000000       0.000000      18.542118       0.000000      86.717854       0.000000       0.000000       0.000000      13.282146       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000
4       0.985600          7  306848.000000    6239.285714    4213.000000       2.033347       1.372992    3615.000000       0.000000       0.000000       0.000000     869.000000    3864.000000       0.000000       0.000000       0.000000     620.000000       0.273440       0.000000       0.000000       0.000000       0.273440       0.204397       0.000000       0.000000       0.000000       0.204397      80.619982       0.000000       0.000000       0.000000      19.380018       0.000000      86.173060       0.000000       0.000000       0.000000      13.826940       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000
5       0.985600          7  306848.000000    5044.285714    3494.285714       1.643904       1.138768    3592.714286       0.000000       0.000000       0.000000     891.285714    3847.428571       0.000000       0.000000       0.000000     636.571429       0.262209       0.000000       0.000000       0.000000       0.262209       0.203118       0.000000       0.000000       0.000000       0.203118      80.122977       0.000000       0.000000       0.000000      19.877023       0.000000      85.803492       0.000000       0.000000       0.000000      14.196508       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000       0.000000       6.633108       0.000000      93.366892       0.000000       0.000000
cutoff 0 Affix  a       0.884193 b       0.790803: N(rules)=      2.421031*N(trainpairs)^0.790803
         Suffix a      -0.706998 b       0.933881: N(rules)=      0.493123*N(trainpairs)^0.933881
cutoff 1 Affix  a       0.242492 b       0.779144: N(rules)=      1.274421*N(trainpairs)^0.779144
         Suffix a      -2.482111 b       0.967669: N(rules)=      0.083567*N(trainpairs)^0.967669
cutoff 2 Affix  a      -0.969721 b       0.794266: N(rules)=      0.379189*N(trainpairs)^0.794266
         Suffix a      -3.356905 b       0.996574: N(rules)=      0.034843*N(trainpairs)^0.996574
cutoff 3 Affix  a      -1.776403 b       0.828451: N(rules)=      0.169246*N(trainpairs)^0.828451
         Suffix a      -3.794027 b       1.003571: N(rules)=      0.022505*N(trainpairs)^1.003571
cutoff 4 Affix  a      -2.351608 b       0.854877: N(rules)=      0.095216*N(trainpairs)^0.854877
         Suffix a      -4.020226 b       1.000475: N(rules)=      0.017949*N(trainpairs)^1.000475
cutoff 5 Affix  a      -2.757390 b       0.871710: N(rules)=      0.063457*N(trainpairs)^0.871710
         Suffix a      -4.205819 b       0.997875: N(rules)=      0.014909*N(trainpairs)^0.997875

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
cutoff                      2
fraction          9856.000000
iterations                  7
trainlines      306848.000000
rules             8029.142857 (  11355.142857)
rules%               2.616652 (      3.700576)
same%stdev           0.301886
ambi1%stdev          0.000000
ambi2%stdev          0.000000
ambi3%stdev          0.000000
diff%stdev           0.301886
same%               87.192558 (     82.445521)
ambi1%               0.000000 (      0.000000)
ambi2%               0.000000 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               12.807442 (     17.554479)
amb.rules%           0.000000 (      0.000000)
false_amb%           0.000000 (      0.000000)
false_not_amb%       6.633108 (      6.633108)
true_amb%            0.000000 (      0.000000)
true_not_amb%       93.366892 (     93.366892)
precision            0.000000 (      0.000000)
recall               0.000000 (      0.000000)
*/

#if 1 // 20090511
    int A1 = a->W__R           + 2*a->R__R + a->R__NA; // good: all words that are lemmatised correctly. bad: all previously right words that got it wrong
    // wr + 2rr + rn - r = ww + rr - rw
    int B1 = b->W__R           + 2*b->R__R + b->R__NA;
    int A2 = a->W__R           +   a->R__R + a->R__NA;
    // wr + rr + rn - r = wr - rw
    int B2 = b->W__R           +   b->R__R + b->R__NA;
    int A3 = a->W__R + a->W__W             + a->R__NA;
    // wr + ww + rn - w = -wn + rn
    int B3 = b->W__R + b->W__W             + b->R__NA;
#else
    int A1 = a->W__R           + a->R__R - a->R__W; // good: all words that are lemmatised correctly. bad: all previously right words that got it wrong
    // wr + 2rr + rn - r = ww + rr - rw
    int B1 = b->W__R           + b->R__R - b->R__W;
    int A2 = a->W__R                     - a->R__W;
    // wr + rr + rn - r = wr - rw
    int B2 = b->W__R                     - b->R__W;
    int A3 = a->W__R + a->W__W - a->R__R - a->R__W;
    // wr + ww + rn - w = -wn + rn
    int B3 = b->W__R + b->W__W - b->R__R - b->R__W;
#endif
/*  int A1 = a->W__R + a->R__R - a->R__W; // good: all words that are lemmatised correctly. bad: all previously right words that got it wrong
    int B1 = b->W__R + b->R__R - b->R__W;
    int A2 = a->W__R - a->R__W;
    int B2 = b->W__R - b->R__W;
    int A3 = a->R__NA - a->W__NA;
    int B3 = b->R__NA - b->W__NA;
*/
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_affiksFEW2org(const vertex * a,const vertex * b)
    {
    // BEST Norwegian 87.494563 +/- 0.217147 at 0.9856 of dataset, 6 iterations, 101814.500000 = 21.570549% rules, cutoff = 0
    // English 88.260606 +/- 0.826699 at 0.9856 of dataset, 15 iterations, 7362.466667 = 9.785960% rules, cutoff = 1
    // Icelandic 70.651411 +/- 1.565857 at 0.9856 of dataset, 17 iterations, 23232.941176 = 40.517860% rules, cutoff = 0
    // German 90.307358 +/- 0.355867 at 0.9856 of dataset, 7 iterations, 50595.857143 = 16.280019% rules, cutoff = 0
    // Dutch 90.274675 +/- 0.462929 at 0.9856 of dataset, 7 iterations, 23452.142857 = 7.881563% rules, cutoff = 1    
    // Slovene 86.417162 +/- 0.540735 at 0.9856 of dataset, 9 iterations, 40847.666667 = 20.878142% rules, cutoff = 0
    // Swedish 91.982663 +/- 0.250703 at 0.9856 of dataset, 6 iterations, 28998.000000 = 6.156204% rules, cutoff = 1
    // Greek 90.258032 +/- 0.234665 at 0.9856 of dataset, 5 iterations, 43156.000000 = 7.975310% rules, cutoff = 1 (but exactly the same as cutoff = 0)
    // Danish 92.425041 +/- 0.374415 at 0.9856 of dataset, 5 iterations, 73177.800000 = 13.414099% rules, cutoff = 0
//AMBI:
    // French ok 84.761146 ambi1 2.015924 ambi2 1.665605 diff 11.557325 rules 7262.500000 2.703935% cutoff 2
    int A1 = a->W__R + a->R__R - a->R__W; // good: all words that are lemmatised correctly. bad: all previously right words that got it wrong
    int B1 = b->W__R + b->R__R - b->R__W;
    int A2 = a->W__R - a->R__W;
    int B2 = b->W__R - b->R__W;
    int A3 = a->R__NA - a->W__NA;
    int B3 = b->R__NA - b->W__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_fixNA(const vertex * a,const vertex * b)
    {
    /*
    Icelandic 47.982267 (at 0.8488 of dataset)
    */
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //_fixNA
//AMBI:
    // French: stopped because of very bad results. (> 25% wrong results) 

    int A1 = a->W__R + a->W__NA           + a->R__NA;
    int B1 = b->W__R + b->W__NA           + b->R__NA;
    int A2 = a->W__R            + a->R__R + a->R__NA;
    int B2 = b->W__R            + b->R__R + b->R__NA;
    int A3 = a->W__R                      + a->R__NA;
    int B3 = b->W__R                      + b->R__NA;

    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_fruit(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //Icelandic 71.344041 at 0.939 of dataset
    //ALMOST BEST Icelandic 71.521831 +/- 1.988737 at 0.9856 of dataset, 17 iterations, 23539.352941 = 41.052237% rules
    //Slovene 85.900276 +/- 0.456532 at 0.9856 of dataset, 9 iterations, 42167.333333 = 21.552652% rules
    //English 87.626771    +/- 0.060148 at 0.4928 (!) of dataset, 3 iterations, 933.000000 = 2.480262% rules
//AMBI:
    // French ok 85.382166 ambi1 1.359873 ambi2 1.089172 diff 12.168790 rules 7259.125000 2.899075% cutoff 2
    int A1 = a->W__R            + a->R__R;
    int B1 = b->W__R            + b->R__R;
    int A2 = a->W__R            + a->R__R + a->R__NA;
    int B2 = b->W__R            + b->R__R + b->R__NA;
    int A3 = a->W__R + a->W__NA           + a->R__NA;
    int B3 = b->W__R + b->W__NA           + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_ice(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //Icelandic 60.242322 at 0.939 of dataset
//AMBI:
    // French ok 82.557325 ambi1 2.522293 ambi2 1.866242 diff 13.054140 rules 8556.625000 3.185757% cutoff 2
    int A1 = a->W__R            + a->R__R + a->R__NA;
    int B1 = b->W__R            + b->R__R + b->R__NA;
    int A2 = a->W__R            + a->R__R;
    int B2 = b->W__R            + b->R__R;
    int A3 = a->W__R + a->W__NA           + a->R__NA;
    int B3 = b->W__R + b->W__NA           + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_pisang(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //Icelandic 71.287687 at 0.939 of dataset
//AMBI:
    // French ok 85.414013 ambi1 1.359873 ambi2 1.085987 diff 12.140127 rules 7848.375000 2.922065% cutoff 2
    int A1 = a->W__R            + a->R__R;
    int B1 = b->W__R            + b->R__R;
    int A2 = a->W__R + a->W__NA           + a->R__NA;
    int B2 = b->W__R + b->W__NA           + b->R__NA;
    int A3 = a->W__R            + a->R__R + a->R__NA;
    int B3 = b->W__R            + b->R__R + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_kiwi(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //Icelandic 70.865032 at 0.939 of dataset
//AMBI:
    // French ok 85.410828 ambi1 1.378981 ambi2 1.035032 diff 12.175159 rules 7676.875000 2.858213% cutoff 2
    int A1 = a->W__R            + a->R__R;
    int B1 = b->W__R            + b->R__R;
    int A2 = a->W__R + a->W__NA + a->R__R + a->R__NA;
    int B2 = b->W__R + b->W__NA + b->R__R + b->R__NA;
    int A3 = a->W__R                      + a->R__NA;
    int B3 = b->W__R                      + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_carrot(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //Icelandic 71.090448 at 0.939 of dataset
//AMBI:
    // French ok 85.060510 ambi1 1.328025 ambi2 1.041401 diff 12.570064 rules 7241.625000 2.696163% cutoff 2
    int A1 = 4*(a->W__R        + a->R__R) + a->R__NA;
    int B1 = 4*(b->W__R        + b->R__R) + b->R__NA;
    int A2 = a->W__R            + a->R__R + a->R__NA;
    int B2 = b->W__R            + b->R__R + b->R__NA;
    int A3 = a->W__R + a->W__NA           + a->R__NA;
    int B3 = b->W__R + b->W__NA           + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_peen(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    // Icelandic 71.344041 at 0.939 of dataset
    // ALMOST BEST Icelandic 71.507792 +/- 1.645702 at 0.9856 of dataset, 17 iterations, 25240.882353 = 44.019676% rules
    // Slovene 86.133458    +/- 0.549185 at 0.9856 of dataset, 9 iterations, 40898.777778 = 20.904266% rules
    // English 87.803261    +/- 0.106156 at 0.4928 (!) of dataset, 3 iterations, 889.333333 = 2.364179% rules
    // Dutch 89.837692 +/- 0.412795 at 0.9856 of dataset, 7 iterations, 56640.285714 = 19.035104% rules, cutoff = 0
    // ALMOST BEST German 91.288892 +/- 0.670828 at 0.9856 of dataset, 7 iterations, 50584.857143 = 16.276480% rules, cutoff = 0
    // Swedish 91.873698 +/- 0.367967 at 0.9856 of dataset, 6 iterations, 9066.166667 = 1.924725% rules, cutoff = 2
    // ALMOST BEST Norwegian 87.535644 +/- 0.344659 at 0.9856 of dataset, 6 iterations, 48468 = 10.268492% rules, cutoff = 1
    // ALMOST BEST Greek 90.414875+/- 0.385254 at 0.9856 of dataset, 5 iterations, 120691.4 = 22.303999% rules, cutoff = 0
    // BEST Danish 92.796387 +/- 0.214267 at 0.9856 of dataset, 5 iterations, 67807 = 12.429587% rules, cutoff = 0
    // ALMOST BEST Russian 80.484806 +/- 0.409391 at 0.9856 of dataset, 6 iterations, 54630 = 14.022614% rules, cutoff = 1
    // BEST Polish 93.880103 +/- 0.077021 at 0.9856 of dataset, 2 iterations, 344944.5    = 10.165818% rules, cutoff = 0
//AMBI:
    // French ok 84.993631 ambi1 1.388535 ambi2 1.085987 diff 12.531847 rules 7318.375000 2.724738% cutoff 2
/*
0    0.985600    2    3490123.000000    0.000000    415069.500000    0.000000    11.892690    0.000000    0.000000    0.000000    0.000000    0.000000    45811.000000    299.500000    168.000000    0.000000    4714.500000    0.000000    0.000000    0.000000    0.000000    0.000000    0.044374    0.026347    0.033280    0.000000    0.015253    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    89.837821    0.587336    0.329457    0.000000    9.245387    1.282529    0.000000    0.000000    0.000000    0.000000    0.000000    nan    1.013865    6.709744    0.268664    92.007727    0.116994    0.038499
1    0.985600    2    3490123.000000    0.000000    198203.500000    0.000000    5.678983    0.000000    0.000000    0.000000    0.000000    0.000000    46176.500000    283.000000    241.000000    0.000000    4292.500000    0.000000    0.000000    0.000000    0.000000    0.000000    0.092907    0.063787    0.016640    0.000000    0.012480    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    90.554586    0.554978    0.472614    0.000000    8.417822    1.257035    0.000000    0.000000    0.000000    0.000000    0.000000    nan    0.784421    6.505795    0.472614    92.237170    0.231508    0.067725
2    0.985600    2    3490123.000000    0.000000    57342.000000    0.000000    1.642979    0.000000    0.000000    0.000000    0.000000    0.000000    46504.500000    192.000000    245.000000    0.000000    4051.500000    0.000000    0.000000    0.000000    0.000000    0.000000    0.134507    0.030507    0.044374    0.000000    0.059627    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    91.197811    0.376522    0.480458    0.000000    7.945208    0.999157    0.000000    0.000000    0.000000    0.000000    0.000000    nan    0.543212    6.522464    0.455945    92.478379    0.295613    0.065337
3    0.985600    2    3490123.000000    0.000000    34167.000000    0.000000    0.978963    0.000000    0.000000    0.000000    0.000000    0.000000    46470.500000    178.500000    210.000000    0.000000    4134.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.051307    0.029120    0.030507    0.000000    0.008320    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    91.131136    0.350048    0.411821    0.000000    8.106995    0.890318    0.000000    0.000000    0.000000    0.000000    0.000000    nan    0.498108    6.586198    0.392211    92.523484    0.282486    0.056203
4    0.985600    2    3490123.000000    0.000000    24896.500000    0.000000    0.713342    0.000000    0.000000    0.000000    0.000000    0.000000    46392.500000    166.500000    180.000000    0.000000    4254.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.059627    0.029120    0.044374    0.000000    0.013867    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    90.978173    0.326515    0.352990    0.000000    8.342321    0.809915    0.000000    0.000000    0.000000    0.000000    0.000000    nan    0.455945    6.624439    0.353970    92.565646    0.279628    0.050724
5    0.985600    2    3490123.000000    0.000000    19778.500000    0.000000    0.566699    0.000000    0.000000    0.000000    0.000000    0.000000    46335.500000    151.500000    180.500000    0.000000    4325.500000    0.000000    0.000000    0.000000    0.000000    0.000000    0.065174    0.006933    0.009707    0.000000    0.048534    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    90.866393    0.297100    0.353970    0.000000    8.482537    0.759908    0.000000    0.000000    0.000000    0.000000    0.000000    nan    0.429471    6.647971    0.330438    92.592120    0.277824    0.047351

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      2
trainlines      3490123.000000
rules           57342.000000
rules%         1.642979
same%stdev     0.134507
ambi1%stdev    0.030507
ambi2%stdev    0.044374
ambi3%stdev    0.000000
diff%stdev     0.059627
same%          91.197811
ambi1%         0.376522
ambi2%         0.480458
ambi3%         0.000000
diff%          7.945208
amb.rules%     0.999157
false_amb%     0.543212
false_not_amb% 6.522464
true_amb%      0.455945
true_not_amb%  92.478379
precision       0.295613
recall          0.065337

Comment: If comparing by diff%, comp_peen is marginally worse than best_pl 
                             (best_pl uses automatically computed parameters)
If compared by same%, comp_peen is 0.3% better. Reason: best_pl has many more
ambiguous rules and generates more false ambiguous results, but also more true
ambiguous results.
*/
    int A1 = 3*(a->W__R        + a->R__R) + a->R__NA;
    // 3wr + 3rr + rn - r = 3wr + 2rr - rw
    int B1 = 3*(b->W__R        + b->R__R) + b->R__NA;
    int A2 = a->W__R            + a->R__R + a->R__NA;
    //wr + rr + rn - r = wr - rw
    int B2 = b->W__R            + b->R__R + b->R__NA;
    int A3 = a->W__R + a->W__NA           + a->R__NA;
    // wr + wn + rn - w = -ww + rn
    int B3 = b->W__R + b->W__NA           + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif
#if _NA
static int comp_sugar(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    // Slovene 86.273367 +/- 0.410931 at 0.9856 of dataset, 9 iterations, 17254.777778 = 8.819297% rules (cutoff = 1)
    // BEST English 89.060606 +/- 1.320829 at 0.9856 of dataset, 3 iterations, 1318.266667 = 1.752199% rules, cutoff=2
    // Icelandic 70.925172 +/- 1.858255 at 0.9856 of dataset, 17 iterations, 27151.294118 = 47.351402% rules, cutoff = 0
    // Dutch 90.172822 +/- 0.307911 at 0.9856 of dataset, 7 iterations, 57761.142857 = 19.411791% rules, cutoff = 0
    // BEST Greek 90.422464 +/- 0.437009 at 0.9856 of dataset, 5 iterations, 132765.6 = 24.535334% rules, cutoff = 0
    // BEST German 91.461918 +/- 0.167574 at 0.9856 of dataset, 7 iterations, 50986 = 16.405554% rules, cutoff = 0
    // BEST Swedish 92.265969 +/- 0.277289 at 0.9856 of dataset, 6 iterations, 25935.333333 = 5.506008% rules, cutoff = 1
    // Norwegian 86.665700 +/- 0.676264 at 0.9856 of dataset, 6 iterations, 46685.5 = 9.890849% rules, cutoff = 1
    // Danish 92.585623 +/- 0.171327 at 0.9856 of dataset, 5 iterations, 30422.400000 = 5.576679% rules, cutoff = 1
    // BEST Russian 80.815622 +/- 0.450500 at 0.9856 of dataset, 6 iterations, 47079.166667 = 12.084440% rules, cutoff = 1        
//AMBI:
    // French ok 75.472316 ambi1 4.615600 ambi2 3.493266 diff 16.418818 rules 4162.909091 3.129560% cutoff 2
#if 1
    // next lines from affixFEW2
    int A1 = a->W__R           + 2*a->R__R + a->R__NA; // good: all words that are lemmatised correctly. bad: all previously right words that got it wrong
    // wr - rw + rr
    int B1 = b->W__R           + 2*b->R__R + b->R__NA;
    int A2 = a->W__R           +   a->R__R + a->R__NA;
    // wr - rw
    int B2 = b->W__R           +   b->R__R + b->R__NA;
    // next lines from peen
    int A3 = a->W__R + a->W__NA           + a->R__NA;
    // -ww + rn 
    int B3 = b->W__R + b->W__NA           + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
#else
    //equivalent with:
    int AA1 = a->W__R - a->R__W + a->R__R;
    int AA2 = - a->R__R;
    int AA3 = - a->W__R - 2*a->W__W;//a->R__NA - a->W__W;
    int BB1 = b->W__R - b->R__W + b->R__R;
    int BB2 = - b->R__R;
    int BB3 = - b->W__R - 2*b->W__W;//b->R__NA - b->W__W;
#endif
    }
#endif

static int comp_honey(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    // (OK) Dutch 90.179393 +/- 0.589662 at 0.9856 of dataset, 7 iterations, 73324.571429 = 24.642193% rules, cutoff = 0
    // (OK) Norwegian 87.272244 +/- 0.267729 at 0.9856 of dataset, 6 iterations, 141038.666667 = 29.880630% rules, cutoff = 0
    // (OK) English 88.315152 +/- 1.097312 at 0.9856 of dataset, 3 iterations, 5285.466667 = 7.025276% rules, cutoff=1
    // (OK) Icelandic 70.742665 +/- 1.686147 at 0.9856 of dataset, 17 iterations, 29857.000000 = 52.070108% rules, cutoff = 0
    // (?) Slovene 86.273367 +/- 0.410931 at 0.9856 of dataset, 9 iterations, 17254.777778 = 8.819297% rules (cutoff = 1)
    // (?) BEST Greek 90.422464 +/- 0.437009 at 0.9856 of dataset, 5 iterations, 132765.6 = 24.535334% rules, cutoff = 0
    // (?) BEST German 91.461918 +/- 0.167574 at 0.9856 of dataset, 7 iterations, 50986 = 16.405554% rules, cutoff = 0
    // (?) BEST Swedish 92.265969 +/- 0.277289 at 0.9856 of dataset, 6 iterations, 25935.333333 = 5.506008% rules, cutoff = 1
    // (?) Danish 92.585623 +/- 0.171327 at 0.9856 of dataset, 5 iterations, 30422.400000 = 5.576679% rules, cutoff = 1
    // (?) BEST Russian 80.815622 +/- 0.450500 at 0.9856 of dataset, 6 iterations, 47079.166667 = 12.084440% rules, cutoff = 1        
//AMBI:
    // French ok 84.477707 ambi1 2.251592 ambi2 1.426752 diff 11.843949 rules 7413.875000 2.760295% cutoff 2
    int A1 = a->W__R           + 2*a->R__R;
    int B1 = b->W__R           + 2*b->R__R;
    int A2 = a->W__R           +   a->R__R;
    int B2 = b->W__R           +   b->R__R;
    int A3 = a->W__R ;
    int B3 = b->W__R ;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }

#if _NA
static int comp_beet(const vertex * a,const vertex * b)
    {
    //const vertex * a = *(const vertex **)A;
    //const vertex * b = *(const vertex **)B;
    //Icelandic 71.034094 at 0.939 of dataset
//AMBI:
    // French ok 85.057325 ambi1 1.283439 ambi2 1.057325 diff 12.601911 rules 7260.375000 2.703144% cutoff 2
    int A1 = 2*(a->W__R        + a->R__R) + a->R__NA;
    int B1 = 2*(b->W__R        + b->R__R) + b->R__NA;
    int A2 = a->W__R            + a->R__R + a->R__NA;
    int B2 = b->W__R            + b->R__R + b->R__NA;
    int A3 = a->W__R + a->W__NA           + a->R__NA;
    int B3 = b->W__R + b->W__NA           + b->R__NA;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }
#endif

static int comp_koud(const vertex * a,const vertex * b)
    {
    // German 91.260578 +/- 0.363285 at 0.9856 of dataset, 7 iterations, 30890.714286 = 9.939577% rules, cutoff = 0
//AMBI:
    // French 86.356688 ambi1 0.996815 ambi2 0.796178 diff 11.850318 rules 3335.250000 1.241763% cutoff 3 (!)
    // French 85.250493 ambi1 2.333057 ambi2 2.161181 diff 10.255268 rules 28520.250000 10.618597% cutoff 0 (!) paradigms+homographs clumped
    // French 85.313973 ambi1 2.050694 ambi2 2.289517 diff 10.345816 rules 28509.250000 10.614432% cutoff 0 (!) homographs clumped
    // Dutch.clumped.ph 85.789838 ambi1 1.086067 ambi2 1.256060 diff 11.868035 rules 37400.142857 12.190637% cutoff 0 paradigms+homographs clumped
    // Dutch.clumped.h 85.818923 ambi1 1.095507 ambi2 1.060476 diff 12.025095 rules 37411.857143 12.192383% cutoff 0 homographs clumped
    // (Dutch.clumped.ph suffix, old algo:
    //                  83.532708 ambi1 1.948624 ambi2 2.719889 ambi3 0.107033 diff 11.691746 rules 73024.571429 23.802477% cutoff 0 paradigms+homographs clumped
    // (Dutch.clumped.h suffix, old algo:
    //                  83.624725 ambi1 1.859813 ambi2 2.611382 ambi3 0.162415 diff 11.741664 rules 72975.428571 23.782417% cutoff 0 paradigms+homographs clumped
    // Russian clumped ph 74.983460 ambi1 0.517762 ambi2 0.558033 diff 23.940745 rules 95077.500000 24.184389% cutoff 0 paradigms+homographs clumped
    // (old algo:)
    //                  79.485114 ambi1 0.218611 ambi2 0.342298 ambi3 0.005753 diff 19.948224 rules 94247.166667 23.973181% cutoff 0 paradigms+homographs clumped
    // The A1 vs B1 condition is pretty close to what later was found as the
    // best factors using automatic factor setting (comp_parms).
    // These factors were found by manual optimizing.
    int A1 = 6*a->W__R - 5*a->R__W + a->W__W;
    int B1 = 6*b->W__R - 5*b->R__W + b->W__W;
    int A2 = a->W__R - 6*a->R__W;
    int B2 = b->W__R - 6*b->R__W;
    int A3 = a->R__R - a->W__W;
    int B3 = b->R__R - b->W__W;
    return (A1>B1)?-1:(A1<B1)?1:(A2>B2)?-1:(A2<B2)?1:(A3>B3)?-1:(A3<B3)?1:0;
    }

int (*comp)(const vertex * a,const vertex * b) = comp_parms;
// returns b > a ? 1 : b < a ? -1 : 0
// (Chosen like this to let qsort sort in descending order.)


// You can find a local optimum for the parameters by using comp_parms as the
// weight function and setting compute_parms = true. The parameters parms.Matrix[] 
// can be seeded with non-zero values by hard coding. The file parms.txt
// will contain the currently best parameters.

// Optimal parameters == resulting in smallest rule set.
// Hypothesis: small rule sets give (almost) best lemmatization results.
// Optimizing for the size of rule sets is computationally MUCH cheaper
// than optimizing for accuracy.

// If you have found a good set of parameters (presumably with a small subset
// of the training data), you can hard code them (as is done below) and run 
// the program with the full set of the training data. In that case,
// set compute_parms = false

//bool compute_parms = false;


struct rotation
    {
    double Matrix[6];
    int ROWPARMS;
    // metadata:
    const char * Language;
    const char * Xparm;
    int NumberOfRules;
    int CutOffForFewestErrors;
    double FractionErroneousPedictionOOV;
    void init(optionStruct * options)
        {
        ROWPARMS = options->numberOfParms();
        for(int i = 0;i < ROWPARMS;++i)
            {
            Matrix[i] = options->parm(i);
            }
        }
    void better(optionStruct * options)
        {
        for(int i = 0;i < ROWPARMS;++i)
            {
            options->setParm(i,Matrix[i]);
            }
        }
    } rotation;

static struct rotation parms = 
   /* R_R   W_R   R_W   W_W  R_NA  W_NA */   
    {{  0.0,  3.0, -2.0,  1.0,  0.0,  0.0},6,"","",0,0,0.0}
    ;

static void normalise(double * ROW)
    {
    double modulus = 0.0;
    for(int i = 0;i < parms.ROWPARMS;++i)
        modulus += ROW[i] * ROW[i];
    modulus = sqrt(modulus);
    for(int i = 0;i < parms.ROWPARMS;++i)
        ROW[i] /= modulus;
    }

static double inner(double * a, double * b)
    {
    double ret = 0;
    for(int i = 0;i < parms.ROWPARMS;++i)
        ret += a[i]*b[i];
    return ret;
    }

static void times(double * a, double f)
    {
    for(int i = 0;i < parms.ROWPARMS;++i)
        a[i] *= f;
    }


struct bestParms
    {
    bool suffixonly;
    const char * langbase;
    int rowss;
    struct rotation val;
    // Each row:
    // R__R W__R R__W W__W R__NA W__NA
    // Generally, good that Wrongs change to Rights (W__R > 0) and that Rights don't change to Wrongs (R__W < 0)
    // But what about rules that don't improve lemmatisation? (R__R > 0 or W__W > 0)
    // Intuitively difficult to decide!
    };







#if 1
static bestParms best_is_suffix =
    {
    true,
    "is",
    1,
//iteration:18.1
/*weight (not  used): 1.41244386452166131e+05 suffix only: yes */
/* number of nodes: 152108, nodes/line: 1.05629895368709495e-01 weight (not  used): 1.41244386452166131e+05 blobs 2809220 lines 2873370 * fraction 5.01187233627272799e-01 = 1440009 lines*/
        {{
        0.00000000000000000e+00,    6.94542434383270568e-01,    -7.18112257666929654e-01,    4.38815704990783637e-02
        }}
    };
#elif 1
static bestParms best_is_suffix =
    {
    true,
    "is",
    1,
/*
0    0.985600    2    2831993.000000    471048.500000    306036.000000    16.633110    10.806383    36063.500000    712.500000    528.500000    83.000000    3989.500000    36948.000000    188.500000    140.500000    0.000000    4100.000000    0.176020    0.158931    0.029052    0.006836    0.052977    0.102536    0.039306    0.042723    0.000000    0.105954    87.158325    1.721971    1.277280    0.200595    9.641830    3.905551    89.295986    0.455567    0.339561    0.000000    9.908887    1.377577    1.841603    33.482369    2.063949    62.612079    0.359125    0.058064    0.424149    34.592890    0.953428    64.029533    0.529175    0.026822
1    0.985600    2    2831993.000000    128789.000000    168125.000000    4.547645    5.936632    36520.000000    0.000000    0.000000    0.000000    4857.000000    37099.500000    273.500000    263.000000    0.000000    3741.000000    0.051268    0.000000    0.000000    0.000000    0.051268    0.032470    0.001709    0.047850    0.000000    0.082029    88.261595    0.000000    0.000000    0.000000    11.738405    0.000000    89.662131    0.660995    0.635619    0.000000    9.041255    1.598714    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.493028    34.440631    1.105687    63.960654    0.528596    0.031106
2    0.985600    2    2831993.000000    84003.000000    52867.500000    2.966215    1.866795    36089.500000    0.000000    0.000000    0.000000    5287.500000    37338.500000    236.500000    207.000000    0.000000    3595.000000    0.063231    0.000000    0.000000    0.000000    0.063231    0.029052    0.049559    0.037597    0.000000    0.041014    87.221162    0.000000    0.000000    0.000000    12.778838    0.000000    90.239747    0.571574    0.500278    0.000000    8.688402    1.203567    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.420523    34.763274    0.783044    64.033159    0.482143    0.022029
3    0.985600    2    2831993.000000    63704.000000    38416.000000    2.249441    1.356501    35649.000000    0.000000    0.000000    0.000000    5728.000000    37204.000000    219.500000    209.500000    0.000000    3744.000000    0.017089    0.000000    0.000000    0.000000    0.017089    0.088865    0.022216    0.025634    0.000000    0.085447    86.156560    0.000000    0.000000    0.000000    13.843440    0.000000    89.914687    0.530488    0.506320    0.000000    9.048505    1.185441    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.422940    34.783817    0.762501    64.030742    0.474080    0.021451
4    0.985600    2    2831993.000000    51623.000000    31338.500000    1.822851    1.106588    35334.000000    0.000000    0.000000    0.000000    6043.000000    37065.500000    224.000000    188.000000    0.000000    3899.500000    0.061522    0.000000    0.000000    0.000000    0.061522    0.097409    0.006836    0.023925    0.000000    0.080320    85.395268    0.000000    0.000000    0.000000    14.604732    0.000000    89.579960    0.541364    0.454359    0.000000    9.424318    1.152814    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.435024    34.828528    0.717790    64.018658    0.452055    0.020193
5    0.985600    2    2831993.000000    43987.500000    27450.500000    1.553235    0.969300    35034.500000    0.000000    0.000000    0.000000    6342.500000    36936.000000    219.500000    179.000000    0.000000    4042.500000    0.093992    0.000000    0.000000    0.000000    0.093992    0.082029    0.015380    0.003418    0.000000    0.100827    84.671436    0.000000    0.000000    0.000000    15.328564    0.000000    89.266984    0.530488    0.432607    0.000000    9.769920    1.125021    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.408439    34.829736    0.716582    64.045243    0.467297    0.020159

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          2
fraction        9856.000000
iterations      2
trainlines      2831993.000000
rules           52867.500000 (84003.000000)
rules%         1.866795 (2.966215)
same%stdev     0.029052
ambi1%stdev    0.049559
ambi2%stdev    0.037597
ambi3%stdev    0.000000
diff%stdev     0.041014
same%          90.239747 (87.221162)
ambi1%         0.571574 (0.000000)
ambi2%         0.500278 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          8.688402 (12.778838)
amb.rules%     1.203567 (0.000000)
false_amb%     0.420523 (0.000000)
false_not_amb% 34.763274 (35.546318)
true_amb%      0.783044 (0.000000)
true_not_amb%  64.033159 (64.453682)
precision       0.482143 (0.000000)
recall          0.022029 (0.000000)

bests[11].suffixonly == [true]
bests[11].langbase == [is]
comp = comp_parms0_off
bests[11].rows == [1]
  R->R  W->R  R->W  W->W

0.0085050.669878-0.7383730.077434
*/
//iteration:14.-1
/* number of nodes: 88858, nodes/line: 1.20004915909585078e-01 weight: 8.44399637102287234e+04 blobs 1 lines 5881633 * fraction 1.25892541179416839e-01 = 740453 lines*/
        {                                                 // # decisions
        8.50547688621742723e-03,    6.69877760720549498e-01,    -7.38373491250877478e-01,    7.74340362692699236e-02, //1177883
        -7.82684292973299223e-01,    4.59948960180274258e-01,    3.92892970820137744e-01,    -1.46585691805502294e-01, //0
        -6.01147073676968957e-01,    -5.82239786374368462e-01,    -5.42665925463000409e-01,    -7.16430060350067843e-02, //0
        -1.61106021629986690e-01,    -2.66002045826051053e-02,    7.71582415526826937e-02,    9.83556752135442247e-01  //0
        }                                                 //(0 unresolved comparisons)
// Same as
//iteration:13.11
/* number of nodes: 74744, nodes/line: 1.42586526923832640e-01 weight: 7.11594670680982090e+04 blobs 1 lines 5881633 * fraction 8.91250938133746201e-02 = 524201 lines*/
    };
#endif

static bestParms best_isC0 =
    {
    false,
    "isC0",
    1,
/* SINGLE SHOT
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          2 2831993.000000  471048.500000   60440.000000      16.633110       2.134186   36063.500000     712.500000     528.500000      83.000000    3989.500000   36357.500000     343.000000     254.000000       0.000000    4422.500000       0.176020       0.158931       0.029052       0.006836       0.052977       0.227289       0.092283       0.003418       0.000000       0.131588      87.158325       1.721971       1.277280       0.200595       9.641830       3.905551      87.868864       0.828963       0.613868       0.000000      10.688305       2.130411       1.841603      33.482369       2.063949      62.612079       0.359125       0.058064       0.775793      34.191701       1.354617      63.677889       0.466112       0.038109
1       0.985600          2 2831993.000000  128789.000000   32445.000000       4.547645       1.145660   36520.000000       0.000000       0.000000       0.000000    4857.000000   36779.500000     343.500000     318.000000       0.000000    3936.000000       0.051268       0.000000       0.000000       0.000000       0.051268       0.278557       0.083738       0.064940       0.000000       0.259758      88.261595       0.000000       0.000000       0.000000      11.738405       0.000000      88.888755       0.830171       0.768543       0.000000       9.512531       1.969693       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.679121      34.255746       1.290572      63.774561       0.487226       0.036307
2       0.985600          2 2831993.000000   84003.000000   14192.000000       2.966215       0.501131   36089.500000       0.000000       0.000000       0.000000    5287.500000   36974.500000     310.500000     289.500000       0.000000    3802.500000       0.063231       0.000000       0.000000       0.000000       0.063231       0.295646       0.090574       0.029052       0.000000       0.234124      87.221162       0.000000       0.000000       0.000000      12.778838       0.000000      89.360031       0.750417       0.699664       0.000000       9.189888       1.688136       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.607826      34.466008       1.080310      63.845856       0.470526       0.030392
3       0.985600          2 2831993.000000   63704.000000   10130.000000       2.249441       0.357699   35649.000000       0.000000       0.000000       0.000000    5728.000000   36863.000000     286.500000     274.000000       0.000000    3953.500000       0.017089       0.000000       0.000000       0.000000       0.017089       0.283684       0.083738       0.013671       0.000000       0.213617      86.156560       0.000000       0.000000       0.000000      13.843440       0.000000      89.090558       0.692414       0.662204       0.000000       9.554825       1.579380       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.566740      34.533678       1.012640      63.886942       0.471847       0.028488
4       0.985600          2 2831993.000000   51623.000000    8154.000000       1.822851       0.287924   35334.000000       0.000000       0.000000       0.000000    6043.000000   36739.500000     269.500000     257.500000       0.000000    4110.500000       0.061522       0.000000       0.000000       0.000000       0.061522       0.275139       0.029052       0.008545       0.000000       0.237542      85.395268       0.000000       0.000000       0.000000      14.604732       0.000000      88.792083       0.651328       0.622326       0.000000       9.934263       1.498417       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.547406      34.595307       0.951011      63.906276       0.464855       0.026754
5       0.985600          2 2831993.000000   43987.500000    6904.000000       1.553235       0.243786   35034.500000       0.000000       0.000000       0.000000    6342.500000   36631.500000     255.000000     239.500000       0.000000    4251.000000       0.093992       0.000000       0.000000       0.000000       0.093992       0.305900       0.017089       0.015380       0.000000       0.273430      84.671436       0.000000       0.000000       0.000000      15.328564       0.000000      88.531068       0.616284       0.578824       0.000000      10.273824       1.407787       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.526863      34.665394       0.880924      63.926819       0.455340       0.024782
cutoff 0 Affix  a      11.009406 b       0.000000: N(rules)=  60440.000000*N(trainpairs)^0.000000
         Suffix a      -0.112324 b       0.878094: N(rules)=      0.893755*N(trainpairs)^0.878094
cutoff 1 Affix  a      10.387302 b       0.000000: N(rules)=  32445.000000*N(trainpairs)^0.000000
         Suffix a      -0.733958 b       0.836893: N(rules)=      0.480005*N(trainpairs)^0.836893
cutoff 2 Affix  a       9.560434 b       0.000000: N(rules)=  14192.000000*N(trainpairs)^0.000000
         Suffix a      -1.605792 b       0.868583: N(rules)=      0.200730*N(trainpairs)^0.868583
cutoff 3 Affix  a       9.223257 b       0.000000: N(rules)=  10130.000000*N(trainpairs)^0.000000
         Suffix a      -2.263068 b       0.896039: N(rules)=      0.104031*N(trainpairs)^0.896039
cutoff 4 Affix  a       9.006264 b       0.000000: N(rules)=   8154.000000*N(trainpairs)^0.000000
         Suffix a      -2.871889 b       0.924806: N(rules)=      0.056592*N(trainpairs)^0.924806
cutoff 5 Affix  a       8.839856 b       0.000000: N(rules)=   6904.000000*N(trainpairs)^0.000000
         Suffix a      -3.315163 b       0.944783: N(rules)=      0.036328*N(trainpairs)^0.944783

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        no
cutoff                      2
fraction          9856.000000
iterations                  2
trainlines     2831993.000000
rules            14192.000000 (  84003.000000)
rules%               0.501131 (      2.966215)
same%stdev           0.295646
ambi1%stdev          0.090574
ambi2%stdev          0.029052
ambi3%stdev          0.000000
diff%stdev           0.234124
same%               89.360031 (     87.221162)
ambi1%               0.750417 (      0.000000)
ambi2%               0.699664 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%                9.189888 (     12.778838)
amb.rules%           1.688136 (      0.000000)
false_amb%           0.607826 (      0.000000)
false_not_amb%      34.466008 (     35.546318)
true_amb%            1.080310 (      0.000000)
true_not_amb%       63.845856 (     64.453682)
precision            0.470526 (      0.000000)
recall               0.030392 (      0.000000)

bests[16].suffixonly == [false]
bests[16].langbase == [isC0]
comp = comp_parms0_off
bests[16].rows == [1]
  R->R     W->R     R->W     W->W

0.247669 0.636764 -0.729558 0.022230 -0.004214 0.020625 
*/
/* REDO
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          2 2831993.000000  471048.500000   60440.000000      16.633110       2.134186   36063.500000     712.500000     528.500000      83.000000    3989.500000   38009.500000     563.500000     543.000000       0.000000    2261.000000       0.176020       0.158931       0.029052       0.006836       0.052977       5.419039       0.845924       0.991183       0.000000       7.256146      87.158325       1.721971       1.277280       0.200595       9.641830       3.905551      91.861421       1.361868       1.312323       0.000000       5.464388       3.111632       1.841603      33.482369       2.063949      62.612079       0.359125       0.058064       0.374604      32.809290       2.737028      64.079078       0.785095       0.076999
1       0.985600          2 2831993.000000  128789.000000   32445.000000       4.547645       1.145660   36520.000000       0.000000       0.000000       0.000000    4857.000000   37620.000000     466.500000     460.000000       0.000000    2830.500000       0.051268       0.000000       0.000000       0.000000       0.051268       2.594166       0.504136       0.420398       0.000000       3.518701      88.261595       0.000000       0.000000       0.000000      11.738405       0.000000      90.920076       1.127438       1.111729       0.000000       6.840757       2.505015       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.546197      33.587500       1.958818      63.907485       0.641980       0.055106
2       0.985600          2 2831993.000000   84003.000000   14192.000000       2.966215       0.501131   36089.500000       0.000000       0.000000       0.000000    5287.500000   37357.500000     328.500000     317.500000       0.000000    3373.500000       0.063231       0.000000       0.000000       0.000000       0.063231       1.013400       0.152095       0.066649       0.000000       1.232143      87.221162       0.000000       0.000000       0.000000      12.778838       0.000000      90.285666       0.793919       0.767335       0.000000       8.153080       1.789642       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.541364      34.298040       1.248278      63.912318       0.535511       0.035117
3       0.985600          2 2831993.000000   63704.000000   10130.000000       2.249441       0.357699   35649.000000       0.000000       0.000000       0.000000    5728.000000   37147.500000     295.000000     281.500000       0.000000    3653.000000       0.017089       0.000000       0.000000       0.000000       0.017089       0.688702       0.112790       0.011963       0.000000       0.813454      86.156560       0.000000       0.000000       0.000000      13.843440       0.000000      89.778138       0.712956       0.680330       0.000000       8.828576       1.604756       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.513570      34.455132       1.091186      63.940112       0.515117       0.030698
4       0.985600          2 2831993.000000   51623.000000    8154.000000       1.822851       0.287924   35334.000000       0.000000       0.000000       0.000000    6043.000000   36962.000000     285.000000     261.500000       0.000000    3868.500000       0.061522       0.000000       0.000000       0.000000       0.061522       0.485338       0.082029       0.022216       0.000000       0.589583      85.395268       0.000000       0.000000       0.000000      14.604732       0.000000      89.329821       0.688788       0.631994       0.000000       9.349397       1.540711       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.525654      34.531261       1.015057      63.928028       0.491228       0.028556
5       0.985600          2 2831993.000000   43987.500000    6904.000000       1.553235       0.243786   35034.500000       0.000000       0.000000       0.000000    6342.500000   36809.500000     273.000000     235.500000       0.000000    4059.000000       0.093992       0.000000       0.000000       0.000000       0.093992       0.302482       0.078611       0.001709       0.000000       0.382802      84.671436       0.000000       0.000000       0.000000      15.328564       0.000000      88.961259       0.659787       0.569157       0.000000       9.809798       1.437997       0.000000      35.546318       0.000000      64.453682       0.000000       0.000000       0.506320      34.614641       0.931677      63.947362       0.479180       0.026210
cutoff 0 Affix  a      11.009406 b       0.000000: N(rules)=  60440.000000*N(trainpairs)^0.000000
         Suffix a      -0.112324 b       0.878094: N(rules)=      0.893755*N(trainpairs)^0.878094
cutoff 1 Affix  a      10.387302 b       0.000000: N(rules)=  32445.000000*N(trainpairs)^0.000000
         Suffix a      -0.733958 b       0.836893: N(rules)=      0.480005*N(trainpairs)^0.836893
cutoff 2 Affix  a       9.560434 b       0.000000: N(rules)=  14192.000000*N(trainpairs)^0.000000
         Suffix a      -1.605792 b       0.868583: N(rules)=      0.200730*N(trainpairs)^0.868583
cutoff 3 Affix  a       9.223257 b       0.000000: N(rules)=  10130.000000*N(trainpairs)^0.000000
         Suffix a      -2.263068 b       0.896039: N(rules)=      0.104031*N(trainpairs)^0.896039
cutoff 4 Affix  a       9.006264 b       0.000000: N(rules)=   8154.000000*N(trainpairs)^0.000000
         Suffix a      -2.871889 b       0.924806: N(rules)=      0.056592*N(trainpairs)^0.924806
cutoff 5 Affix  a       8.839856 b       0.000000: N(rules)=   6904.000000*N(trainpairs)^0.000000
         Suffix a      -3.315163 b       0.944783: N(rules)=      0.036328*N(trainpairs)^0.944783

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        yes
cutoff                      0
fraction          9856.000000
iterations                  2
trainlines     2831993.000000
rules            60440.000000 ( 471048.500000)
rules%               2.134186 (     16.633110)
same%stdev           5.419039
ambi1%stdev          0.845924
ambi2%stdev          0.991183
ambi3%stdev          0.000000
diff%stdev           7.256146
same%               91.861421 (     87.158325)
ambi1%               1.361868 (      1.721971)
ambi2%               1.312323 (      1.277280)
ambi3%               0.000000 (      0.200595)
diff%                5.464388 (      9.641830)
amb.rules%           3.111632 (      3.905551)
false_amb%           0.374604 (      1.841603)
false_not_amb%      32.809290 (     33.482369)
true_amb%            2.737028 (      2.063949)
true_not_amb%       64.079078 (     62.612079)
precision            0.785095 (      0.359125)
recall               0.076999 (      0.058064)

bests[16].suffixonly == [false]
bests[16].langbase == [isC0]
comp = comp_parms0_off
bests[16].rows == [1]
  R->R     W->R     R->W     W->W

0.247669 0.636764 -0.729558 0.022230 -0.004214 0.020625 
*/
//iteration:20.6
/*weight ( used): 1.05436295090904787e+04 suffix only: no */
/* number of nodes: 336797, nodes/line: inf weight ( used): 1.05436295090904787e+04 blobs 1 lines 0 * fraction 1.00000000000000000e+00 = 0 lines*/
        {{
        2.47669087481595079e-01,    6.36764047976876468e-01,    -7.29557569755324042e-01,    2.22303428808458027e-02,    -4.21447897318842114e-03,    2.06245665463890698e-02
        }}
    };



#if 1
static bestParms best_is =
    {
    false,
    "is",
    1,
//iteration:18.2
/*weight (not  used): 1.34340843669173279e+05 suffix only: no */
/* number of nodes: 145852, nodes/line: 1.01285478076873131e-01 weight (not  used): 1.34340843669173279e+05 blobs 2809220 lines 2873370 * fraction 5.01187233627272799e-01 = 1440009 lines*/
        {{
        0.00000000000000000e+00,    6.96451349087997107e-01,    -7.13849249589145862e-01,    7.33128038921041919e-02
        }}
    };
#elif 1
static bestParms best_is =
    {
    false,
    "is",
    1,
//iteration:12.0
/*weight (not  used): 3.09072013007138085e+04 suffix only: no */
/* number of nodes: 32528, nodes/line: 1.79481661728272457e-01 weight (not  used): 3.09072013007138085e+04 blobs 2809220 lines 2873370 * fraction 6.30957344480193721e-02 = 181233 lines*/
        {
        0.00000000000000000e+00,        7.04722964067779345e-01,        -7.05629153519126029e-01,       7.38447128737378389e-02
        }
/*
0    0.985600    2    2831993.000000    471048.500000    279488.000000    16.633110    9.868951    36063.500000    712.500000    528.500000    83.000000    3989.500000    36331.500000    392.000000    307.500000    0.000000    4346.000000    0.176020    0.158931    0.029052    0.006836    0.052977    0.312735    0.047850    0.025634    0.000000    0.239251    87.158325    1.721971    1.277280    0.200595    9.641830    3.905551    87.806028    0.947386    0.743166    0.000000    10.503420    2.399884    3.385939    3.707374    0.519612    92.387075    0.071263    0.122927    1.883897    3.710999    0.515987    93.889117    0.120451    0.122070
cutoff 0    Affix a 1.518217 b 0.742218    
cutoff 0    Suffix a -0.112324 b 0.878094    
1    0.985600    2    2831993.000000    128789.000000    145060.000000    4.547645    5.122188    36520.000000    0.000000    0.000000    0.000000    4857.000000    36701.000000    424.000000    408.000000    0.000000    3844.000000    0.051268    0.000000    0.000000    0.000000    0.051268    0.218744    0.092283    0.068357    0.000000    0.194819    88.261595    0.000000    0.000000    0.000000    11.738405    0.000000    88.699036    1.024724    0.986055    0.000000    9.290185    2.352756    0.000000    4.226986    0.000000    95.773014    0.000000    0.000000    1.678469    3.552698    0.674288    94.094545    0.167266    0.159520
cutoff 1    Affix a 1.181955 b 0.721542    
cutoff 1    Suffix a -0.733958 b 0.836893    
2    0.985600    2    2831993.000000    84003.000000    55547.000000    2.966215    1.961410    36089.500000    0.000000    0.000000    0.000000    5287.500000    36932.500000    351.500000    319.500000    0.000000    3773.500000    0.063231    0.000000    0.000000    0.000000    0.063231    0.124752    0.008545    0.052977    0.000000    0.169185    87.221162    0.000000    0.000000    0.000000    12.778838    0.000000    89.258525    0.849506    0.772168    0.000000    9.119801    1.813810    0.000000    4.226986    0.000000    95.773014    0.000000    0.000000    1.340116    3.753293    0.473693    94.432898    0.150192    0.112064
cutoff 2    Affix a 0.232484 b 0.725173    
cutoff 2    Suffix a -1.605792 b 0.868583    
3    0.985600    2    2831993.000000    63704.000000    38924.500000    2.249441    1.374456    35649.000000    0.000000    0.000000    0.000000    5728.000000    36852.000000    316.500000    290.000000    0.000000    3918.500000    0.017089    0.000000    0.000000    0.000000    0.017089    0.143551    0.008545    0.068357    0.000000    0.203363    86.156560    0.000000    0.000000    0.000000    13.843440    0.000000    89.063973    0.764918    0.700872    0.000000    9.470237    1.657926    0.000000    4.226986    0.000000    95.773014    0.000000    0.000000    1.236194    3.805254    0.421732    94.536820    0.145720    0.099771
cutoff 3    Affix a -0.574206 b 0.758604    
cutoff 3    Suffix a -2.263068 b 0.896039    
4    0.985600    2    2831993.000000    51623.000000    31706.000000    1.822851    1.119565    35334.000000    0.000000    0.000000    0.000000    6043.000000    36747.500000    303.000000    251.000000    0.000000    4075.500000    0.061522    0.000000    0.000000    0.000000    0.061522    0.152095    0.017089    0.044432    0.000000    0.179438    85.395268    0.000000    0.000000    0.000000    14.604732    0.000000    88.811417    0.732291    0.606617    0.000000    9.849675    1.534669    0.000000    4.226986    0.000000    95.773014    0.000000    0.000000    1.155231    3.847548    0.379438    94.617783    0.141060    0.089766
cutoff 4    Affix a -1.201505 b 0.788680    
cutoff 4    Suffix a -2.871889 b 0.924806    
5    0.985600    2    2831993.000000    43987.500000    27360.000000    1.553235    0.966104    35034.500000    0.000000    0.000000    0.000000    6342.500000    36633.500000    293.500000    228.500000    0.000000    4221.500000    0.093992    0.000000    0.000000    0.000000    0.093992    0.073484    0.049559    0.052977    0.000000    0.176020    84.671436    0.000000    0.000000    0.000000    15.328564    0.000000    88.535902    0.709331    0.552239    0.000000    10.202528    1.451289    0.000000    4.226986    0.000000    95.773014    0.000000    0.000000    1.096020    3.871716    0.355270    94.676994    0.139469    0.084048
cutoff 5    Affix a -1.734091 b 0.815366    
cutoff 5    Suffix a -3.315163 b 0.944783    
New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      2
trainlines      2831993.000000
rules           55547.000000 (84003.000000)
rules%         1.961410 (2.966215)
same%stdev     0.124752
ambi1%stdev    0.008545
ambi2%stdev    0.052977
ambi3%stdev    0.000000
diff%stdev     0.169185
same%          89.258525 (87.221162)
ambi1%         0.849506 (0.000000)
ambi2%         0.772168 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          9.119801 (12.778838)
amb.rules%     1.813810 (0.000000)
false_amb%     1.340116 (0.000000)
false_not_amb% 3.753293 (4.226986)
true_amb%      0.473693 (0.000000)
true_not_amb%  94.432898 (95.773014)
precision       0.150192 (0.000000)
recall          0.112064 (0.000000)

bests[10].suffixonly == [false]
bests[10].langbase == [is]
comp = comp_parms0_off
bests[10].rows == [1]
  R->R     W->R     R->W     W->W

0.000000 0.704723 -0.705629 0.073845
*/
    };
#elif 1
static bestParms best_is =
    {
    false,
    "is",
    1,
//iteration:15.-1
/* number of nodes: 111593, nodes/line: 1.06693832594907057e-01 weight: 1.05561353388409800e+05 blobs 1 lines 5881633 * fraction 1.77827941003892431e-01 = 1045918 lines*/
        {                                                 // # decisions
        5.66585282075018903e-03,    5.64699254907257475e-01,    -8.21269355775872678e-01,    8.12360442321313353e-02, //4842840
        -3.53836789873742896e-01,    6.35502799370625149e-01,    3.77862909934449431e-01,    -5.72848443674078389e-01, //0
        2.85433068754089581e-01,    -3.75433738261612693e-01,    -3.36787965276463597e-01,    -8.14954807263245429e-01, //0
        -8.90671312833910434e-01,    -3.69188910551492833e-01,    -2.63268176968452838e-01,    -3.30760913081912244e-02  //0
        }                                                 //(0 unresolved comparisons)
/*
0    0.985600    2    2831993.000000    471048.500000    289516.500000    16.633110    10.223066    36063.500000    712.500000    528.500000    83.000000    3989.500000    36264.000000    362.000000    305.500000    0.000000    4445.500000    0.176020    0.158931    0.029052    0.006836    0.052977    0.276848    0.061522    0.029052    0.000000    0.244378    87.158325    1.721971    1.277280    0.200595    9.641830    3.905551    87.642893    0.874882    0.738333    0.000000    10.743892    2.341881    1.841603    33.482369    2.063949    62.612079    0.359125    0.058064    0.964304    34.168741    1.377577    63.489378    0.416667    0.038754
1    0.985600    2    2831993.000000    128789.000000    148030.500000    4.547645    5.227079    36520.000000    0.000000    0.000000    0.000000    4857.000000    36661.500000    382.500000    410.500000    0.000000    3922.500000    0.051268    0.000000    0.000000    0.000000    0.051268    0.220453    0.015380    0.073484    0.000000    0.278557    88.261595    0.000000    0.000000    0.000000    11.738405    0.000000    88.603572    0.924427    0.992097    0.000000    9.479904    2.309254    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.856756    34.093820    1.452498    63.596926    0.458779    0.040862
2    0.985600    2    2831993.000000    84003.000000    57418.500000    2.966215    2.027494    36089.500000    0.000000    0.000000    0.000000    5287.500000    36885.000000    307.500000    310.500000    0.000000    3874.000000    0.063231    0.000000    0.000000    0.000000    0.063231    0.140133    0.022216    0.059813    0.000000    0.177729    87.221162    0.000000    0.000000    0.000000    12.778838    0.000000    89.143727    0.743166    0.750417    0.000000    9.362689    1.737680    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.694830    34.503468    1.042850    63.758852    0.428713    0.029338
3    0.985600    2    2831993.000000    63704.000000    40149.000000    2.249441    1.417694    35649.000000    0.000000    0.000000    0.000000    5728.000000    36788.000000    286.000000    271.500000    0.000000    4031.500000    0.017089    0.000000    0.000000    0.000000    0.017089    0.256341    0.058104    0.008545    0.000000    0.206781    86.156560    0.000000    0.000000    0.000000    13.843440    0.000000    88.909297    0.691205    0.656162    0.000000    9.743336    1.590255    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.631994    34.588056    0.958262    63.821688    0.431213    0.026958
4    0.985600    2    2831993.000000    51623.000000    32297.000000    1.822851    1.140434    35334.000000    0.000000    0.000000    0.000000    6043.000000    36652.500000    283.000000    234.000000    0.000000    4207.500000    0.061522    0.000000    0.000000    0.000000    0.061522    0.251214    0.061522    0.006836    0.000000    0.182856    85.395268    0.000000    0.000000    0.000000    14.604732    0.000000    88.581821    0.683955    0.565532    0.000000    10.168693    1.494792    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.600575    34.652101    0.894217    63.853107    0.426759    0.025156
5    0.985600    2    2831993.000000    43987.500000    27785.500000    1.553235    0.981129    35034.500000    0.000000    0.000000    0.000000    6342.500000    36539.000000    272.500000    218.500000    0.000000    4347.000000    0.093992    0.000000    0.000000    0.000000    0.093992    0.215326    0.035888    0.005127    0.000000    0.174312    84.671436    0.000000    0.000000    0.000000    15.328564    0.000000    88.307514    0.658578    0.528071    0.000000    10.505837    1.436789    0.000000    35.546318    0.000000    64.453682    0.000000    0.000000    0.592116    34.701646    0.844672    63.861566    0.416319    0.023763

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      2
trainlines      2831993.000000
rules           57418.500000 (84003.000000)
rules%         2.027494 (2.966215)
same%stdev     0.140133
ambi1%stdev    0.022216
ambi2%stdev    0.059813
ambi3%stdev    0.000000
diff%stdev     0.177729
same%          89.143727 (87.221162)
ambi1%         0.743166 (0.000000)
ambi2%         0.750417 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          9.362689 (12.778838)
amb.rules%     1.737680 (0.000000)
false_amb%     0.694830 (0.000000)
false_not_amb% 34.503468 (35.546318)
true_amb%      1.042850 (0.000000)
true_not_amb%  63.758852 (64.453682)
precision       0.428713 (0.000000)
recall          0.029338 (0.000000)

bests[10].suffixonly == [false]
bests[10].langbase == [is]
comp = comp_parms0_off
bests[10].rows == [1]
  R->R     W->R      R->W     W->W
0.005666 0.564699 -0.821269 0.081236
*/
    };
#elif 1
static bestParms best_is =
    {
    false,
    "is",
    4,
//iteration:12.-1
/* number of nodes: 60151, nodes/line: 0.162086 weight: 57935.959065 blobs 1 lines 5881633 * fraction 0.063096 = 371105 lines*/
        {                                                         // # decisions
        0.013188,       0.712282,       -0.695044,      0.096930, //2283425
        -0.235599,      0.552619,       0.471857,       -0.645334, //43635
        0.170415,       -0.386280,      -0.498239,      -0.757301, //2353
        -0.956702,      -0.195077,      -0.214532,      0.025362  //0
        }                                                         //(0 unresolved comparisons)
// Same parameters since:
//iteration:1.41
/* number of nodes: 3149, nodes/line: 0.379032 weight: 3077.123342 blobs 1 lines 5881633 * fraction 0.001413 = 8308 lines*/
    /*
        {                                                         // # decisions
        0.013188,       0.712282,       -0.695044,      0.096930, //100450
        -0.235599,      0.552619,       0.471857,       -0.645334, //2918
        0.170415,       -0.386280,      -0.498239,      -0.757301, //124
        -0.956702,      -0.195077,      -0.214532,      0.025362  //0
        }                                                         //(0 unresolved comparisons)
    */

    };
#endif

static bestParms best_nlC0 = // Dutch, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*rcount*exp(-2.0*rcount/1.0)
// Aiming at cutoff == 0, because maximum penalty is for 1 word/lemma pairs for a rule.
    {
    false,
    "nlC0", 
    1,
/* Figures for SLOW 0
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0    0.985600    7    306848.000000    73132.000000    51919.000000    23.833299    16.920104    3635.857143    78.000000    117.000000    4.714286    648.428571    3797.285714    47.285714    28.857143    0.000000    610.571429    0.368578    0.192276    0.216604    0.027956    0.556177    0.613962    0.144366    0.137131    0.000000    0.612630    81.085128    1.739518    2.609277    0.105136    14.460940    4.982796    84.685230    1.054543    0.643558    0.000000    13.616669    2.182363    4.610042    5.428826    0.372754    89.588378    0.038858    0.064250    1.815981    5.435198    0.366382    92.382439    0.091633    0.063152
1    0.985600    7    306848.000000    18391.714286    24506.428571    5.993754    7.986504    3743.714286    0.000000    0.000000    0.000000    740.285714    3857.285714    44.714286    38.571429    0.000000    543.428571    0.341564    0.000000    0.000000    0.000000    0.341564    0.260942    0.167244    0.162285    0.000000    0.313508    83.490506    0.000000    0.000000    0.000000    16.509494    0.000000    86.023321    0.997196    0.860201    0.000000    12.119281    2.121830    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.535619    5.215369    0.586211    92.662801    0.160279    0.101043
2    0.985600    7    306848.000000    11355.142857    8512.857143    3.700576    2.774291    3696.857143    0.000000    0.000000    0.000000    787.142857    3886.285714    36.571429    31.000000    0.000000    530.142857    0.396739    0.000000    0.000000    0.000000    0.396739    0.408533    0.152814    0.077255    0.000000    0.357852    82.445521    0.000000    0.000000    0.000000    17.554479    0.000000    86.670065    0.815598    0.691347    0.000000    11.822990    1.710845    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.268001    5.358736    0.442844    92.930419    0.148663    0.076332
3    0.985600    7    306848.000000    8112.285714    5634.571429    2.643747    1.836274    3652.571429    0.000000    0.000000    0.000000    831.428571    3872.285714    33.000000    27.000000    0.000000    551.714286    0.255206    0.000000    0.000000    0.000000    0.255206    0.329961    0.112985    0.082445    0.000000    0.356293    81.457882    0.000000    0.000000    0.000000    18.542118    0.000000    86.357844    0.735950    0.602141    0.000000    12.304065    1.532433    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.156493    5.425640    0.375940    93.041927    0.139810    0.064800
4    0.985600    7    306848.000000    6239.285714    4347.571429    2.033347    1.416849    3615.000000    0.000000    0.000000    0.000000    869.000000    3850.571429    35.714286    22.571429    0.000000    575.142857    0.273440    0.000000    0.000000    0.000000    0.273440    0.311652    0.095241    0.028377    0.000000    0.290613    80.619982    0.000000    0.000000    0.000000    19.380018    0.000000    85.873582    0.796483    0.503377    0.000000    12.826558    1.491016    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.156493    5.467057    0.334523    93.041927    0.126354    0.057661
5    0.985600    7    306848.000000    5044.285714    3567.571429    1.643904    1.162651    3592.714286    0.000000    0.000000    0.000000    891.285714    3839.428571    28.428571    20.857143    0.000000    595.285714    0.262209    0.000000    0.000000    0.000000    0.262209    0.328738    0.117908    0.037381    0.000000    0.347817    80.122977    0.000000    0.000000    0.000000    19.877023    0.000000    85.625080    0.634000    0.465146    0.000000    13.275774    1.325347    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.022684    5.498917    0.302663    93.175736    0.128901    0.052169
cutoff 0    Affix a 0.534867 b 0.817364
            Suffix a -0.706998 b 0.933881
cutoff 1    Affix a 0.011389 b 0.800182
            Suffix a -2.482111 b 0.967669
cutoff 2    Affix a -1.095794 b 0.809484
            Suffix a -3.356905 b 0.996574
cutoff 3    Affix a -1.997422 b 0.849935
            Suffix a -3.794027 b 1.003571
cutoff 4    Affix a -2.581689 b 0.876345
            Suffix a -4.020226 b 1.000475
cutoff 5    Affix a -2.983694 b 0.892699
            Suffix a -4.205819 b 0.997875

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      7
trainlines      306848.000000
rules           8512.857143 (11355.142857)
rules%         2.774291 (3.700576)
same%stdev     0.408533
ambi1%stdev    0.152814
ambi2%stdev    0.077255
ambi3%stdev    0.000000
diff%stdev     0.357852
same%          86.670065 (82.445521)
ambi1%         0.815598 (0.000000)
ambi2%         0.691347 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          11.822990 (17.554479)
amb.rules%     1.710845 (0.000000)
false_amb%     1.268001 (0.000000)
false_not_amb% 5.358736 (5.801580)
true_amb%      0.442844 (0.000000)
true_not_amb%  92.930419 (94.198420)
precision       0.148663 (0.000000)
recall          0.076332 (0.000000)

bests[11].suffixonly == [false]
bests[11].langbase == [nlC0]
comp = comp_parms0_off
bests[11].rows == [1]
  R->R     W->R     R->W     W->W

0.280881 0.707997 -0.646540 0.002056 0.042723 0.001589 
*/
/* Figures for SLOW 1
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          7  306848.000000   73132.000000   45134.714286      23.833299      14.709144    3635.857143      78.000000     117.000000       4.714286     648.428571    3830.714286     108.142857     442.714286       0.000000     102.428571       0.368578       0.192276       0.216604       0.027956       0.556177       1.713095       1.518002       2.627062       0.000000       2.779675      81.085128       1.739518       2.609277       0.105136      14.460940       4.982796      85.430738       2.411750       9.873200       0.000000       2.284312      12.648146       4.610042       5.428826       0.372754      89.588378       0.038858       0.064250       8.780426       1.933860       3.867720      85.417994       0.180494       0.666667
1       0.985600          7  306848.000000   18391.714286   21055.714286       5.993754       6.861936    3743.714286       0.000000       0.000000       0.000000     740.285714    3788.857143      93.857143     342.571429       0.000000     258.714286       0.341564       0.000000       0.000000       0.000000       0.341564       0.875026       0.607954       1.583138       0.000000       1.711304      83.490506       0.000000       0.000000       0.000000      16.509494       0.000000      84.497260       2.093157       7.639862       0.000000       5.769721      10.198165       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       7.795973       3.399388       2.402192      86.402447       0.133499       0.414058
2       0.985600          7  306848.000000   11355.142857    8544.571429       3.700576       2.784627    3696.857143       0.000000       0.000000       0.000000     787.142857    3744.428571      67.571429     268.428571       0.000000     403.571429       0.396739       0.000000       0.000000       0.000000       0.396739       0.467952       0.257791       0.758127       0.000000       0.783612      82.445521       0.000000       0.000000       0.000000      17.554479       0.000000      83.506436       1.506945       5.986364       0.000000       9.000255       8.270677       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       7.114184       4.645087       1.156493      87.084236       0.075171       0.199341
3       0.985600          7  306848.000000    8112.285714    5744.714286       2.643747       1.872169    3652.571429       0.000000       0.000000       0.000000     831.428571    3696.571429      62.714286     263.000000       0.000000     461.714286       0.255206       0.000000       0.000000       0.000000       0.255206       0.324679       0.181961       0.539402       0.000000       0.600171      81.457882       0.000000       0.000000       0.000000      18.542118       0.000000      82.439149       1.398624       5.865299       0.000000      10.296929       8.066777       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       7.050465       4.785268       1.016312      87.147955       0.067229       0.175178
4       0.985600          7  306848.000000    6239.285714    4401.142857       2.033347       1.434307    3615.000000       0.000000       0.000000       0.000000     869.000000    3655.142857      64.428571     270.000000       0.000000     494.428571       0.273440       0.000000       0.000000       0.000000       0.273440       0.327619       0.260351       0.471155       0.000000       0.436422      80.619982       0.000000       0.000000       0.000000      19.380018       0.000000      81.515229       1.436855       6.021409       0.000000      11.026507       8.471390       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       7.423219       4.753409       1.048171      86.775201       0.065945       0.180670
5       0.985600          7  306848.000000    5044.285714    3640.000000       1.643904       1.186255    3592.714286       0.000000       0.000000       0.000000     891.285714    3620.428571      61.285714     278.714286       0.000000     523.571429       0.262209       0.000000       0.000000       0.000000       0.262209       0.257148       0.212241       0.442619       0.000000       0.381709      80.122977       0.000000       0.000000       0.000000      19.877023       0.000000      80.741048       1.366764       6.215751       0.000000      11.676437       8.640245       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       7.620747       4.782082       1.019498      86.577673       0.062696       0.175728
cutoff 0 Affix  a       0.798721 b       0.789038: N(rules)=      2.222696*N(trainpairs)^0.789038
         Suffix a      -0.706998 b       0.933881: N(rules)=      0.493123*N(trainpairs)^0.933881
cutoff 1 Affix  a       0.293562 b       0.770051: N(rules)=      1.341196*N(trainpairs)^0.770051
         Suffix a      -2.482111 b       0.967669: N(rules)=      0.083567*N(trainpairs)^0.967669
cutoff 2 Affix  a      -1.114487 b       0.811816: N(rules)=      0.328084*N(trainpairs)^0.811816
         Suffix a      -3.356905 b       0.996574: N(rules)=      0.034843*N(trainpairs)^0.996574
cutoff 3 Affix  a      -2.040671 b       0.854564: N(rules)=      0.129941*N(trainpairs)^0.854564
         Suffix a      -3.794027 b       1.003571: N(rules)=      0.022505*N(trainpairs)^1.003571
cutoff 4 Affix  a      -2.607640 b       0.879272: N(rules)=      0.073708*N(trainpairs)^0.879272
         Suffix a      -4.020226 b       1.000475: N(rules)=      0.017949*N(trainpairs)^1.000475
cutoff 5 Affix  a      -3.025073 b       0.897062: N(rules)=      0.048554*N(trainpairs)^0.897062
         Suffix a      -4.205819 b       0.997875: N(rules)=      0.014909*N(trainpairs)^0.997875

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
cutoff                      0
fraction          9856.000000
iterations                  7
trainlines      306848.000000
rules            45134.714286 (  73132.000000)
rules%              14.709144 (     23.833299)
same%stdev           1.713095
ambi1%stdev          1.518002
ambi2%stdev          2.627062
ambi3%stdev          0.000000
diff%stdev           2.779675
same%               85.430738 (     81.085128)
ambi1%               2.411750 (      1.739518)
ambi2%               9.873200 (      2.609277)
ambi3%               0.000000 (      0.105136)
diff%                2.284312 (     14.460940)
amb.rules%          12.648146 (      4.982796)
false_amb%           8.780426 (      4.610042)
false_not_amb%       1.933860 (      5.428826)
true_amb%            3.867720 (      0.372754)
true_not_amb%       85.417994 (     89.588378)
precision            0.180494 (      0.038858)
recall               0.666667 (      0.064250)

bests[11].suffixonly == [false]
bests[11].langbase == [nlC0]
comp = comp_parms0_off
bests[11].rows == [1]
  R->R     W->R     R->W     W->W

0.280881 0.707997 -0.646540 0.002056 0.042723 0.001589 
*/
//iteration:20.2
/*weight ( used): 1.44449673909528588e+03 suffix only: no */
/* number of nodes: 45281, nodes/line: 1.34238713139272547e-01 weight ( used): 1.44449673909528588e+03 blobs 1 lines 337317 * fraction 1.00000000000000000e+00 = 337317 lines*/
        {{ // These were computed using dict_nl_without_doubles_UTF8, a file with doublets when only looking a word and lemma
        2.80881366831663093e-01,    7.07996695016342659e-01,    -6.46540271422676893e-01,    2.05573578451381843e-03,    4.27231086744624400e-02,    1.58901761416129893e-03
        }}
    };

static bestParms best_nlC1 = // Dutch, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*rcount*exp(-2.0*rcount/2.0)
// Aiming at cutoff == 1, because maximum penalty is for 2 word/lemma pairs for a rule.
    {
    false,
    "nlC1", 
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0    0.985600    7    306848.000000    73132.000000    46062.571429    23.833299    15.011527    3635.857143    78.000000    117.000000    4.714286    648.428571    3755.857143    64.000000    40.142857    0.000000    624.000000    0.368578    0.192276    0.216604    0.027956    0.556177    0.444807    0.183452    0.135916    0.000000    0.461737    81.085128    1.739518    2.609277    0.105136    14.460940    4.982796    83.761310    1.427297    0.895247    0.000000    13.916146    2.905569    4.610042    5.428826    0.372754    89.588378    0.038858    0.064250    2.532815    5.428826    0.372754    91.665605    0.068541    0.064250
1    0.985600    7    306848.000000    18391.714286    22782.285714    5.993754    7.424616    3743.714286    0.000000    0.000000    0.000000    740.285714    3821.428571    69.285714    47.714286    0.000000    545.571429    0.341564    0.000000    0.000000    0.000000    0.341564    0.277609    0.142051    0.110548    0.000000    0.298889    83.490506    0.000000    0.000000    0.000000    16.509494    0.000000    85.223652    1.545177    1.064101    0.000000    12.167070    3.080795    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    2.312986    5.033771    0.767809    91.885434    0.142351    0.132345
2    0.985600    7    306848.000000    11355.142857    8143.714286    3.700576    2.653990    3696.857143    0.000000    0.000000    0.000000    787.142857    3857.285714    57.857143    34.285714    0.000000    534.571429    0.396739    0.000000    0.000000    0.000000    0.396739    0.298017    0.105729    0.157018    0.000000    0.306015    82.445521    0.000000    0.000000    0.000000    17.554479    0.000000    86.023321    1.290302    0.764623    0.000000    11.921754    2.363961    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.892443    5.330062    0.471518    92.305977    0.110778    0.081274
3    0.985600    7    306848.000000    8112.285714    5475.714286    2.643747    1.784504    3652.571429    0.000000    0.000000    0.000000    831.428571    3845.714286    50.571429    30.857143    0.000000    556.857143    0.255206    0.000000    0.000000    0.000000    0.255206    0.205091    0.174112    0.113299    0.000000    0.281926    81.457882    0.000000    0.000000    0.000000    18.542118    0.000000    85.765261    1.127820    0.688161    0.000000    12.418759    2.121830    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.701287    5.381037    0.420543    92.497133    0.110000    0.072488
4    0.985600    7    306848.000000    6239.285714    4244.428571    2.033347    1.383235    3615.000000    0.000000    0.000000    0.000000    869.000000    3826.142857    51.714286    28.285714    0.000000    577.857143    0.273440    0.000000    0.000000    0.000000    0.273440    0.224655    0.192584    0.151100    0.000000    0.308597    80.619982    0.000000    0.000000    0.000000    19.380018    0.000000    85.328788    1.153307    0.630814    0.000000    12.887091    2.096343    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.723589    5.428826    0.372754    92.474831    0.097581    0.064250
5    0.985600    7    306848.000000    5044.285714    3439.571429    1.643904    1.120937    3592.714286    0.000000    0.000000    0.000000    891.285714    3807.571429    49.571429    27.428571    0.000000    599.428571    0.262209    0.000000    0.000000    0.000000    0.262209    0.228471    0.106065    0.100445    0.000000    0.321600    80.122977    0.000000    0.000000    0.000000    19.877023    0.000000    84.914617    1.105518    0.611699    0.000000    13.368166    2.042182    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.720403    5.479801    0.321779    92.478017    0.085521    0.055464
cutoff 0    Affix a 0.384392 b 0.819628
            Suffix a -0.706998 b 0.933881
cutoff 1    Affix a -0.001749 b 0.794660
            Suffix a -2.482111 b 0.967669
cutoff 2    Affix a -1.188392 b 0.811560
            Suffix a -3.356905 b 0.996574
cutoff 3    Affix a -2.070862 b 0.850851
            Suffix a -3.794027 b 1.003571
cutoff 4    Affix a -2.610700 b 0.874822
            Suffix a -4.020226 b 1.000475
cutoff 5    Affix a -2.980595 b 0.887978
            Suffix a -4.205819 b 0.997875

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      7
trainlines      306848.000000
rules           8143.714286 (11355.142857)
rules%         2.653990 (3.700576)
same%stdev     0.298017
ambi1%stdev    0.105729
ambi2%stdev    0.157018
ambi3%stdev    0.000000
diff%stdev     0.306015
same%          86.023321 (82.445521)
ambi1%         1.290302 (0.000000)
ambi2%         0.764623 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          11.921754 (17.554479)
amb.rules%     2.363961 (0.000000)
false_amb%     1.892443 (0.000000)
false_not_amb% 5.330062 (5.801580)
true_amb%      0.471518 (0.000000)
true_not_amb%  92.305977 (94.198420)
precision       0.110778 (0.000000)
recall          0.081274 (0.000000)

bests[12].suffixonly == [false]
bests[12].langbase == [nlC1]
comp = comp_parms0_off
bests[12].rows == [1]
  R->R     W->R     R->W     W->W

0.069568 0.656948 -0.726409 0.186749 -0.004976 0.031780 
*/
//iteration:20.3
/*weight ( used): 7.72334598212476067e+03 suffix only: no */
/* number of nodes: 40017, nodes/line: 1.28534811712255725e-01 weight ( used): 7.72334598212476067e+03 blobs 1 lines 311332 * fraction 1.00000000000000000e+00 = 311332 lines*/
        {{
        6.95681858201713382e-02,    6.56948240807808115e-01,    -7.26408565910014570e-01,    1.86748935639343000e-01,    -4.97642864640904150e-03,    3.17795844729270860e-02
        }}

    };

static bestParms best_nlC2 = // Dutch, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*rcount*exp(-2.0*rcount/3.0)
// Aiming at cutoff == 2, because maximum penalty is for 3 word/lemma pairs for a rule.
    {
    false,
    "nlC2", 
    1,
/* SLOW 0
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0    0.985600    7    306848.000000    73132.000000    43253.857143    23.833299    14.096183    3635.857143    78.000000    117.000000    4.714286    648.428571    3766.428571    54.714286    38.142857    0.000000    624.714286    0.368578    0.192276    0.216604    0.027956    0.556177    0.521089    0.122634    0.117606    0.000000    0.496464    81.085128    1.739518    2.609277    0.105136    14.460940    4.982796    83.997069    1.220212    0.850644    0.000000    13.932076    2.606091    4.610042    5.428826    0.372754    89.588378    0.038858    0.064250    2.233338    5.428826    0.372754    91.965082    0.077024    0.064250
1    0.985600    7    306848.000000    18391.714286    21760.000000    5.993754    7.091459    3743.714286    0.000000    0.000000    0.000000    740.285714    3831.857143    55.428571    47.714286    0.000000    549.000000    0.341564    0.000000    0.000000    0.000000    0.341564    0.357620    0.186079    0.099498    0.000000    0.409401    83.490506    0.000000    0.000000    0.000000    16.509494    0.000000    85.456225    1.236141    1.064101    0.000000    12.243533    2.669810    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.902001    5.033771    0.767809    92.296419    0.167944    0.132345
2    0.985600    7    306848.000000    11355.142857    8030.142857    3.700576    2.616977    3696.857143    0.000000    0.000000    0.000000    787.142857    3868.428571    43.142857    32.285714    0.000000    540.142857    0.396739    0.000000    0.000000    0.000000    0.396739    0.345939    0.117606    0.131758    0.000000    0.327113    82.445521    0.000000    0.000000    0.000000    17.554479    0.000000    86.271824    0.962151    0.720020    0.000000    12.046005    1.876513    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.465528    5.390595    0.410985    92.732892    0.122974    0.070840
3    0.985600    7    306848.000000    8112.285714    5347.000000    2.643747    1.742557    3652.571429    0.000000    0.000000    0.000000    831.428571    3857.571429    37.428571    28.571429    0.000000    560.428571    0.255206    0.000000    0.000000    0.000000    0.255206    0.209546    0.163810    0.088138    0.000000    0.257791    81.457882    0.000000    0.000000    0.000000    18.542118    0.000000    86.029693    0.834714    0.637186    0.000000    12.498407    1.663056    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.293488    5.432012    0.369568    92.904932    0.125000    0.063701
4    0.985600    7    306848.000000    6239.285714    4134.000000    2.033347    1.347247    3615.000000    0.000000    0.000000    0.000000    869.000000    3836.285714    36.142857    26.714286    0.000000    584.857143    0.273440    0.000000    0.000000    0.000000    0.273440    0.171646    0.125780    0.106733    0.000000    0.262435    80.619982    0.000000    0.000000    0.000000    19.380018    0.000000    85.554989    0.806041    0.595769    0.000000    13.043201    1.615267    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.277558    5.463872    0.337709    92.920861    0.116740    0.058210
5    0.985600    7    306848.000000    5044.285714    3407.285714    1.643904    1.110415    3592.714286    0.000000    0.000000    0.000000    891.285714    3818.142857    33.714286    24.428571    0.000000    607.714286    0.262209    0.000000    0.000000    0.000000    0.262209    0.201303    0.120589    0.088138    0.000000    0.283769    80.122977    0.000000    0.000000    0.000000    19.877023    0.000000    85.150376    0.751880    0.544794    0.000000    13.552950    1.516503    0.000000    5.801580    0.000000    94.198420    0.000000    0.000000    1.217026    5.502103    0.299478    92.981394    0.109557    0.051620
cutoff 0    Affix a 0.416960 b 0.814680
            Suffix a -0.706998 b 0.933881
cutoff 1    Affix a 0.031775 b 0.790236
            Suffix a -2.482111 b 0.967669
cutoff 2    Affix a -1.153585 b 0.808514
            Suffix a -3.356905 b 0.996574
cutoff 3    Affix a -2.028061 b 0.846547
            Suffix a -3.794027 b 1.003571
cutoff 4    Affix a -2.571438 b 0.870756
            Suffix a -4.020226 b 1.000475
cutoff 5    Affix a -2.961484 b 0.885966
            Suffix a -4.205819 b 0.997875

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      7
trainlines      306848.000000
rules           8030.142857 (11355.142857)
rules%         2.616977 (3.700576)
same%stdev     0.345939
ambi1%stdev    0.117606
ambi2%stdev    0.131758
ambi3%stdev    0.000000
diff%stdev     0.327113
same%          86.271824 (82.445521)
ambi1%         0.962151 (0.000000)
ambi2%         0.720020 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.046005 (17.554479)
amb.rules%     1.876513 (0.000000)
false_amb%     1.465528 (0.000000)
false_not_amb% 5.390595 (5.801580)
true_amb%      0.410985 (0.000000)
true_not_amb%  92.732892 (94.198420)
precision       0.122974 (0.000000)
recall          0.070840 (0.000000)

bests[13].suffixonly == [false]
bests[13].langbase == [nlC2]
comp = comp_parms0_off
bests[13].rows == [1]
  R->R     W->R     R->W     W->W

0.052945 0.603264 -0.766378 0.211174 0.017177 0.032318 
*/
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          7  306848.000000   73132.000000   38529.714286      23.833299      12.556612    3635.857143      78.000000     117.000000       4.714286     648.428571    3281.428571       0.000000       0.000000       0.000000    1202.571429       0.368578       0.192276       0.216604       0.027956       0.556177       0.528513       0.000000       0.000000       0.000000       0.528513      81.085128       1.739518       2.609277       0.105136      14.460940       4.982796      73.180833       0.000000       0.000000       0.000000      26.819167       0.000000       4.610042       5.428826       0.372754      89.588378       0.038858       0.064250       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000
1       0.985600          7  306848.000000   18391.714286   18986.142857       5.993754       6.187475    3743.714286       0.000000       0.000000       0.000000     740.285714    3283.428571       0.000000       0.000000       0.000000    1200.571429       0.341564       0.000000       0.000000       0.000000       0.341564       0.617750       0.000000       0.000000       0.000000       0.617750      83.490506       0.000000       0.000000       0.000000      16.509494       0.000000      73.225436       0.000000       0.000000       0.000000      26.774564       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000
2       0.985600          7  306848.000000   11355.142857    7903.142857       3.700576       2.575589    3696.857143       0.000000       0.000000       0.000000     787.142857    3218.000000       0.000000       0.000000       0.000000    1266.000000       0.396739       0.000000       0.000000       0.000000       0.396739       0.605573       0.000000       0.000000       0.000000       0.605573      82.445521       0.000000       0.000000       0.000000      17.554479       0.000000      71.766280       0.000000       0.000000       0.000000      28.233720       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000
3       0.985600          7  306848.000000    8112.285714    5368.714286       2.643747       1.749633    3652.571429       0.000000       0.000000       0.000000     831.428571    3149.428571       0.000000       0.000000       0.000000    1334.571429       0.255206       0.000000       0.000000       0.000000       0.255206       0.594222       0.000000       0.000000       0.000000       0.594222      81.457882       0.000000       0.000000       0.000000      18.542118       0.000000      70.237033       0.000000       0.000000       0.000000      29.762967       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000
4       0.985600          7  306848.000000    6239.285714    4172.571429       2.033347       1.359817    3615.000000       0.000000       0.000000       0.000000     869.000000    3087.428571       0.000000       0.000000       0.000000    1396.571429       0.273440       0.000000       0.000000       0.000000       0.273440       0.489086       0.000000       0.000000       0.000000       0.489086      80.619982       0.000000       0.000000       0.000000      19.380018       0.000000      68.854339       0.000000       0.000000       0.000000      31.145661       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000
5       0.985600          7  306848.000000    5044.285714    3460.571429       1.643904       1.127780    3592.714286       0.000000       0.000000       0.000000     891.285714    3048.142857       0.000000       0.000000       0.000000    1435.857143       0.262209       0.000000       0.000000       0.000000       0.262209       0.504696       0.000000       0.000000       0.000000       0.504696      80.122977       0.000000       0.000000       0.000000      19.877023       0.000000      67.978208       0.000000       0.000000       0.000000      32.021792       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000       0.000000       5.801580       0.000000      94.198420       0.000000       0.000000
cutoff 0 Affix  a       0.661338 b       0.786985: N(rules)=      1.937383*N(trainpairs)^0.786985
         Suffix a      -0.706998 b       0.933881: N(rules)=      0.493123*N(trainpairs)^0.933881
cutoff 1 Affix  a       0.321587 b       0.758378: N(rules)=      1.379315*N(trainpairs)^0.758378
         Suffix a      -2.482111 b       0.967669: N(rules)=      0.083567*N(trainpairs)^0.967669
cutoff 2 Affix  a      -1.091314 b       0.802674: N(rules)=      0.335775*N(trainpairs)^0.802674
         Suffix a      -3.356905 b       0.996574: N(rules)=      0.034843*N(trainpairs)^0.996574
cutoff 3 Affix  a      -2.002932 b       0.844252: N(rules)=      0.134939*N(trainpairs)^0.844252
         Suffix a      -3.794027 b       1.003571: N(rules)=      0.022505*N(trainpairs)^1.003571
cutoff 4 Affix  a      -2.546591 b       0.868571: N(rules)=      0.078348*N(trainpairs)^0.868571
         Suffix a      -4.020226 b       1.000475: N(rules)=      0.017949*N(trainpairs)^1.000475
cutoff 5 Affix  a      -2.964140 b       0.886552: N(rules)=      0.051605*N(trainpairs)^0.886552
         Suffix a      -4.205819 b       0.997875: N(rules)=      0.014909*N(trainpairs)^0.997875

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        yes                                                     (SLOW 1)
cutoff                      1
fraction          9856.000000
iterations                  7
trainlines      306848.000000
rules            18986.142857 (  18391.714286)
rules%               6.187475 (      5.993754)
same%stdev           0.617750
ambi1%stdev          0.000000
ambi2%stdev          0.000000
ambi3%stdev          0.000000
diff%stdev           0.617750
same%               73.225436 (     83.490506)
ambi1%               0.000000 (      0.000000)
ambi2%               0.000000 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               26.774564 (     16.509494)
amb.rules%           0.000000 (      0.000000)
false_amb%           0.000000 (      0.000000)
false_not_amb%       5.801580 (      5.801580)
true_amb%            0.000000 (      0.000000)
true_not_amb%       94.198420 (     94.198420)
precision            0.000000 (      0.000000)
recall               0.000000 (      0.000000)

bests[13].suffixonly == [false]
bests[13].langbase == [nlC2]
comp = comp_parms0_off
bests[13].rows == [1]
  R->R     W->R     R->W     W->W

0.052945 0.603264 -0.766378 0.211174 0.017177 0.032318 
*/
//iteration:20.1
/*weight ( used): 1.74855637296820641e+04 suffix only: no */
/* number of nodes: 37996, nodes/line: 1.22043349222052344e-01 weight ( used): 1.74855637296820641e+04 blobs 1 lines 311332 * fraction 1.00000000000000000e+00 = 311332 lines*/
        {{
        5.29451965787332626e-02,    6.03263920561431743e-01,    -7.66378180909070839e-01,    2.11173926307903825e-01,    1.71771467887680100e-02,    3.23179611128934577e-02
        }}
    };



static bestParms best_daC0 =
    {
    false,
    "daC0", 
    1,
//iteration:20.3
/*weight ( used): 1.75942185857700269e+03 suffix only: no */
/* number of nodes: 59989, nodes/line: inf weight ( used): 1.75942185857700269e+03 blobs 1 lines 0 * fraction 1.00000000000000000e+00 = 0 lines*/
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          5  555065.000000   81861.400000   65106.200000      14.748075      11.729473    6963.600000      61.000000     134.400000       6.400000     944.400000    7061.200000      39.200000      28.400000       0.000000     981.000000       6.489145       0.130178       0.211105       0.014057       6.794025       6.658768       0.038417       0.018684       0.000000       6.656800      85.866483       0.752176       1.657254       0.078917      11.645170       2.924856      87.069965       0.483366       0.350194       0.000000      12.096476       1.141828       2.525340       7.714124       0.399517      89.361020       0.073303       0.049240       0.855755       7.827567       0.286074      91.030605       0.143210       0.035258
1       0.985600          5  555065.000000   25122.800000   32643.000000       4.526101       5.880933    7021.000000       0.000000       0.000000       0.000000    1088.800000    7092.400000      38.000000      39.800000       0.000000     939.600000       6.473833       0.000000       0.000000       0.000000       6.473833       6.727077       0.077506       0.052589       0.000000       6.732606      86.574268       0.000000       0.000000       0.000000      13.425732       0.000000      87.454684       0.468569       0.490764       0.000000      11.585982       1.181287       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.791635       7.723988       0.389652      91.094725       0.197500       0.048024
2       0.985600          5  555065.000000   16595.200000   13356.000000       2.989776       2.406205    6952.200000       0.000000       0.000000       0.000000    1157.600000    7105.800000      37.600000      36.400000       0.000000     930.000000       6.380172       0.000000       0.000000       0.000000       6.380172       6.751026       0.046619       0.109794       0.000000       6.819038      85.725912       0.000000       0.000000       0.000000      14.274088       0.000000      87.619917       0.463637       0.448840       0.000000      11.467607       1.053047       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.695455       7.756048       0.357592      91.190905       0.204513       0.044073
3       0.985600          5  555065.000000   12638.600000    9456.800000       2.276959       1.703728    6886.600000       0.000000       0.000000       0.000000    1223.200000    7076.000000      36.400000      30.400000       0.000000     967.000000       6.294505       0.000000       0.000000       0.000000       6.294505       6.669820       0.057565       0.117165       0.000000       6.738610      84.917014       0.000000       0.000000       0.000000      15.082986       0.000000      87.252460       0.448840       0.374855       0.000000      11.923845       0.956867       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.638733       7.795507       0.318134      91.247626       0.199382       0.039210
4       0.985600          5  555065.000000   10397.800000    7534.400000       1.873258       1.357391    6834.600000       0.000000       0.000000       0.000000    1275.200000    7049.200000      34.600000      26.800000       0.000000     999.200000       6.317295       0.000000       0.000000       0.000000       6.317295       6.733371       0.056237       0.063712       0.000000       6.736659      84.275814       0.000000       0.000000       0.000000      15.724186       0.000000      86.921996       0.426644       0.330464       0.000000      12.320896       0.892747       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.594343       7.815236       0.298404      91.292017       0.200663       0.036778
5       0.985600          5  555065.000000    8830.000000    6346.200000       1.590805       1.143326    6782.200000       0.000000       0.000000       0.000000    1327.600000    7010.000000      31.800000      26.200000       0.000000    1041.800000       6.337696       0.000000       0.000000       0.000000       6.337696       6.787495       0.081038       0.057433       0.000000       6.836868      83.629683       0.000000       0.000000       0.000000      16.370317       0.000000      86.438630       0.392118       0.323066       0.000000      12.846186       0.853289       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.586944       7.847296       0.266344      91.299416       0.184932       0.032827
cutoff 0 Affix  a       1.073154 b       0.757703: N(rules)=      2.924588*N(trainpairs)^0.757703
         Suffix a      -0.095489 b       0.854699: N(rules)=      0.908928*N(trainpairs)^0.854699
cutoff 1 Affix  a      -0.254162 b       0.813340: N(rules)=      0.775566*N(trainpairs)^0.813340
         Suffix a      -1.274755 b       0.858995: N(rules)=      0.279499*N(trainpairs)^0.858995
cutoff 2 Affix  a      -1.775778 b       0.868243: N(rules)=      0.169352*N(trainpairs)^0.868243
         Suffix a      -2.211908 b       0.899654: N(rules)=      0.109492*N(trainpairs)^0.899654
cutoff 3 Affix  a      -2.803938 b       0.922054: N(rules)=      0.060571*N(trainpairs)^0.922054
         Suffix a      -2.901611 b       0.932627: N(rules)=      0.054935*N(trainpairs)^0.932627
cutoff 4 Affix  a      -3.366096 b       0.946465: N(rules)=      0.034524*N(trainpairs)^0.946465
         Suffix a      -3.423357 b       0.957609: N(rules)=      0.032603*N(trainpairs)^0.957609
cutoff 5 Affix  a      -3.846574 b       0.971006: N(rules)=      0.021353*N(trainpairs)^0.971006
         Suffix a      -3.811529 b       0.975132: N(rules)=      0.022114*N(trainpairs)^0.975132

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        no
cutoff                      2
fraction          9856.000000
iterations                  5
trainlines      555065.000000
rules            13356.000000 (  16595.200000)
rules%               2.406205 (      2.989776)
same%stdev           6.751026
ambi1%stdev          0.046619
ambi2%stdev          0.109794
ambi3%stdev          0.000000
diff%stdev           6.819038
same%               87.619917 (     85.725912)
ambi1%               0.463637 (      0.000000)
ambi2%               0.448840 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               11.467607 (     14.274088)
amb.rules%           1.053047 (      0.000000)
false_amb%           0.695455 (      0.000000)
false_not_amb%       7.756048 (      8.113640)
true_amb%            0.357592 (      0.000000)
true_not_amb%       91.190905 (     91.886360)
precision            0.204513 (      0.000000)
recall               0.044073 (      0.000000)

bests[1].suffixonly == [false]
bests[1].langbase == [daC0]
comp = comp_parms0_off
bests[1].rows == [1]
  R->R     W->R     R->W     W->W

0.337801 0.729767 -0.592110 0.028323 0.043129 0.008601 
*/
        {{
        3.37800842415481084e-01,    7.29766766436709680e-01,    -5.92110476598071700e-01,    2.83227179535092063e-02,    4.31287955142629492e-02,    8.60067531762726684e-03
        }}
    };

static bestParms best_daC1 =
    {
    false,
    "daC1",
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          5  555065.000000   81861.400000   60375.000000      14.748075      10.877104    6963.600000      61.000000     134.400000       6.400000     944.400000    7025.800000      38.400000      24.400000       0.000000    1021.200000       6.489145       0.130178       0.211105       0.014057       6.794025       6.654557       0.059535       0.048222       0.000000       6.654642      85.866483       0.752176       1.657254       0.078917      11.645170       2.924856      86.633456       0.473501       0.300871       0.000000      12.592172       1.107302       2.525340       7.714124       0.399517      89.361020       0.073303       0.049240       0.816296       7.822634       0.291006      91.070063       0.151282       0.035866
1       0.985600          5  555065.000000   25122.800000   32443.200000       4.526101       5.844937    7021.000000       0.000000       0.000000       0.000000    1088.800000    7072.000000      33.000000      31.000000       0.000000     973.800000       6.473833       0.000000       0.000000       0.000000       6.473833       6.710209       0.048542       0.046133       0.000000       6.724254      86.574268       0.000000       0.000000       0.000000      13.425732       0.000000      87.203137       0.406915       0.382254       0.000000      12.007694       1.011122       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.685590       7.788108       0.325532      91.200769       0.191860       0.040122
2       0.985600          5  555065.000000   16595.200000   14189.000000       2.989776       2.556277    6952.200000       0.000000       0.000000       0.000000    1157.600000    7093.000000      33.200000      24.800000       0.000000     958.800000       6.380172       0.000000       0.000000       0.000000       6.380172       6.759828       0.066623       0.061275       0.000000       6.808558      85.725912       0.000000       0.000000       0.000000      14.274088       0.000000      87.462083       0.409381       0.305803       0.000000      11.822733       0.865619       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.601741       7.849762       0.263878      91.284619       0.179832       0.032523
3       0.985600          5  555065.000000   12638.600000   10146.400000       2.276959       1.827966    6886.600000       0.000000       0.000000       0.000000    1223.200000    7067.200000      27.200000      19.400000       0.000000     996.000000       6.294505       0.000000       0.000000       0.000000       6.294505       6.663662       0.048064       0.039564       0.000000       6.730579      84.917014       0.000000       0.000000       0.000000      15.082986       0.000000      87.143949       0.335397       0.239217       0.000000      12.281437       0.727515       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.532689       7.918814       0.194826      91.353671       0.154599       0.024012
4       0.985600          5  555065.000000   10397.800000    8153.400000       1.873258       1.468909    6834.600000       0.000000       0.000000       0.000000    1275.200000    7042.600000      23.000000      17.800000       0.000000    1026.400000       6.317295       0.000000       0.000000       0.000000       6.317295       6.694804       0.048550       0.048851       0.000000       6.721411      84.275814       0.000000       0.000000       0.000000      15.724186       0.000000      86.840613       0.283607       0.219488       0.000000      12.656292       0.675726       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.495697       7.933611       0.180029      91.390663       0.153684       0.022188
5       0.985600          5  555065.000000    8830.000000    6905.400000       1.590805       1.244071    6782.200000       0.000000       0.000000       0.000000    1327.600000    7003.200000      20.600000      15.400000       0.000000    1070.600000       6.337696       0.000000       0.000000       0.000000       6.337696       6.733081       0.043238       0.020628       0.000000       6.768454      83.629683       0.000000       0.000000       0.000000      16.370317       0.000000      86.354781       0.254014       0.189894       0.000000      13.201312       0.611606       0.000000       8.113640       0.000000      91.886360       0.000000       0.000000       0.463637       7.965671       0.147969      91.422723       0.137615       0.018237
cutoff 0 Affix  a       1.036474 b       0.756051: N(rules)=      2.819260*N(trainpairs)^0.756051
         Suffix a      -0.095489 b       0.854699: N(rules)=      0.908928*N(trainpairs)^0.854699
cutoff 1 Affix  a       0.279046 b       0.768751: N(rules)=      1.321868*N(trainpairs)^0.768751
         Suffix a      -1.274755 b       0.858995: N(rules)=      0.279499*N(trainpairs)^0.858995
cutoff 2 Affix  a      -1.140317 b       0.820294: N(rules)=      0.319718*N(trainpairs)^0.820294
         Suffix a      -2.211908 b       0.899654: N(rules)=      0.109492*N(trainpairs)^0.899654
cutoff 3 Affix  a      -2.173315 b       0.875709: N(rules)=      0.113800*N(trainpairs)^0.875709
         Suffix a      -2.901611 b       0.932627: N(rules)=      0.054935*N(trainpairs)^0.932627
cutoff 4 Affix  a      -2.860454 b       0.911966: N(rules)=      0.057243*N(trainpairs)^0.911966
         Suffix a      -3.423357 b       0.957609: N(rules)=      0.032603*N(trainpairs)^0.957609
cutoff 5 Affix  a      -3.347758 b       0.936701: N(rules)=      0.035163*N(trainpairs)^0.936701
         Suffix a      -3.811529 b       0.975132: N(rules)=      0.022114*N(trainpairs)^0.975132

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        no
cutoff                      2
fraction          9856.000000
iterations                  5
trainlines      555065.000000
rules            14189.000000 (  16595.200000)
rules%               2.556277 (      2.989776)
same%stdev           6.759828
ambi1%stdev          0.066623
ambi2%stdev          0.061275
ambi3%stdev          0.000000
diff%stdev           6.808558
same%               87.462083 (     85.725912)
ambi1%               0.409381 (      0.000000)
ambi2%               0.305803 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               11.822733 (     14.274088)
amb.rules%           0.865619 (      0.000000)
false_amb%           0.601741 (      0.000000)
false_not_amb%       7.849762 (      8.113640)
true_amb%            0.263878 (      0.000000)
true_not_amb%       91.284619 (     91.886360)
precision            0.179832 (      0.000000)
recall               0.032523 (      0.000000)

bests[2].suffixonly == [false]
bests[2].langbase == [daC1]
comp = comp_parms0_off
bests[2].rows == [1]
  R->R     W->R     R->W     W->W

0.138401 0.770176 -0.594523 0.183745 -0.001895 0.021206 
*/
//iteration:20.1
/*weight ( used): 1.10455190216752599e+04 suffix only: no */
/* number of nodes: 59917, nodes/line: 1.00941233155262197e-01 weight ( used): 1.10455190216752599e+04 blobs 1 lines 593583 * fraction 1.00000000000000000e+00 = 593583 lines*/
        {{
        1.38401137522408346e-01,    7.70176363116090945e-01,    -5.94523323442289975e-01,    1.83744971823973091e-01,    -1.89503087444323529e-03,    2.12062938841620120e-02
        }}
    };


static bestParms best_da3 =
    {
    false,
    "da3",
    1,
//iteration:20.5
/*weight (not  used): 5.00470146882884801e+04 suffix only: no */
/* number of nodes: 54059, nodes/line: 9.10723521394649072e-02 weight (not  used): 5.00470146882884801e+04 blobs 1 lines 593583 * fraction 1.00000000000000000e+00 = 593583 lines*/
        {{ // LOWEST fraction 0.001
        0.00000000000000000e+00,        7.35136672243838163e-01,        -6.17008258102692664e-01,       2.80846724309429585e-01
        }}
    };

static bestParms best_en4 = // English, ambiguous training pairs in training set derived from CELEX
    {
    false,
    "en4",
    1,
//iteration:20.5
/*weight (not  used): 1.02754092132184778e+04 suffix only: no */
/* number of nodes: 10801, nodes/line: 1.37949091280636538e-01 weight (not  used): 1.02754092132184778e+04 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        1.06322242479528462e-03,        7.78077313282915517e-01,        -5.98834680747542203e-01,       1.89714494033809300e-01
        }}
    };

static bestParms best_en4_suffix = // English, ambiguous training pairs in training set derived from CELEX
    {
    true,
    "en4",
    1,
//iteration:20.5
/*weight (not  used): 1.23994408913745046e+04 suffix only: yes */
/* number of nodes: 13281, nodes/line: 1.69623357216751591e-01 weight (not  used): 1.23994408913745046e+04 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        7.47138752010872206e-03,        7.20070629992721201e-01,        -6.91457516803205441e-01,       5.76972152426369414e-02
        }}
    };

static bestParms best_en4W = // English, ambiguous training pairs in training set derived from CELEX
    {
    false,
    "en4W",
    1,
/*
0    0.985600    14    77169.000000    14451.285714    11095.214286    18.726802    14.377813    935.285714    5.428571    15.500000    0.000000    171.785714    934.928571    14.428571    7.500000    0.000000    171.142857    1.122684    0.142160    0.476458    0.000000    1.182352    1.080444    0.338409    0.207908    0.000000    0.959778    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    82.883739    1.279129    0.664894    0.000000    15.172239    2.121327    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    2.070669    4.939210    0.050659    92.939463    0.012085    0.010152
1    0.985600    14    77169.000000    2563.357143    5913.000000    3.321745    7.662403    975.428571    0.000000    0.000000    0.000000    152.571429    955.357143    14.357143    10.428571    0.000000    147.857143    1.082261    0.000000    0.000000    0.000000    1.082261    1.168686    0.270705    0.213544    0.000000    1.031431    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    84.694782    1.272796    0.924519    0.000000    13.107903    2.368288    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.652736    4.274316    0.715552    93.357396    0.177953    0.143401
2    0.985600    14    77169.000000    1454.357143    1465.785714    1.884639    1.899449    971.785714    0.000000    0.000000    0.000000    156.214286    979.500000    6.285714    3.785714    0.000000    138.428571    1.146528    0.000000    0.000000    0.000000    1.146528    0.938049    0.244469    0.130933    0.000000    0.918487    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.835106    0.557244    0.335613    0.000000    12.272036    1.006839    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.797872    4.780902    0.208967    94.212259    0.115789    0.041878
3    0.985600    14    77169.000000    1004.857143    766.642857    1.302151    0.993460    965.928571    0.000000    0.000000    0.000000    162.071429    980.000000    4.714286    3.357143    0.000000    139.928571    1.079324    0.000000    0.000000    0.000000    1.079324    0.901403    0.194938    0.195491    0.000000    0.937864    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    86.879433    0.417933    0.297619    0.000000    12.405015    0.861196    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.690223    4.818896    0.170973    94.319909    0.110204    0.034264
4    0.985600    14    77169.000000    782.785714    538.428571    1.014378    0.697727    962.642857    0.000000    0.000000    0.000000    165.357143    979.357143    3.428571    2.857143    0.000000    142.357143    0.977102    0.000000    0.000000    0.000000    0.977102    0.893052    0.146351    0.134030    0.000000    0.896431    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.822442    0.303951    0.253293    0.000000    12.620314    0.677558    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.557244    4.869554    0.120314    94.452888    0.097436    0.024112
5    0.985600    14    77169.000000    653.928571    398.500000    0.847398    0.516399    960.142857    0.000000    0.000000    0.000000    167.857143    978.285714    2.214286    2.285714    0.000000    145.214286    1.007113    0.000000    0.000000    0.000000    1.007113    0.798576    0.126232    0.073176    0.000000    0.758786    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.727457    0.196302    0.202634    0.000000    12.873607    0.481256    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.367275    4.875887    0.113982    94.642857    0.134328    0.022843
cutoff 0    Affix a -0.766864 b 0.896467
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.202355 b 0.878363
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -2.151304 b 0.835420
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -2.479762 b 0.808472
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -2.713938 b 0.798032
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -2.807226 b 0.783347
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1465.785714 (1454.357143)
rules%         1.899449 (1.884639)
same%stdev     0.938049
ambi1%stdev    0.244469
ambi2%stdev    0.130933
ambi3%stdev    0.000000
diff%stdev     0.918487
same%          86.835106 (86.151216)
ambi1%         0.557244 (0.000000)
ambi2%         0.335613 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.272036 (13.848784)
amb.rules%     1.006839 (0.000000)
false_amb%     0.797872 (0.000000)
false_not_amb% 4.780902 (4.989868)
true_amb%      0.208967 (0.000000)
true_not_amb%  94.212259 (95.010132)
precision       0.115789 (0.000000)
recall          0.041878 (0.000000)

bests[5].suffixonly == [false]
bests[5].langbase == [en4W]
comp = comp_parms0_off
bests[5].rows == [1]
  R->R     W->R     R->W     W->W

0.004029 0.680122 -0.724874 0.109431 

*/
//iteration:20.9
/*weight ( used): 1.13023631391451454e+04 suffix only: no */
/* number of nodes: 9850, nodes/line: 1.25803032044650487e-01 weight ( used): 1.13023631391451454e+04 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        4.02884374087577289e-03,        6.80122085729639836e-01,        -7.24874245486841984e-01,       1.09430549440079072e-01
        }}
    };


/*
To compare: comp_sugar, with ambiguous training set instead of disambiguated set from ACL:

0    0.985600    14    77169.000000    14451.285714    18381.357143    18.726802    23.819613    935.285714    5.428571    15.500000    0.000000    171.785714    943.428571    10.071429    4.285714    0.000000    170.214286    1.122684    0.142160    0.476458    0.000000    1.182352    0.951464    0.215657    0.106772    0.000000    0.918793    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    83.637285    0.892857    0.379939    0.000000    15.089919    1.443769    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    1.405775    4.951874    0.037994    93.604357    0.013333    0.007614
1    0.985600    14    77169.000000    2563.357143    7041.000000    3.321745    9.124130    975.428571    0.000000    0.000000    0.000000    152.571429    966.928571    7.500000    6.785714    0.000000    146.785714    1.082261    0.000000    0.000000    0.000000    1.082261    0.890922    0.284093    0.241358    0.000000    0.767500    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    85.720618    0.664894    0.601570    0.000000    13.012918    1.494428    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.038501    4.533941    0.455927    93.971631    0.180000    0.091371
2    0.985600    14    77169.000000    1454.357143    1365.714286    1.884639    1.769771    971.785714    0.000000    0.000000    0.000000    156.214286    980.357143    3.642857    4.428571    0.000000    139.571429    1.146528    0.000000    0.000000    0.000000    1.146528    0.871120    0.204558    0.165723    0.000000    0.792714    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.911094    0.322948    0.392604    0.000000    12.373354    0.861196    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.658561    4.787234    0.202634    94.351570    0.133333    0.040609
3    0.985600    14    77169.000000    1004.857143    716.642857    1.302151    0.928667    965.928571    0.000000    0.000000    0.000000    162.071429    979.714286    2.071429    3.285714    0.000000    142.928571    1.079324    0.000000    0.000000    0.000000    1.079324    0.778644    0.101165    0.122587    0.000000    0.766937    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    86.854103    0.183637    0.291287    0.000000    12.670973    0.639564    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.500253    4.850557    0.139311    94.509878    0.122222    0.027919
4    0.985600    14    77169.000000    782.785714    455.928571    1.014378    0.590818    962.642857    0.000000    0.000000    0.000000    165.357143    977.857143    2.142857    2.571429    0.000000    145.428571    0.977102    0.000000    0.000000    0.000000    0.977102    0.701383    0.119736    0.128772    0.000000    0.692585    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.689463    0.189970    0.227964    0.000000    12.892604    0.538247    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.436930    4.888551    0.101317    94.573202    0.103896    0.020305
5    0.985600    14    77169.000000    653.928571    342.071429    0.847398    0.443276    960.142857    0.000000    0.000000    0.000000    167.857143    976.500000    1.571429    2.357143    0.000000    147.571429    1.007113    0.000000    0.000000    0.000000    1.007113    0.704730    0.102648    0.118467    0.000000    0.708976    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.569149    0.139311    0.208967    0.000000    13.082573    0.455927    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.373607    4.907548    0.082320    94.636525    0.099237    0.016497
cutoff 0    Affix a -0.301365 b 0.902523
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.463083 b 0.912868
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -1.744098 b 0.787651
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -1.561877 b 0.715037
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -1.408205 b 0.667043
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -1.279623 b 0.632865  
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1365.714286 (1454.357143)
rules%         1.769771 (1.884639)
same%stdev     0.871120
ambi1%stdev    0.204558
ambi2%stdev    0.165723
ambi3%stdev    0.000000
diff%stdev     0.792714
same%          86.911094 (86.151216)
ambi1%         0.322948 (0.000000)
ambi2%         0.392604 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.373354 (13.848784)
amb.rules%     0.861196 (0.000000)
false_amb%     0.658561 (0.000000)
false_not_amb% 4.787234 (4.989868)
true_amb%      0.202634 (0.000000)
true_not_amb%  94.351570 (95.010132)
precision       0.133333 (0.000000)
recall          0.040609 (0.000000)
*/

static bestParms best_en6W = // English, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) 5.0*(rcount-1.0)*exp(-0.9*rcount)+1
    {
    false,
    "en6W", /* Better than comp_sugar 20140915 */
    1,
/*
0    0.985600    14    77169.000000    14451.285714    11093.357143    18.726802    14.375406    935.285714    5.428571    15.500000    0.000000    171.785714    941.428571    12.785714    6.714286    0.000000    167.071429    1.122684    0.142160    0.476458    0.000000    1.182352    0.930909    0.311360    0.212733    0.000000    0.831971    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    83.459980    1.133485    0.595238    0.000000    14.811297    1.918693    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    1.880699    4.951874    0.037994    93.129433    0.010000    0.007614
1    0.985600    14    77169.000000    2563.357143    5831.642857    3.321745    7.556976    975.428571    0.000000    0.000000    0.000000    152.571429    959.500000    12.071429    9.642857    0.000000    146.785714    1.082261    0.000000    0.000000    0.000000    1.082261    0.909914    0.298038    0.201581    0.000000    0.781550    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    85.062057    1.070162    0.854863    0.000000    13.012918    2.165653    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.481763    4.305978    0.683891    93.528369    0.187500    0.137056
2    0.985600    14    77169.000000    1454.357143    1516.928571    1.884639    1.965723    971.785714    0.000000    0.000000    0.000000    156.214286    981.000000    5.928571    3.428571    0.000000    137.642857    1.146528    0.000000    0.000000    0.000000    1.146528    0.861625    0.242072    0.096579    0.000000    0.842186    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.968085    0.525583    0.303951    0.000000    12.202381    1.006839    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.721884    4.704914    0.284954    94.288247    0.164835    0.057107
3    0.985600    14    77169.000000    1004.857143    781.857143    1.302151    1.013175    965.928571    0.000000    0.000000    0.000000    162.071429    981.428571    4.000000    2.285714    0.000000    140.285714    1.079324    0.000000    0.000000    0.000000    1.079324    0.778865    0.180683    0.100951    0.000000    0.826112    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    87.006079    0.354610    0.202634    0.000000    12.436677    0.728217    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.563576    4.825228    0.164640    94.446555    0.127451    0.032995
4    0.985600    14    77169.000000    782.785714    548.285714    1.014378    0.710500    962.642857    0.000000    0.000000    0.000000    165.357143    980.071429    2.714286    2.285714    0.000000    142.928571    0.977102    0.000000    0.000000    0.000000    0.977102    0.842801    0.145166    0.117552    0.000000    0.829060    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.885765    0.240628    0.202634    0.000000    12.670973    0.595238    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.493921    4.888551    0.101317    94.516211    0.093023    0.020305
5    0.985600    14    77169.000000    653.928571    400.928571    0.847398    0.519546    960.142857    0.000000    0.000000    0.000000    167.857143    979.500000    2.357143    1.785714    0.000000    144.357143    1.007113    0.000000    0.000000    0.000000    1.007113    0.759241    0.128268    0.079130    0.000000    0.772547    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.835106    0.208967    0.158308    0.000000    12.797619    0.455927    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.341945    4.875887    0.113982    94.668186    0.142857    0.022843
cutoff 0    Affix a -0.751238 b 0.894845
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.168673 b 0.873893
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -2.186093 b 0.839892
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -2.451129 b 0.806491
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -2.687102 b 0.796650
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -2.783724 b 0.782550
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1516.928571 (1454.357143)
rules%         1.965723 (1.884639)
same%stdev     0.861625
ambi1%stdev    0.242072
ambi2%stdev    0.096579
ambi3%stdev    0.000000
diff%stdev     0.842186
same%          86.968085 (86.151216)
ambi1%         0.525583 (0.000000)
ambi2%         0.303951 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.202381 (13.848784)
amb.rules%     1.006839 (0.000000)
false_amb%     0.721884 (0.000000)
false_not_amb% 4.704914 (4.989868)
true_amb%      0.284954 (0.000000)
true_not_amb%  94.288247 (95.010132)
precision       0.164835 (0.000000)
recall          0.057107 (0.000000)

bests[5].suffixonly == [false]
bests[5].langbase == [en6W]
comp = comp_parms0_off
bests[5].rows == [1]
  R->R     W->R     R->W     W->W

0.010390 0.693128 -0.691156 0.204301 0.005382 -0.000404 
*/
//iteration:20.3
/*weight ( used): 1.12271613326426923e+04 suffix only: no */
/* number of nodes: 9840, nodes/line: 1.25675313230391966e-01 weight ( used): 1.12271613326426923e+04 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        1.03898721887023399e-02,        6.93128471853588057e-01,        -6.91156252727169851e-01,       2.04300945832150471e-01,        5.38207657061103278e-03,        -4.03932947156769432e-04
        }}
    };


static bestParms best_enC0 = // English, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*rcount*exp(-2.0*rcount/1.0)
// Aiming at cutoff == 0, because maximum penalty is for 1 word/lemma pairs for a rule.

// Not good, compared with nlC0
    {
    false,
    "enC0", 
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600         14   77169.000000   14451.285714   14440.214286      18.726802      18.712455     935.285714       5.428571      15.500000       0.000000     171.785714     915.857143       0.000000       0.000000       0.000000     212.142857       1.122684       0.142160       0.476458       0.000000       1.182352       1.195915       0.000000       0.000000       0.000000       1.195915      82.915400       0.481256       1.374113       0.000000      15.229230       1.994681      81.193009       0.000000       0.000000       0.000000      18.806991       0.000000       1.956687       4.951874       0.037994      93.053445       0.009615       0.007614       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000
1       0.985600         14   77169.000000    2563.357143    6563.214286       3.321745       8.504988     975.428571       0.000000       0.000000       0.000000     152.571429     930.928571       0.000000       0.000000       0.000000     197.071429       1.082261       0.000000       0.000000       0.000000       1.082261       0.944929       0.000000       0.000000       0.000000       0.944929      86.474164       0.000000       0.000000       0.000000      13.525836       0.000000      82.529129       0.000000       0.000000       0.000000      17.470871       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000
2       0.985600         14   77169.000000    1454.357143    1599.142857       1.884639       2.072261     971.785714       0.000000       0.000000       0.000000     156.214286     931.500000       0.000000       0.000000       0.000000     196.500000       1.146528       0.000000       0.000000       0.000000       1.146528       0.895852       0.000000       0.000000       0.000000       0.895852      86.151216       0.000000       0.000000       0.000000      13.848784       0.000000      82.579787       0.000000       0.000000       0.000000      17.420213       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000
3       0.985600         14   77169.000000    1004.857143     845.428571       1.302151       1.095555     965.928571       0.000000       0.000000       0.000000     162.071429     928.500000       0.000000       0.000000       0.000000     199.500000       1.079324       0.000000       0.000000       0.000000       1.079324       0.951487       0.000000       0.000000       0.000000       0.951487      85.631966       0.000000       0.000000       0.000000      14.368034       0.000000      82.313830       0.000000       0.000000       0.000000      17.686170       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000
4       0.985600         14   77169.000000     782.785714     602.357143       1.014378       0.780569     962.642857       0.000000       0.000000       0.000000     165.357143     918.357143       0.000000       0.000000       0.000000     209.642857       0.977102       0.000000       0.000000       0.000000       0.977102       0.972761       0.000000       0.000000       0.000000       0.972761      85.340679       0.000000       0.000000       0.000000      14.659321       0.000000      81.414640       0.000000       0.000000       0.000000      18.585360       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000
5       0.985600         14   77169.000000     653.928571     465.142857       0.847398       0.602759     960.142857       0.000000       0.000000       0.000000     167.857143     909.357143       0.000000       0.000000       0.000000     218.642857       1.007113       0.000000       0.000000       0.000000       1.007113       0.908489       0.000000       0.000000       0.000000       0.908489      85.119048       0.000000       0.000000       0.000000      14.880952       0.000000      80.616768       0.000000       0.000000       0.000000      19.383232       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000       0.000000       4.989868       0.000000      95.010132       0.000000       0.000000
cutoff 0 Affix  a      -0.366782 b       0.887585: N(rules)=      0.692961*N(trainpairs)^0.887585
         Suffix a      -1.386611 b       0.965747: N(rules)=      0.249921*N(trainpairs)^0.965747
cutoff 1 Affix  a      -0.918458 b       0.862552: N(rules)=      0.399134*N(trainpairs)^0.862552
         Suffix a      -2.548348 b       0.910826: N(rules)=      0.078211*N(trainpairs)^0.910826
cutoff 2 Affix  a      -2.127536 b       0.843785: N(rules)=      0.119131*N(trainpairs)^0.843785
         Suffix a      -2.774658 b       0.877415: N(rules)=      0.062371*N(trainpairs)^0.877415
cutoff 3 Affix  a      -2.457774 b       0.817145: N(rules)=      0.085625*N(trainpairs)^0.817145
         Suffix a      -2.982081 b       0.866045: N(rules)=      0.050687*N(trainpairs)^0.866045
cutoff 4 Affix  a      -2.524994 b       0.789327: N(rules)=      0.080059*N(trainpairs)^0.789327
         Suffix a      -3.164283 b       0.861559: N(rules)=      0.042244*N(trainpairs)^0.861559
cutoff 5 Affix  a      -2.538105 b       0.765681: N(rules)=      0.079016*N(trainpairs)^0.765681
         Suffix a      -3.194672 b       0.847348: N(rules)=      0.040980*N(trainpairs)^0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        yes
cutoff                      2
fraction          9856.000000
iterations                 14
trainlines       77169.000000
rules             1599.142857 (   1454.357143)
rules%               2.072261 (      1.884639)
same%stdev           0.895852
ambi1%stdev          0.000000
ambi2%stdev          0.000000
ambi3%stdev          0.000000
diff%stdev           0.895852
same%               82.579787 (     86.151216)
ambi1%               0.000000 (      0.000000)
ambi2%               0.000000 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               17.420213 (     13.848784)
amb.rules%           0.000000 (      0.000000)
false_amb%           0.000000 (      0.000000)
false_not_amb%       4.989868 (      4.989868)
true_amb%            0.000000 (      0.000000)
true_not_amb%       95.010132 (     95.010132)
precision            0.000000 (      0.000000)
recall               0.000000 (      0.000000)

bests[11].suffixonly == [false]
bests[11].langbase == [enC0]
comp = comp_parms0_off
bests[11].rows == [1]
  R->R     W->R     R->W     W->W

0.107837 0.622363 -0.720085 0.215238 -0.009468 0.189990 
*/
//iteration:20.6
/*weight ( used): 5.65781320119657494e+02 suffix only: no */
/* number of nodes: 14310, nodes/line: inf weight ( used): 5.65781320119657494e+02 blobs 1 lines 0 * fraction 1.00000000000000000e+00 = 0 lines*/
        {{
        1.07837025520926860e-01,    6.22362687879378429e-01,    -7.20085301339376205e-01,    2.15237712457007291e-01,    -9.46806103021151715e-03,    1.89989742917709065e-01
        }}
    };

static bestParms best_en6WS1 = // English, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*exp(-rcount/2.0)
// Aiming at cutoff == 2, because maximum penalty is for 2 word/lemma pairs for a rule.
    {
    false,
    "en6WS1", 
    1,
/*
0    0.985600    14    77169.000000    14451.285714    12738.785714    18.726802    16.507646    935.285714    5.428571    15.500000    0.000000    171.785714    940.571429    10.357143    5.000000    0.000000    172.071429    1.122684    0.142160    0.476458    0.000000    1.182352    0.773413    0.327317    0.278179    0.000000    0.742912    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    83.383992    0.918186    0.443262    0.000000    15.254559    1.519757    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    1.481763    4.951874    0.037994    93.528369    0.012658    0.007614
1    0.985600    14    77169.000000    2563.357143    6752.357143    3.321745    8.750090    975.428571    0.000000    0.000000    0.000000    152.571429    960.857143    10.428571    8.000000    0.000000    148.714286    1.082261    0.000000    0.000000    0.000000    1.082261    0.678602    0.284017    0.235838    0.000000    0.706169    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    85.182371    0.924519    0.709220    0.000000    13.183891    1.855370    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.348784    4.483283    0.506586    93.661348    0.158103    0.101523
2    0.985600    14    77169.000000    1454.357143    1770.214286    1.884639    2.293945    971.785714    0.000000    0.000000    0.000000    156.214286    975.357143    5.285714    3.071429    0.000000    144.285714    1.146528    0.000000    0.000000    0.000000    1.146528    0.881469    0.221095    0.207076    0.000000    0.907704    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.467832    0.468592    0.272290    0.000000    12.791287    0.930851    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.778875    4.837893    0.151976    94.231256    0.088889    0.030457
3    0.985600    14    77169.000000    1004.857143    954.571429    1.302151    1.236988    965.928571    0.000000    0.000000    0.000000    162.071429    975.214286    4.642857    2.142857    0.000000    146.000000    1.079324    0.000000    0.000000    0.000000    1.079324    0.926655    0.157849    0.166244    0.000000    0.867220    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    86.455167    0.411601    0.189970    0.000000    12.943262    0.804205    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.690223    4.875887    0.113982    94.319909    0.076271    0.022843
4    0.985600    14    77169.000000    782.785714    675.214286    1.014378    0.874981    962.642857    0.000000    0.000000    0.000000    165.357143    973.785714    3.428571    1.857143    0.000000    148.928571    0.977102    0.000000    0.000000    0.000000    0.977102    0.993406    0.154392    0.119736    0.000000    0.988089    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.328521    0.303951    0.164640    0.000000    13.202888    0.658561    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.557244    4.888551    0.101317    94.452888    0.083333    0.020305
5    0.985600    14    77169.000000    653.928571    510.214286    0.847398    0.661165    960.142857    0.000000    0.000000    0.000000    167.857143    971.714286    3.357143    2.000000    0.000000    150.928571    1.007113    0.000000    0.000000    0.000000    1.007113    0.941701    0.192374    0.115327    0.000000    0.916345    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.144883    0.297619    0.177305    0.000000    13.380193    0.633232    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.550912    4.907548    0.082320    94.459220    0.069519    0.016497
cutoff 0    Affix a -0.611390 b 0.895684
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.154468 b 0.886261
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -2.266144 b 0.864527
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -2.617687 b 0.842355
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -2.846212 b 0.831893
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -2.862054 b 0.808661
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1770.214286 (1454.357143)
rules%         2.293945 (1.884639)
same%stdev     0.881469
ambi1%stdev    0.221095
ambi2%stdev    0.207076
ambi3%stdev    0.000000
diff%stdev     0.907704
same%          86.467832 (86.151216)
ambi1%         0.468592 (0.000000)
ambi2%         0.272290 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.791287 (13.848784)
amb.rules%     0.930851 (0.000000)
false_amb%     0.778875 (0.000000)
false_not_amb% 4.837893 (4.989868)
true_amb%      0.151976 (0.000000)
true_not_amb%  94.231256 (95.010132)
precision       0.088889 (0.000000)
recall          0.030457 (0.000000)

bests[8].suffixonly == [false]
bests[8].langbase == [en6WS1]
comp = comp_parms0_off
bests[8].rows == [1]
  R->R     W->R     R->W     W->W

0.099984 0.698758 -0.658037 0.246512 0.028495 0.084549 
*/
//iteration:20.7
/*weight ( used): 4.15520713829023862e+03 suffix only: no */
/* number of nodes: 11740, nodes/line: 1.49941887939512380e-01 weight ( used): 4.15520713829023862e+03 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        9.99843486891121014e-02,        6.98757557409581676e-01,        -6.58036595438773131e-01,       2.46512251235911778e-01,        2.84952417970883061e-02,        8.45492532733789354e-02
        }}
    };

static bestParms best_en6WS2 = // English, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*exp(-rcount/3.0)
// Aiming at cutoff == 2, because maximum penalty is for 3 word/lemma pairs for a rule.
    {
    false,
    "en6WS2", 
    1,
/*
0    0.985600    14    77169.000000    14451.285714    11399.714286    18.726802    14.772401    935.285714    5.428571    15.500000    0.000000    171.785714    944.285714    11.000000    5.214286    0.000000    167.500000    1.122684    0.142160    0.476458    0.000000    1.182352    0.773971    0.243407    0.200292    0.000000    0.705588    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    83.713273    0.975177    0.462259    0.000000    14.849291    1.621074    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    1.583080    4.951874    0.037994    93.427052    0.011858    0.007614
1    0.985600    14    77169.000000    2563.357143    5901.071429    3.321745    7.646946    975.428571    0.000000    0.000000    0.000000    152.571429    963.714286    10.000000    7.500000    0.000000    146.785714    1.082261    0.000000    0.000000    0.000000    1.082261    0.701875    0.307101    0.154532    0.000000    0.566244    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    85.435664    0.886525    0.664894    0.000000    13.012918    1.804711    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.329787    4.514944    0.474924    93.680344    0.151515    0.095178
2    0.985600    14    77169.000000    1454.357143    1501.357143    1.884639    1.945544    971.785714    0.000000    0.000000    0.000000    156.214286    980.785714    5.071429    2.571429    0.000000    139.571429    1.146528    0.000000    0.000000    0.000000    1.146528    0.858236    0.207076    0.137842    0.000000    0.818972    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.949088    0.449595    0.227964    0.000000    12.373354    0.848531    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.645897    4.787234    0.202634    94.364235    0.135593    0.040609
3    0.985600    14    77169.000000    1004.857143    741.285714    1.302151    0.960600    965.928571    0.000000    0.000000    0.000000    162.071429    981.428571    2.714286    2.000000    0.000000    141.857143    1.079324    0.000000    0.000000    0.000000    1.079324    0.794999    0.112291    0.098351    0.000000    0.790532    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    87.006079    0.240628    0.177305    0.000000    12.575988    0.563576    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.436930    4.863222    0.126646    94.573202    0.126582    0.025381
4    0.985600    14    77169.000000    782.785714    524.000000    1.014378    0.679029    962.642857    0.000000    0.000000    0.000000    165.357143    980.928571    2.142857    1.785714    0.000000    143.142857    0.977102    0.000000    0.000000    0.000000    0.977102    0.852784    0.114576    0.110937    0.000000    0.842365    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.961753    0.189970    0.158308    0.000000    12.689970    0.481256    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.367275    4.875887    0.113982    94.642857    0.134328    0.022843
5    0.985600    14    77169.000000    653.928571    401.428571    0.847398    0.520194    960.142857    0.000000    0.000000    0.000000    167.857143    979.857143    2.071429    1.214286    0.000000    144.857143    1.007113    0.000000    0.000000    0.000000    1.007113    0.740788    0.122763    0.079130    0.000000    0.725116    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.866768    0.183637    0.107649    0.000000    12.841945    0.360942    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.246960    4.875887    0.113982    94.763171    0.187500    0.022843
cutoff 0    Affix a -0.743918 b 0.896681
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.197041 b 0.877837
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -2.161264 b 0.837071
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -2.305112 b 0.789647
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -2.562766 b 0.782704
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -2.797124 b 0.784175
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1501.357143 (1454.357143)
rules%         1.945544 (1.884639)
same%stdev     0.858236
ambi1%stdev    0.207076
ambi2%stdev    0.137842
ambi3%stdev    0.000000
diff%stdev     0.818972
same%          86.949088 (86.151216)
ambi1%         0.449595 (0.000000)
ambi2%         0.227964 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.373354 (13.848784)
amb.rules%     0.848531 (0.000000)
false_amb%     0.645897 (0.000000)
false_not_amb% 4.787234 (4.989868)
true_amb%      0.202634 (0.000000)
true_not_amb%  94.364235 (95.010132)
precision       0.135593 (0.000000)
recall          0.040609 (0.000000)

bests[7].suffixonly == [false]
bests[7].langbase == [en6WS2]
comp = comp_parms0_off
bests[7].rows == [1]
  R->R     W->R     R->W     W->W

0.040853 0.673407 -0.708980 0.204693 0.017152 -0.002812 
*/
//iteration:20.3
/*weight ( used): 5.80760738773845060e+03 suffix only: no */
/* number of nodes: 10161, nodes/line: 1.29775087168090719e-01 weight ( used): 5.80760738773845060e+03 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        4.08527421595744900e-02,        6.73406780072449362e-01,        -7.08980349688113631e-01,       2.04692793670780471e-01,        1.71516852929795953e-02,        -2.81170366166549109e-03
        }}
    };

static bestParms best_en6WS3 = // English, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*exp(-rcount/4.0)
// Aiming at cutoff == 3, because maximum penalty is for 4 word/lemma pairs for a rule.
    {
    false,
    "en6WS3", 
    1,
/*
0    0.985600    14    77169.000000    14451.285714    11318.285714    18.726802    14.666881    935.285714    5.428571    15.500000    0.000000    171.785714    942.142857    11.285714    6.071429    0.000000    168.500000    1.122684    0.142160    0.476458    0.000000    1.182352    1.015482    0.279109    0.175469    0.000000    0.819736    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    83.523303    1.000507    0.538247    0.000000    14.937943    1.728723    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    1.690729    4.951874    0.037994    93.319402    0.011111    0.007614
1    0.985600    14    77169.000000    2563.357143    5873.857143    3.321745    7.611680    975.428571    0.000000    0.000000    0.000000    152.571429    962.071429    11.571429    8.642857    0.000000    145.714286    1.082261    0.000000    0.000000    0.000000    1.082261    0.862652    0.279727    0.216057    0.000000    0.687956    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    85.290020    1.025836    0.766211    0.000000    12.917933    2.045339    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.481763    4.426292    0.563576    93.528369    0.159785    0.112944
2    0.985600    14    77169.000000    1454.357143    1506.500000    1.884639    1.952209    971.785714    0.000000    0.000000    0.000000    156.214286    980.714286    5.285714    3.357143    0.000000    138.642857    1.146528    0.000000    0.000000    0.000000    1.146528    0.840620    0.229151    0.132898    0.000000    0.881469    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.942756    0.468592    0.297619    0.000000    12.291033    0.930851    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.728217    4.787234    0.202634    94.281915    0.122137    0.040609
3    0.985600    14    77169.000000    1004.857143    772.357143    1.302151    1.000865    965.928571    0.000000    0.000000    0.000000    162.071429    980.642857    3.428571    2.571429    0.000000    141.357143    1.079324    0.000000    0.000000    0.000000    1.079324    0.806994    0.165723    0.108378    0.000000    0.874583    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    86.936424    0.303951    0.227964    0.000000    12.531662    0.677558    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.550912    4.863222    0.126646    94.459220    0.103093    0.025381
4    0.985600    14    77169.000000    782.785714    543.642857    1.014378    0.704483    962.642857    0.000000    0.000000    0.000000    165.357143    979.571429    2.642857    2.214286    0.000000    143.571429    0.977102    0.000000    0.000000    0.000000    0.977102    0.841545    0.149994    0.110937    0.000000    0.876729    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.841439    0.234296    0.196302    0.000000    12.727964    0.563576    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.462259    4.888551    0.101317    94.547872    0.098765    0.020305
5    0.985600    14    77169.000000    653.928571    400.714286    0.847398    0.519268    960.142857    0.000000    0.000000    0.000000    167.857143    978.928571    2.071429    1.857143    0.000000    145.142857    1.007113    0.000000    0.000000    0.000000    1.007113    0.723954    0.122763    0.109172    0.000000    0.743232    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.784448    0.183637    0.164640    0.000000    12.867275    0.417933    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.316616    4.888551    0.101317    94.693516    0.137931    0.020305
cutoff 0    Affix a -0.755377 b 0.897222
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.217108 b 0.879540
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -2.116345 b 0.832979
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -2.356162 b 0.796254
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -2.473309 b 0.774972
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -2.713573 b 0.775252
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1506.500000 (1454.357143)
rules%         1.952209 (1.884639)
same%stdev     0.840620
ambi1%stdev    0.229151
ambi2%stdev    0.132898
ambi3%stdev    0.000000
diff%stdev     0.881469
same%          86.942756 (86.151216)
ambi1%         0.468592 (0.000000)
ambi2%         0.297619 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.291033 (13.848784)
amb.rules%     0.930851 (0.000000)
false_amb%     0.728217 (0.000000)
false_not_amb% 4.787234 (4.989868)
true_amb%      0.202634 (0.000000)
true_not_amb%  94.281915 (95.010132)
precision       0.122137 (0.000000)
recall          0.040609 (0.000000)

bests[9].suffixonly == [false]
bests[9].langbase == [en6WS3]
comp = comp_parms0_off
bests[9].rows == [1]
  R->R     W->R     R->W     W->W

0.020806 0.683019 -0.702896 0.197456 -0.000308 0.001105 
*/
//iteration:20.9
/*weight ( used): 7.45427103780143261e+03 suffix only: no */
/* number of nodes: 10102, nodes/line: 1.29021546163965412e-01 weight ( used): 7.45427103780143261e+03 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        2.08057417098401302e-02,    6.83018660111897935e-01,    -7.02895839518911436e-01,    1.97455702463298050e-01,    -3.08022357388210524e-04,    1.10477364594036725e-03
        }}

    };

static bestParms best_enS2 = // English, ambiguous training pairs in training set derived from CELEX
// weight function (graph.h) rcount*rcount*exp(-2.0*rcount/3.0)
// Aiming at cutoff == 2, because maximum penalty is for 3 word/lemma pairs for a rule.
    {
    false,
    "enS2", 
    1,
/*
0    0.985600    14    77169.000000    14451.285714    11784.142857    18.726802    15.270566    935.285714    5.428571    15.500000    0.000000    171.785714    945.071429    9.500000    4.928571    0.000000    168.500000    1.122684    0.142160    0.476458    0.000000    1.182352    1.057826    0.308574    0.188748    0.000000    0.849436    82.915400    0.481256    1.374113    0.000000    15.229230    1.994681    83.782928    0.842199    0.436930    0.000000    14.937943    1.437437    1.956687    4.951874    0.037994    93.053445    0.009615    0.007614    1.399443    4.951874    0.037994    93.610689    0.013393    0.007614
1    0.985600    14    77169.000000    2563.357143    6187.785714    3.321745    8.018486    975.428571    0.000000    0.000000    0.000000    152.571429    965.642857    9.428571    6.642857    0.000000    146.285714    1.082261    0.000000    0.000000    0.000000    1.082261    1.018306    0.294467    0.198559    0.000000    0.846354    86.474164    0.000000    0.000000    0.000000    13.525836    0.000000    85.606636    0.835866    0.588906    0.000000    12.968592    1.621074    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    1.127153    4.495947    0.493921    93.882979    0.179724    0.098985
2    0.985600    14    77169.000000    1454.357143    1461.000000    1.884639    1.893247    971.785714    0.000000    0.000000    0.000000    156.214286    980.928571    4.500000    2.928571    0.000000    139.642857    1.146528    0.000000    0.000000    0.000000    1.146528    0.981951    0.179844    0.149417    0.000000    1.001545    86.151216    0.000000    0.000000    0.000000    13.848784    0.000000    86.961753    0.398936    0.259625    0.000000    12.379686    0.816869    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.595238    4.768237    0.221631    94.414894    0.156951    0.044416
3    0.985600    14    77169.000000    1004.857143    740.500000    1.302151    0.959582    965.928571    0.000000    0.000000    0.000000    162.071429    981.214286    3.000000    1.785714    0.000000    142.000000    1.079324    0.000000    0.000000    0.000000    1.079324    0.956827    0.155507    0.061994    0.000000    0.977345    85.631966    0.000000    0.000000    0.000000    14.368034    0.000000    86.987082    0.265957    0.158308    0.000000    12.588652    0.607903    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.455927    4.837893    0.151976    94.554205    0.142857    0.030457
4    0.985600    14    77169.000000    782.785714    542.357143    1.014378    0.702817    962.642857    0.000000    0.000000    0.000000    165.357143    980.142857    2.357143    1.357143    0.000000    144.142857    0.977102    0.000000    0.000000    0.000000    0.977102    0.955992    0.145908    0.102015    0.000000    0.992602    85.340679    0.000000    0.000000    0.000000    14.659321    0.000000    86.892097    0.208967    0.120314    0.000000    12.778622    0.487589    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.360942    4.863222    0.126646    94.649189    0.149254    0.025381
5    0.985600    14    77169.000000    653.928571    418.785714    0.847398    0.542686    960.142857    0.000000    0.000000    0.000000    167.857143    978.928571    2.642857    1.285714    0.000000    145.142857    1.007113    0.000000    0.000000    0.000000    1.007113    0.893632    0.149994    0.117552    0.000000    0.880954    85.119048    0.000000    0.000000    0.000000    14.880952    0.000000    86.784448    0.234296    0.113982    0.000000    12.867275    0.468592    0.000000    4.989868    0.000000    95.010132    0.000000    0.000000    0.354610    4.875887    0.113982    94.655522    0.138462    0.022843
cutoff 0    Affix a -0.666235 b 0.892776
            Suffix a -1.386611 b 0.965747
cutoff 1    Affix a -1.147043 b 0.876909
            Suffix a -2.548348 b 0.910826
cutoff 2    Affix a -1.806319 b 0.801931
            Suffix a -2.774658 b 0.877415
cutoff 3    Affix a -1.747592 b 0.737090
            Suffix a -2.982081 b 0.866045
cutoff 4    Affix a -1.922767 b 0.724786
            Suffix a -3.164283 b 0.861559
cutoff 5    Affix a -1.957773 b 0.707745
            Suffix a -3.194672 b 0.847348

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      14
trainlines      77169.000000
rules           1461.000000 (1454.357143)
rules%         1.893247 (1.884639)
same%stdev     0.981951
ambi1%stdev    0.179844
ambi2%stdev    0.149417
ambi3%stdev    0.000000
diff%stdev     1.001545
same%          86.961753 (86.151216)
ambi1%         0.398936 (0.000000)
ambi2%         0.259625 (0.000000)
ambi3%         0.000000 (0.000000)
diff%          12.379686 (13.848784)
amb.rules%     0.816869 (0.000000)
false_amb%     0.595238 (0.000000)
false_not_amb% 4.768237 (4.989868)
true_amb%      0.221631 (0.000000)
true_not_amb%  94.414894 (95.010132)
precision       0.156951 (0.000000)
recall          0.044416 (0.000000)

bests[10].suffixonly == [false]
bests[10].langbase == [enS2]
comp = comp_parms0_off
bests[10].rows == [1]
  R->R     W->R     R->W     W->W

0.179015 0.692631 -0.665544 0.161160 0.137396 -0.020423 
*/
//iteration:20.4
/*weight ( used): 4.66981464254577713e+03 suffix only: no */
/* number of nodes: 10378, nodes/line: 1.32546585437500808e-01 weight ( used): 4.66981464254577713e+03 blobs 1 lines 78297 * fraction 1.00000000000000000e+00 = 78297 lines*/
        {{
        1.79014987751130839e-01,    6.92630686252189487e-01,    -6.65544221444446804e-01,    1.61160030170824808e-01,    1.37395843876723128e-01,    -2.04226023055548052e-02
        }}
    };



static bestParms best_deC0 =
    {
    false,
    "deC0",
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          7  313549.000000   42991.285714   34219.714286      13.711186      10.913674    4070.857143      15.000000      17.142857       0.428571     478.571429    4116.142857      19.142857      16.714286       0.000000     430.000000       0.487918       0.087298       0.075000       0.017171       0.477629       0.326014       0.159099       0.106706       0.000000       0.267592      88.844547       0.327368       0.374135       0.009353      10.444597       0.841803      89.832886       0.417784       0.364781       0.000000       9.384548       0.957162       0.779448       1.686724       0.062356      97.471472       0.038462       0.035651       0.894806       1.686724       0.062356      97.356114       0.033670       0.035651
1       0.985600          7  313549.000000   12232.142857   18745.285714       3.901190       5.978423    4059.571429       0.000000       0.000000       0.000000     522.428571    4123.428571      15.285714      19.000000       0.000000     424.285714       0.422039       0.000000       0.000000       0.000000       0.422039       0.387932       0.083310       0.110568       0.000000       0.367732      88.598242       0.000000       0.000000       0.000000      11.401758       0.000000      89.991894       0.333604       0.414666       0.000000       9.259837       0.882335       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.760741       1.627486       0.121594      97.490179       0.074004       0.069519
2       0.985600          7  313549.000000    7700.142857    8184.428571       2.455802       2.610255    4015.285714       0.000000       0.000000       0.000000     566.714286    4115.142857      11.714286      12.857143       0.000000     442.285714       0.433335       0.000000       0.000000       0.000000       0.433335       0.510655       0.053033       0.079121       0.000000       0.478270      87.631727       0.000000       0.000000       0.000000      12.368273       0.000000      89.811062       0.255659       0.280601       0.000000       9.652678       0.654736       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.579909       1.674253       0.074827      97.671011       0.060606       0.042781
3       0.985600          7  313549.000000    5708.857143    5542.571429       1.820722       1.767689    3974.571429       0.000000       0.000000       0.000000     607.428571    4092.428571       8.857143       9.857143       0.000000     470.857143       0.478957       0.000000       0.000000       0.000000       0.478957       0.484232       0.042597       0.107236       0.000000       0.460635      86.743156       0.000000       0.000000       0.000000      13.256844       0.000000      89.315333       0.193303       0.215128       0.000000      10.276236       0.514435       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.458315       1.692960       0.056120      97.792605       0.057692       0.032086
4       0.985600          7  313549.000000    4622.571429    4432.714286       1.474274       1.413723    3931.571429       0.000000       0.000000       0.000000     650.428571    4069.000000       9.571429       7.857143       0.000000     495.571429       0.479454       0.000000       0.000000       0.000000       0.479454       0.539468       0.050176       0.092102       0.000000       0.502100      85.804702       0.000000       0.000000       0.000000      14.195298       0.000000      88.804016       0.208892       0.171478       0.000000      10.815614       0.473904       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.420902       1.696078       0.053002      97.830018       0.059233       0.030303
5       0.985600          7  313549.000000    3814.857143    3738.571429       1.216670       1.192340    3900.714286       0.000000       0.000000       0.000000     681.285714    4050.428571      10.000000       8.142857       0.000000     513.428571       0.522987       0.000000       0.000000       0.000000       0.522987       0.520727       0.073472       0.105745       0.000000       0.453064      85.131259       0.000000       0.000000       0.000000      14.868741       0.000000      88.398703       0.218245       0.177714       0.000000      11.205338       0.498846       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.445844       1.696078       0.053002      97.805076       0.056106       0.030303
cutoff 0 Affix  a       0.696894 b       0.771507: N(rules)=      2.007508*N(trainpairs)^0.771507
         Suffix a       0.255148 b       0.819252: N(rules)=      1.290653*N(trainpairs)^0.819252
cutoff 1 Affix  a       0.424828 b       0.745235: N(rules)=      1.529327*N(trainpairs)^0.745235
         Suffix a      -1.476990 b       0.863095: N(rules)=      0.228324*N(trainpairs)^0.863095
cutoff 2 Affix  a      -0.660663 b       0.770330: N(rules)=      0.516509*N(trainpairs)^0.770330
         Suffix a      -2.565046 b       0.914071: N(rules)=      0.076916*N(trainpairs)^0.914071
cutoff 3 Affix  a      -1.199171 b       0.781329: N(rules)=      0.301444*N(trainpairs)^0.781329
         Suffix a      -3.132748 b       0.934644: N(rules)=      0.043598*N(trainpairs)^0.934644
cutoff 4 Affix  a      -1.435560 b       0.778725: N(rules)=      0.237982*N(trainpairs)^0.778725
         Suffix a      -3.498437 b       0.945496: N(rules)=      0.030245*N(trainpairs)^0.945496
cutoff 5 Affix  a      -1.957380 b       0.805932: N(rules)=      0.141228*N(trainpairs)^0.805932
         Suffix a      -3.614366 b       0.938124: N(rules)=      0.026934*N(trainpairs)^0.938124

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        no
cutoff                      1
fraction          9856.000000
iterations                  7
trainlines      313549.000000
rules            18745.285714 (  12232.142857)
rules%               5.978423 (      3.901190)
same%stdev           0.387932
ambi1%stdev          0.083310
ambi2%stdev          0.110568
ambi3%stdev          0.000000
diff%stdev           0.367732
same%               89.991894 (     88.598242)
ambi1%               0.333604 (      0.000000)
ambi2%               0.414666 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%                9.259837 (     11.401758)
amb.rules%           0.882335 (      0.000000)
false_amb%           0.760741 (      0.000000)
false_not_amb%       1.627486 (      1.749080)
true_amb%            0.121594 (      0.000000)
true_not_amb%       97.490179 (     98.250920)
precision            0.074004 (      0.000000)
recall               0.069519 (      0.000000)

bests[2].suffixonly == [false]
bests[2].langbase == [deC0]
comp = comp_parms0_off
bests[2].rows == [1]
  R->R     W->R     R->W     W->W

0.070270 0.738908 -0.655311 0.134438 0.039578 0.002103 
*/
//iteration:20.9
/*weight ( used): 1.57191628625387852e+03 suffix only: no */
/* number of nodes: 33777, nodes/line: inf weight ( used): 1.57191628625387852e+03 blobs 1 lines 0 * fraction 1.00000000000000000e+00 = 0 lines*/
        {{
            7.02703009236997217e-02,    7.38908457216426395e-01,    -6.55310616504478860e-01,    1.34437786206644400e-01,    3.95781745400746385e-02,    2.10292619579706538e-03
        }}
    };

static bestParms best_deC1 =
    {
    false,
    "deC1",
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          7  313549.000000   42991.285714   41378.428571      13.711186      13.196798    4070.857143      15.000000      17.142857       0.428571     478.571429    4113.285714      24.285714      18.857143       0.000000     425.571429       0.487918       0.087298       0.075000       0.017171       0.477629       0.351168       0.068686       0.095488       0.000000       0.394426      88.844547       0.327368       0.374135       0.009353      10.444597       0.841803      89.770531       0.530024       0.411548       0.000000       9.287897       1.147347       0.779448       1.686724       0.062356      97.471472       0.038462       0.035651       1.078755       1.680489       0.068591      97.172164       0.030812       0.039216
1       0.985600          7  313549.000000   12232.142857   21711.000000       3.901190       6.924277    4059.571429       0.000000       0.000000       0.000000     522.428571    4113.428571      22.857143      21.857143       0.000000     423.857143       0.422039       0.000000       0.000000       0.000000       0.422039       0.329750       0.079121       0.089479       0.000000       0.321107      88.598242       0.000000       0.000000       0.000000      11.401758       0.000000      89.773648       0.498846       0.477022       0.000000       9.250483       1.150465       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.991457       1.590073       0.159007      97.259462       0.074236       0.090909
2       0.985600          7  313549.000000    7700.142857    9295.285714       2.455802       2.964540    4015.285714       0.000000       0.000000       0.000000     566.714286    4109.714286      16.571429      17.285714       0.000000     438.428571       0.433335       0.000000       0.000000       0.000000       0.433335       0.426131       0.117432       0.060053       0.000000       0.397034      87.631727       0.000000       0.000000       0.000000      12.368273       0.000000      89.692586       0.361664       0.377253       0.000000       9.568498       0.876099       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.773212       1.646193       0.102887      97.477708       0.062382       0.058824
3       0.985600          7  313549.000000    5708.857143    6627.142857       1.820722       2.113591    3974.571429       0.000000       0.000000       0.000000     607.428571    4081.857143      16.714286      14.571429       0.000000     468.857143       0.478957       0.000000       0.000000       0.000000       0.478957       0.367454       0.099786       0.041519       0.000000       0.406434      86.743156       0.000000       0.000000       0.000000      13.256844       0.000000      89.084617       0.364781       0.318015       0.000000      10.232587       0.826214       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.723327       1.646193       0.102887      97.527592       0.066398       0.058824
4       0.985600          7  313549.000000    4622.571429    5333.857143       1.474274       1.701124    3931.571429       0.000000       0.000000       0.000000     650.428571    4052.000000      16.857143      14.428571       0.000000     498.714286       0.479454       0.000000       0.000000       0.000000       0.479454       0.434120       0.084928       0.046905       0.000000       0.429101      85.804702       0.000000       0.000000       0.000000      14.195298       0.000000      88.432999       0.367899       0.314897       0.000000      10.884205       0.835568       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.713974       1.627486       0.121594      97.536946       0.078471       0.069519
5       0.985600          7  313549.000000    3814.857143    4496.857143       1.216670       1.434180    3900.714286       0.000000       0.000000       0.000000     681.285714    4033.714286      15.857143      14.000000       0.000000     518.428571       0.522987       0.000000       0.000000       0.000000       0.522987       0.378162       0.070639       0.035639       0.000000       0.350295      85.131259       0.000000       0.000000       0.000000      14.868741       0.000000      88.033922       0.346075       0.305543       0.000000      11.314460       0.810625       0.000000       1.749080       0.000000      98.250920       0.000000       0.000000       0.682796       1.621251       0.127829      97.568124       0.085595       0.073084
cutoff 0 Affix  a       0.919873 b       0.765671: N(rules)=      2.508971*N(trainpairs)^0.765671
         Suffix a       0.255148 b       0.819252: N(rules)=      1.290653*N(trainpairs)^0.819252
cutoff 1 Affix  a       0.525258 b       0.746897: N(rules)=      1.690895*N(trainpairs)^0.746897
         Suffix a      -1.476990 b       0.863095: N(rules)=      0.228324*N(trainpairs)^0.863095
cutoff 2 Affix  a      -1.113654 b       0.815502: N(rules)=      0.328357*N(trainpairs)^0.815502
         Suffix a      -2.565046 b       0.914071: N(rules)=      0.076916*N(trainpairs)^0.914071
cutoff 3 Affix  a      -1.594497 b       0.824227: N(rules)=      0.203011*N(trainpairs)^0.824227
         Suffix a      -3.132748 b       0.934644: N(rules)=      0.043598*N(trainpairs)^0.934644
cutoff 4 Affix  a      -2.182396 b       0.851581: N(rules)=      0.112771*N(trainpairs)^0.851581
         Suffix a      -3.498437 b       0.945496: N(rules)=      0.030245*N(trainpairs)^0.945496
cutoff 5 Affix  a      -2.728401 b       0.880592: N(rules)=      0.065324*N(trainpairs)^0.880592
         Suffix a      -3.614366 b       0.938124: N(rules)=      0.026934*N(trainpairs)^0.938124

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        no
cutoff                      1
fraction          9856.000000
iterations                  7
trainlines      313549.000000
rules            21711.000000 (  12232.142857)
rules%               6.924277 (      3.901190)
same%stdev           0.329750
ambi1%stdev          0.079121
ambi2%stdev          0.089479
ambi3%stdev          0.000000
diff%stdev           0.321107
same%               89.773648 (     88.598242)
ambi1%               0.498846 (      0.000000)
ambi2%               0.477022 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%                9.250483 (     11.401758)
amb.rules%           1.150465 (      0.000000)
false_amb%           0.991457 (      0.000000)
false_not_amb%       1.590073 (      1.749080)
true_amb%            0.159007 (      0.000000)
true_not_amb%       97.259462 (     98.250920)
precision            0.074236 (      0.000000)
recall               0.090909 (      0.000000)

bests[4].suffixonly == [false]
bests[4].langbase == [deC1]
comp = comp_parms0_off
bests[4].rows == [1]
  R->R     W->R     R->W     W->W

0.155770 0.629146 -0.737028 0.165459 0.079044 0.055455 
*/
//iteration:20.8
/*weight ( used): 6.78787945680193934e+03 suffix only: no */
/* number of nodes: 36815, nodes/line: 1.15722768293564607e-01 weight ( used): 6.78787945680193934e+03 blobs 1 lines 318131 * fraction 1.00000000000000000e+00 = 318131 lines*/
        {{
        1.55770018417734107e-01,    6.29146340759315748e-01,    -7.37028188136959139e-01,    1.65459089224974654e-01,    7.90444951020151526e-02,    5.54553040724012833e-02
        }}
    };





static bestParms best_elC0 =
    {
    false,
    "elC0",
    1,
/* First time with ambiguity > 2 (2015.02.16)
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          5  556568.000000   91660.600000   77694.400000      16.468895      13.959552    6965.200000      17.600000      11.000000       0.200000    1138.000000    6947.200000      51.200000      25.400000       0.000000    1108.200000       0.210385       0.033225       0.035852       0.005499       0.228739       0.347749       0.047148       0.053319       0.000000       0.309804      85.651746       0.216429       0.135268       0.002459      13.994097       0.450074      85.430398       0.629611       0.312346       0.000000      13.627644       1.160846       0.356616       5.627152       0.093458      93.922774       0.115854       0.016337       1.032956       5.592720       0.127890      93.246434       0.058296       0.022356
1       0.985600          5  556568.000000   20985.400000   40104.800000       3.770501       7.205732    7066.800000       0.000000       0.000000       0.000000    1065.200000    6970.600000      96.600000     111.000000       0.000000     953.800000       0.276261       0.000000       0.000000       0.000000       0.276261       0.407162       0.078548       0.066222       0.000000       0.431136      86.901131       0.000000       0.000000       0.000000      13.098869       0.000000      85.718151       1.187900       1.364978       0.000000      11.728972       2.843089       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.976390       3.853910       1.866699      93.303000       0.488731       0.326311
2       0.985600          5  556568.000000   12589.600000   13686.600000       2.262006       2.459107    6981.000000       0.000000       0.000000       0.000000    1151.000000    6986.200000     134.000000     193.800000       0.200000     817.800000       0.350629       0.000000       0.000000       0.000000       0.350629       0.325745       0.188912       0.167575       0.005499       0.256093      85.846040       0.000000       0.000000       0.000000      14.153960       0.000000      85.909985       1.647811       2.383178       0.002459      10.056567       4.323660       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.782095       2.179046       3.541564      93.497295       0.693642       0.619089
3       0.985600          5  556568.000000    8924.600000    8586.000000       1.603506       1.542669    6909.800000       0.000000       0.000000       0.000000    1222.200000    6959.800000     148.200000     187.200000       0.000000     836.800000       0.399930       0.000000       0.000000       0.000000       0.399930       0.325397       0.194123       0.148942       0.000000       0.269472      84.970487       0.000000       0.000000       0.000000      15.029513       0.000000      85.585342       1.822430       2.302017       0.000000      10.290212       4.402361       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.686178       2.004427       3.716183      93.593212       0.730304       0.649613
4       0.985600          5  556568.000000    6978.800000    6381.800000       1.253899       1.146634    6859.600000       0.000000       0.000000       0.000000    1272.400000    6931.800000     155.600000     181.000000       0.200000     863.400000       0.294951       0.000000       0.000000       0.000000       0.294951       0.305255       0.158818       0.170172       0.005499       0.233580      84.353173       0.000000       0.000000       0.000000      15.646827       0.000000      85.241023       1.913428       2.225775       0.002459      10.617314       4.404820       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.673881       1.989670       3.730939      93.605509       0.734625       0.652193
5       0.985600          5  556568.000000    5790.000000    5232.800000       1.040304       0.940191    6802.000000       0.000000       0.000000       0.000000    1330.000000    6908.400000     161.800000     172.000000       0.000000     889.800000       0.253063       0.000000       0.000000       0.000000       0.253063       0.295079       0.147923       0.123584       0.000000       0.229300      83.644860       0.000000       0.000000       0.000000      16.355140       0.000000      84.953271       1.989670       2.115101       0.000000      10.941958       4.387605       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.659124       1.992130       3.728480      93.620266       0.738791       0.651763
cutoff 0 Affix  a       0.596157 b       0.808684: N(rules)=      1.815130*N(trainpairs)^0.808684
         Suffix a      -0.054097 b       0.862983: N(rules)=      0.947341*N(trainpairs)^0.862983
cutoff 1 Affix  a       0.030996 b       0.801290: N(rules)=      1.031482*N(trainpairs)^0.801290
         Suffix a      -1.391083 b       0.856475: N(rules)=      0.248806*N(trainpairs)^0.856475
cutoff 2 Affix  a      -0.559577 b       0.760960: N(rules)=      0.571451*N(trainpairs)^0.760960
         Suffix a      -1.801956 b       0.843674: N(rules)=      0.164976*N(trainpairs)^0.843674
cutoff 3 Affix  a      -0.787383 b       0.742457: N(rules)=      0.455034*N(trainpairs)^0.742457
         Suffix a      -1.870778 b       0.819764: N(rules)=      0.154004*N(trainpairs)^0.819764
cutoff 4 Affix  a      -0.907361 b       0.730001: N(rules)=      0.403588*N(trainpairs)^0.730001
         Suffix a      -1.821582 b       0.794623: N(rules)=      0.161770*N(trainpairs)^0.794623
cutoff 5 Affix  a      -0.992606 b       0.721078: N(rules)=      0.370610*N(trainpairs)^0.721078
         Suffix a      -1.762567 b       0.774137: N(rules)=      0.171604*N(trainpairs)^0.774137

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        yes
cutoff                      2
fraction          9856.000000
iterations                  5
trainlines      556568.000000
rules            13686.600000 (  12589.600000)
rules%               2.459107 (      2.262006)
same%stdev           0.325745
ambi1%stdev          0.188912
ambi2%stdev          0.167575
ambi3%stdev          0.005499
diff%stdev           0.256093
same%               85.909985 (     85.846040)
ambi1%               1.647811 (      0.000000)
ambi2%               2.383178 (      0.000000)
ambi3%               0.002459 (      0.000000)
diff%               10.056567 (     14.153960)
amb.rules%           4.323660 (      0.000000)
false_amb%           0.782095 (      0.000000)
false_not_amb%       2.179046 (      5.720610)
true_amb%            3.541564 (      0.000000)
true_not_amb%       93.497295 (     94.279390)
precision            0.693642 (      0.000000)
recall               0.619089 (      0.000000)

bests[5].suffixonly == [false]
bests[5].langbase == [elC0]
comp = comp_parms0_off
bests[5].rows == [1]
  R->R     W->R     R->W     W->W

0.162988 0.723313 -0.629186 0.197114 0.003088 0.124558
*/
//iteration:20.5
/*weight ( used): 2.88408936339958746e+03 suffix only: no */
/* number of nodes: 71647, nodes/line: 1.26876217460598545e-01 weight ( used): 2.88408936339958746e+03 blobs 1 lines 564700 * fraction 1.00000000000000000e+00 = 564700 lines*/
        {{
        1.62988138483171396e-01,    7.23313127184479043e-01,    -6.29185831644276772e-01,    1.97113838764327726e-01,    3.08765150599594205e-03,    1.24558327663259708e-01
        }}
    };


static bestParms best_elC1 =
    {
    false,
    "elC1",
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%
s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_tr
ue_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          5  556568.000000   91660.600000   78373.000000      16.468895      14.081478    6965.200000      17.600000      11.000000       0.200000    1138.000000    6959.600000      46.400000      24.400000       0.000000    1101.600000       0.210385       0.033225       0.035852       0.005499
  0.228739       0.406233       0.093894       0.079027       0.000000       0.330721      85.651746       0.216429       0.135268       0.002459      13.994097       0.450074      85.582882       0.570585       0.300049       0.000000      13.546483       1.082145       0.356616       5.627152       0.093458      9
3.922774       0.115854       0.016337       0.954255       5.592720       0.127890      93.325135       0.062802       0.022356
1       0.985600          5  556568.000000   20985.400000   41554.600000       3.770501       7.466222    7066.800000       0.000000       0.000000       0.000000    1065.200000    6973.400000      95.600000     108.600000       0.000000     954.400000       0.276261       0.000000       0.000000       0.000000
  0.276261       0.371832       0.084573       0.066564       0.000000       0.403150      86.901131       0.000000       0.000000       0.000000      13.098869       0.000000      85.752582       1.175603       1.335465       0.000000      11.736350       2.796360       0.000000       5.720610       0.000000      9
4.279390       0.000000       0.000000       0.951795       3.876045       1.844565      93.327595       0.492126       0.322442
2       0.985600          5  556568.000000   12589.600000   14185.800000       2.262006       2.548799    6981.000000       0.000000       0.000000       0.000000    1151.000000    6986.600000     140.800000     184.000000       0.000000     820.600000       0.350629       0.000000       0.000000       0.000000       0.350629       0.291340       0.108047       0.102147       0.000000       0.274916      85.846040       0.000000       0.000000       0.000000      14.153960       0.000000      85.914904       1.731431       2.262666       0.000000      10.090999       4.272012       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.762420       2.211018       3.509592      93.516970       0.697118       0.613500
3       0.985600          5  556568.000000    8924.600000    8917.600000       1.603506       1.602248    6909.800000       0.000000       0.000000       0.000000    1222.200000    6957.800000     151.000000     175.400000       0.000000     847.800000       0.399930       0.000000       0.000000       0.000000       0.399930       0.284752       0.147821       0.074089       0.000000       0.243066      84.970487       0.000000       0.000000       0.000000      15.029513       0.000000      85.560748       1.856862       2.156911       0.000000      10.425480       4.294147       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.688637       2.115101       3.605509      93.590753       0.723593       0.630267
4       0.985600          5  556568.000000    6978.800000    6777.000000       1.253899       1.217641    6859.600000       0.000000       0.000000       0.000000    1272.400000    6928.000000     157.800000     172.400000       0.000000     873.800000       0.294951       0.000000       0.000000       0.000000       0.294951       0.305948       0.129090       0.107062       0.000000       0.244462      84.353173       0.000000       0.000000       0.000000      15.646827       0.000000      85.194294       1.940482       2.120020       0.000000      10.745204       4.340876       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.696016       2.075750       3.644860      93.583374       0.723633       0.637145
5       0.985600          5  556568.000000    5790.000000    5582.600000       1.040304       1.003040    6802.000000       0.000000       0.000000       0.000000    1330.000000    6904.800000     160.800000     163.200000       0.000000     903.200000       0.253063       0.000000       0.000000       0.000000       0.253063       0.281413       0.155401       0.114831       0.000000       0.207854      83.644860       0.000000       0.000000       0.000000      16.355140       0.000000      84.909001       1.977373       2.006886       0.000000      11.106739       4.294147       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       0.664043       2.090507       3.630103      93.615347       0.732143       0.634566
cutoff 0 Affix  a       0.581590 b       0.810821: N(rules)=      1.788881*N(trainpairs)^0.810821
         Suffix a      -0.054097 b       0.862983: N(rules)=      0.947341*N(trainpairs)^0.862983
cutoff 1 Affix  a      -0.024716 b       0.808521: N(rules)=      0.975587*N(trainpairs)^0.808521
         Suffix a      -1.391083 b       0.856475: N(rules)=      0.248806*N(trainpairs)^0.856475
cutoff 2 Affix  a      -0.755380 b       0.779266: N(rules)=      0.469832*N(trainpairs)^0.779266
         Suffix a      -1.801956 b       0.843674: N(rules)=      0.164976*N(trainpairs)^0.843674
cutoff 3 Affix  a      -0.965364 b       0.760000: N(rules)=      0.380844*N(trainpairs)^0.760000
         Suffix a      -1.870778 b       0.819764: N(rules)=      0.154004*N(trainpairs)^0.819764
cutoff 4 Affix  a      -1.117902 b       0.750784: N(rules)=      0.326965*N(trainpairs)^0.750784
         Suffix a      -1.821582 b       0.794623: N(rules)=      0.161770*N(trainpairs)^0.794623
cutoff 5 Affix  a      -1.252385 b       0.745853: N(rules)=      0.285822*N(trainpairs)^0.745853
         Suffix a      -1.762567 b       0.774137: N(rules)=      0.171604*N(trainpairs)^0.774137

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        yes
cutoff                      2
fraction          9856.000000
iterations                  5
trainlines      556568.000000
rules            14185.800000 (  12589.600000)
rules%               2.548799 (      2.262006)
same%stdev           0.291340
ambi1%stdev          0.108047
ambi2%stdev          0.102147
ambi3%stdev          0.000000
diff%stdev           0.274916
same%               85.914904 (     85.846040)
ambi1%               1.731431 (      0.000000)
ambi2%               2.262666 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               10.090999 (     14.153960)
amb.rules%           4.272012 (      0.000000)
false_amb%           0.762420 (      0.000000)
false_not_amb%       2.211018 (      5.720610)
true_amb%            3.509592 (      0.000000)
true_not_amb%       93.516970 (     94.279390)
precision            0.697118 (      0.000000)
recall               0.613500 (      0.000000)

bests[5].suffixonly == [false]
bests[5].langbase == [elC1]
comp = comp_parms0_off
bests[5].rows == [1]
  R->R     W->R     R->W     W->W

0.208419 0.719204 -0.602475 0.220196 0.054780 0.157617
*/
//iteration:20.1
/*weight ( used): 1.40401105178530815e+04 suffix only: no */
/* number of nodes: 72056, nodes/line: 1.27600495838498307e-01 weight ( used): 1.40401105178530815e+04 blobs 1 lines 564700 * fraction 1.00000000000000000e+00 = 564700 lines*/
        {{
        2.08418641675330557e-01,    7.19204329987296109e-01,    -6.02475173735412373e-01,    2.20196492945528233e-01,    5.47799882952482026e-02,    1.57617016654553105e-01
        }}
    };

static bestParms best_elC2 =
    {
    false,
    "elC2",
    1,
/*
cutoff  fraction  iterations    trainlines    suffixrules     affixrules        suffix%         affix%      s-same          s-ambiguous                                s-different     a-same        a-ambiguous                                 a-different    s-same-stddev% s-amb-stddev%                                 s-diff-stddev% a-same-stddev% a-amb-stddev%                               a-diff-stddev%   s-same%      s-ambiguous%                                    s-different%    s-amb.rules%  a-same%       a-ambiguous%                                a-different%   a-amb.rules%    s_false_amb  s_false_not_amb  s_true_amb  s_true_not_amb    s_precision       s_recall    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall
0       0.985600          5  556568.000000   91660.600000   76515.400000      16.468895      13.747718    6965.200000      17.600000      11.000000       0.200000    1138.000000    6840.800000     104.200000      61.600000       0.000000    1125.400000       0.210385       0.033225       0.035852       0.005499       0.228739       0.411576       0.061115       0.108814       0.000000       0.398054      85.651746       0.216429       0.135268       0.002459      13.994097       0.450074      84.121987       1.281358       0.757501       0.000000      13.839154       2.528283       0.356616       5.627152       0.093458      93.922774       0.115854       0.016337       2.400394       5.592720       0.127890      91.878997       0.025948       0.022356
1       0.985600          5  556568.000000   20985.400000   43830.800000       3.770501       7.875192    7066.800000       0.000000       0.000000       0.000000    1065.200000    6864.600000     158.400000     145.000000       0.000000     964.000000       0.276261       0.000000       0.000000       0.000000       0.276261       0.562490       0.123768       0.139939       0.000000       0.617860      86.901131       0.000000       0.000000       0.000000      13.098869       0.000000      84.414658       1.947860       1.783079       0.000000      11.854402       4.289228       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       2.365962       3.797344       1.923266      91.913428       0.288987       0.336199
2       0.985600          5  556568.000000   12589.600000   16235.200000       2.262006       2.917020    6981.000000       0.000000       0.000000       0.000000    1151.000000    6899.600000     180.800000     198.800000       0.000000     852.800000       0.350629       0.000000       0.000000       0.000000       0.350629       0.420769       0.143776       0.088419       0.000000       0.344582      85.846040       0.000000       0.000000       0.000000      14.153960       0.000000      84.845057       2.223315       2.444663       0.000000      10.486965       5.241023       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       2.090507       2.570093       3.150516      92.188883       0.429722       0.550731
3       0.985600          5  556568.000000    8924.600000   10638.200000       1.603506       1.911393    6909.800000       0.000000       0.000000       0.000000    1222.200000    6866.600000     197.400000     194.800000       0.000000     873.200000       0.399930       0.000000       0.000000       0.000000       0.399930       0.430627       0.105998       0.103398       0.000000       0.411117      84.970487       0.000000       0.000000       0.000000      15.029513       0.000000      84.439252       2.427447       2.395475       0.000000      10.737826       5.373832       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       2.090507       2.437285       3.283325      92.188883       0.439868       0.573947
4       0.985600          5  556568.000000    6978.800000    8064.000000       1.253899       1.448880    6859.600000       0.000000       0.000000       0.000000    1272.400000    6839.000000     200.400000     188.800000       0.000000     903.800000       0.294951       0.000000       0.000000       0.000000       0.294951       0.379720       0.161416       0.069236       0.000000       0.288051      84.353173       0.000000       0.000000       0.000000      15.646827       0.000000      84.099852       2.464338       2.321692       0.000000      11.114117       5.371372       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       2.058534       2.407772       3.312838      92.220856       0.445879       0.579106
5       0.985600          5  556568.000000    5790.000000    6704.400000       1.040304       1.204597    6802.000000       0.000000       0.000000       0.000000    1330.000000    6807.000000     202.600000     184.600000       0.000000     937.800000       0.253063       0.000000       0.000000       0.000000       0.253063       0.315079       0.154229       0.063064       0.000000       0.202324      83.644860       0.000000       0.000000       0.000000      16.355140       0.000000      83.706345       2.491392       2.270044       0.000000      11.532218       5.388588       0.000000       5.720610       0.000000      94.279390       0.000000       0.000000       2.031481       2.363502       3.357108      92.247909       0.452436       0.586844
cutoff 0 Affix  a       0.721174 b       0.797865: N(rules)=      2.056846*N(trainpairs)^0.797865
         Suffix a      -0.054097 b       0.862983: N(rules)=      0.947341*N(trainpairs)^0.862983
cutoff 1 Affix  a       0.108771 b       0.802291: N(rules)=      1.114907*N(trainpairs)^0.802291
         Suffix a      -1.391083 b       0.856475: N(rules)=      0.248806*N(trainpairs)^0.856475
cutoff 2 Affix  a      -0.518561 b       0.771959: N(rules)=      0.595377*N(trainpairs)^0.771959
         Suffix a      -1.801956 b       0.843674: N(rules)=      0.164976*N(trainpairs)^0.843674
cutoff 3 Affix  a      -0.784815 b       0.760477: N(rules)=      0.456204*N(trainpairs)^0.760477
         Suffix a      -1.870778 b       0.819764: N(rules)=      0.154004*N(trainpairs)^0.819764
cutoff 4 Affix  a      -0.999690 b       0.756766: N(rules)=      0.367993*N(trainpairs)^0.756766
         Suffix a      -1.821582 b       0.794623: N(rules)=      0.161770*N(trainpairs)^0.794623
cutoff 5 Affix  a      -1.143049 b       0.753603: N(rules)=      0.318845*N(trainpairs)^0.753603
         Suffix a      -1.762567 b       0.774137: N(rules)=      0.171604*N(trainpairs)^0.774137

New (old) algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only          no
Redo training        yes
cutoff                      2
fraction          9856.000000
iterations                  5
trainlines      556568.000000
rules            16235.200000 (  12589.600000)
rules%               2.917020 (      2.262006)
same%stdev           0.420769
ambi1%stdev          0.143776
ambi2%stdev          0.088419
ambi3%stdev          0.000000
diff%stdev           0.344582
same%               84.845057 (     85.846040)
ambi1%               2.223315 (      0.000000)
ambi2%               2.444663 (      0.000000)
ambi3%               0.000000 (      0.000000)
diff%               10.486965 (     14.153960)
amb.rules%           5.241023 (      0.000000)
false_amb%           2.090507 (      0.000000)
false_not_amb%       2.570093 (      5.720610)
true_amb%            3.150516 (      0.000000)
true_not_amb%       92.188883 (     94.279390)
precision            0.429722 (      0.000000)
recall               0.550731 (      0.000000)

bests[6].suffixonly == [false]
bests[6].langbase == [elC2]
comp = comp_parms0_off
bests[6].rows == [1]
  R->R     W->R     R->W     W->W

0.081085 0.737793 -0.622803 0.225382 0.019282 0.100173
*/
//iteration:20.4
/*weight ( used): 3.47608497189166446e+04 suffix only: no */
/* number of nodes: 77922, nodes/line: 1.37988312378253930e-01 weight ( used): 3.47608497189166446e+04 blobs 1 lines 564700 * fraction 1.00000000000000000e+00 = 564700 lines*/
        {{
        8.10850759088081185e-02,        7.37792532174295523e-01,        -6.22803219705343025e-01,       2.25381946520170823e-01,        1.92821878529992295e-02,        1.00173423979330933e-01
        }}
    };


static struct bestParms bests[] =    
    {
     best_da3
    ,best_daC0
    ,best_daC1
    ,best_deC0
    ,best_deC1
    ,best_elC0
    ,best_elC1
    ,best_elC2
    ,best_en4
    ,best_en4_suffix
    ,best_en4W
    ,best_en6W
    ,best_en6WS1
    ,best_en6WS2
    ,best_en6WS3
    ,best_enS2
    ,best_enC0
    ,best_isC0
    ,best_nlC0
    ,best_nlC1
    ,best_nlC2
    ,best_is
    ,best_is_suffix
    };

static struct rotation best    =
    /* R_R   W_R   R_W   W_W  R_NA  W_NA */   
    {{ 0.0,  3.0, -3.0,  1.0,  0.0,  0.0 }};

//iteration:19.63
/* 387 381.726058 */
/*
0.000000,    10.000000,    -10.000000,    2.000000,
-13.581532,    -5.623026,    -4.283725,    6.696503,
-2.889162,    1.349916,    0.464059,    -4.429286,
-0.305579,    0.715326,    0.815889,    0.502817
*/
//iteration:19.63
/* 391 387.069724 */
/*
0.000000,    6.000000,    -6.000000,    2.000000,
-6.337474,    0.399136,    -1.608900,    -6.024108,
-3.947893,    0.883151,    2.100111,    3.650880,
0.463499,    1.732265,    1.476537,    -0.767184
*/

// Orthogonalised:

//iteration:19.39
/* number of nodes: 386, nodes/line: 0.387550 weight: 381.768147 blobs 1 lines 1254 * fraction 0.794328 = 996 lines*/
/*
        {                                                 // # decisions
        0.000000,    4.100000,    -6.000000,    1.000000, //9886
        -4.167700,    0.043205,    -0.140436,    -1.019756, //784
        0.040775,    0.183468,    0.096672,    -0.172187, //17
        0.001182,    -0.003050,    -0.002846,    -0.004568  //0
        }                                                 //(0 unresolved comparisons)
*/

static void plus(double * dest, double * term,int cols)
    {
    for(int col = 0;col < cols;++col)
        dest[col] += term[col];
    }


static int improvements = 0;
static double previous[6] = {0,0,0,0,0,0};

void betterfound(int Nnodes,double weight,int swath,int iterations,int blobs,int lines,double fraction,int fraclines,bool improvement,optionStruct * options)
    {
    if(improvement)
        {
        ++improvements;
        FILE * f = fopen(options->currentParms(),"a");
        ++openfiles;
        assert(f);
        fprintf(f,"//-> IMPROVEMENT #%d\n",improvements);
        --openfiles;
        fclose(f);
        }
    else
        {/*
        times(previous,-1.0);
        plus(parms.Matrix,previous,parms.ROWPARMS);
        normalise(parms.Matrix);*/
        }
    best = parms;
    parms.better(options);
    options->setSwath(swath);
    options->setSwathIteration(iterations);
    options->setNumberOfNodes(Nnodes);
    options->setTrainingPairsLines(lines);
    options->setWeight(weight);

    options->printArgFile();
    printf("%d.%d %d  \tparms ",swath,iterations,Nnodes);
    int i = 0;
    for(;i < NPARMS;++i)
        {
        printf("%7.6f",parms.Matrix[i]);
        if(((i+1) % parms.ROWPARMS) == 0)
            {
            if(i == NPARMS - 1)
                printf("\n");
            else
                printf(";\n");
            }
        else
            printf(";");
        }

    FILE * f = fopen(options->bestParms(),"a");
    ++openfiles;
    if(f)
        {
        fprintf(f,"//iteration:%d.%d\n",swath,iterations);
        fprintf(f
               ,"/*weight (%s used): %.*e suffix only: %s */\n"
               , (options->getWeightFunction() == esupport || options->getWeightFunction() == eentropy) ? "" : "not "
               ,DBL_DIG+2
               ,weight
               ,options->suffixOnly() ? "yes" : "no"
               );

        fprintf(f
               ,"/* number of nodes: %d, nodes/line: %.*e weight (%s used): %.*e blobs %d lines %d * fraction %.*e = %d lines*/\n"
               ,Nnodes
               ,DBL_DIG+2,(double)Nnodes/(double)fraclines
               , (options->getWeightFunction() == esupport || options->getWeightFunction() == eentropy) ? "" : "not "
               ,DBL_DIG+2,weight
               ,blobs
               ,lines
               ,DBL_DIG+2,fraction
               ,fraclines
               );
        fprintf(f,"        {{\n        ");
        int i = 0;
        for(;i < NPARMS;++i)
            {
            fprintf(f,"%.*e", DBL_DIG+2,parms.Matrix[i]);
            if(((i+1) % parms.ROWPARMS) == 0)
                {
                if(i == NPARMS - 1)
                    fprintf(f,"\n        ");
                else
                    fprintf(f,",\n        ");
                }
            else
                fprintf(f,",\t");
            }
        fprintf(f,"}}\n\n");
        --openfiles;
        fclose(f);
        }
    }

void worsefound()
    {
    parms = best;
    }

void copybest()
    {
    parms = best; // go on with best result so far.
    }


static double integral(double x,int n)
        {
        if(n > 1)
            return (-1.0/(double)n) * pow(sin(x),n-1) * cos(x) + ((n-1)/(double)n) * integral(x,n-2);
        else if(n == 1)
            return 1 - cos(x);
        else
            return x;
        }

static double maxintegral;
static double minintegral;

#define pi 3.14159265358979323846

void setMinMaxIntegral(int dimensions)
        {
        assert(dimensions >= 2); // dimension of embedding space. So dimensions == 2 is for a circle.
        maxintegral = integral(pi,dimensions - 2);
        minintegral = integral(0.0,dimensions - 2);
        }

static double expectationValueA(double angle,int dimensions)
        {
        double ret;
        if(dimensions == 3)
            ret = cos(angle);
        else if(dimensions == 2)
            {
            ret = 1.0 - (angle / (pi*0.5));
            }
        else
            {
            double inp = integral(angle,dimensions - 2);
            ret = 1.0 - 2.0*(inp - minintegral)/(maxintegral - minintegral); 
            }
        return ret;
        }

static double angle(double expectectationvalue,int dim)
        {
        double mina = 0.0;
        double maxa = pi;
        double mine = expectationValueA(mina,dim);
        double maxe = expectationValueA(maxa,dim);
        if(expectectationvalue == mine)
            return mina;
        else if(expectectationvalue == maxe)
            return maxa;
        else
            {
            double h = 0.0;
            int i;
            for(i = 0;i < 100 && maxa - mina > 0.0000000001;++i)
                {
                h = (mina+maxa)/2.0;
                double e = expectationValueA(h,dim);
                if(e > expectectationvalue)
                    {
                    mine = e;
                    mina = h;
                    }
                else if(e < expectectationvalue)
                    {
                    maxe = e;
                    maxa = h;
                    }
                else
                    break;
                }
            return h;
            }
        }

void testAngle()
    {
    for(int d = 2;d <= 4;++d)
        {
        setMinMaxIntegral(d);
        for(double z = 1.0;z >= -1.0;z -= 0.25)
            {
            double ang = angle(z,d);
            printf("%d %2.3f -> %2.3f: %2.2f\n",d,z,cos(ang),180.0/pi*ang);
            }
        }
    getchar();
    }

static struct rotation goodParms[] = 
    {
        {{ 0.0238310217,-0.6213812731, 0.7764218951, 0.0162447104, 0.0240052761, 0.0982155421},6,"bg","S", 9419,1,0.19208886},
        {{ 0.013192    ,-0.681737    , 0.722663    ,-0.108477    , 0.015373    , 0.028548    },6,"bg","C", 9407,1,0.19651404},
        {{-0.029092    ,-0.729701    , 0.679400    ,-0.062773    ,-0.005901    , 0.033617    },6,"bg","E", 9614,1,0.19904272},

        {{ 0.0251096524,-0.6115198550, 0.7900264600,-0.0240712845, 0.0258880462,-0.0046418833},6,"cs","S",10956,1,0.17227674},
        {{ 0.014974    ,-0.603475    , 0.774494    , 0.070238    , 0.018401    , 0.174589    },6,"cs","C",10825,1,0.17314619},
        {{ 0.119333    ,-0.433123    , 0.815527    , 0.168871    , 0.139466    , 0.291738    },6,"cs","E",11024,1,0.17314619},

        {{-0.002807    ,-0.748938    , 0.628884    ,-0.204355    ,-0.000636    , 0.042753    },6,"da","C",53832,1,0.07773620},
        {{-0.0010481855,-0.7067872619, 0.6965400975,-0.1214246708, 0.0000826859, 0.0232078163},6,"da","S",53558,1,0.07834370},
        {{-0.010943    ,-0.756421    , 0.639771    ,-0.106260    ,-0.004888    , 0.084177    },6,"da","E",54781,1,0.07921851},

        {{ 0.072141    ,-0.562967    , 0.815173    ,-0.043071    , 0.073552    , 0.078053    },6,"de","C",33033,0,0.08559601},
        {{ 0.0172003300,-0.6998264583, 0.7114994769,-0.0376781049, 0.0169872277, 0.0448036365},6,"de","S",32607,0,0.08580908},
        {{-0.006618    ,-0.584949    , 0.803968    ,-0.026080    , 0.000991    , 0.103662    },6,"de","E",33638,0,0.08666139},
        {{ 0.0223945186,-0.5571704335, 0.8262356319, 0.0296707768, 0.0226460927, 0.0707182696},6,"de","D",33279,0,0.08711798},

        //{{ 0.005600    ,-0.438006    , 0.891624    , 0.072352    , 0.006273    , 0.088605    },6}, // Greek     (-XC)
        {{-0.0029958277,-0.4289888183, 0.8632170809, 0.1740615451,-0.0002318057, 0.2012918646},6,"el","E",62501,1,0.09432719},
        {{ 0.039035    ,-0.441248    , 0.889026    , 0.066674    , 0.039309    , 0.086126    },6,"el","C",59706,1,0.09850769},

        {{-0.0227055905,-0.6269595867, 0.7208545657, 0.0364387780,-0.0049863206, 0.2922707307},6,"en","E",11283,1,0.12178388},
        {{ 0.0132135437,-0.6752102320, 0.7316886509,-0.0246358080, 0.0132845242, 0.0881184515},6,"en","S",11163,2,0.12191093},
        {{ 0.1186315828,-0.5514749040, 0.7814237869, 0.1016008409, 0.1187735190, 0.2162150858},6,"en","C",11149,2,0.12267327},

        {{ 0.012820    ,-0.706711    , 0.698674    ,-0.034104    , 0.011938    , 0.104611    },6,"et","C",20219,0,0.18455896},
        {{ 0.0043017280,-0.6706417662, 0.7359937667, 0.0171824972, 0.0049466980, 0.0906343418},6,"et","S",20317,0,0.18636166},
        {{-0.033463    ,-0.750747    , 0.652380    ,-0.072970    ,-0.006969    , 0.065474    },6,"et","E",20568,0,0.18648598},

        {{-0.001161    ,-0.739382    , 0.661059    ,-0.127638    ,-0.000914    ,-0.004777    },6,"fr","C",32620,0,0.09453368},
        {{ 0.051583    ,-0.564670    , 0.769960    , 0.131552    , 0.059553    , 0.254550    },6,"fr","E",33159,0,0.09512586},
        {{ 0.0661175662,-0.5948519753, 0.7849882938, 0.0502248828, 0.0659135351, 0.1367692967},6,"fr","S",33085,1,0.09567646},

        {{-0.032724    ,-0.547237    , 0.728254    , 0.212133    ,-0.015476    , 0.351947    },6,"hu","E",9096,0,0.11809593},
        {{ 0.0334409657,-0.6849089134, 0.7240905990,-0.0096012049, 0.0336626711, 0.0651838026},6,"hu","S",8757,0,0.11904070},
        {{ 0.042048    ,-0.682636    , 0.724419    ,-0.001560    , 0.043343    , 0.074674    },6,"hu","C",8794,0,0.12013081},
        {{ 0.175030    ,-0.437740    , 0.851017    , 0.105433    , 0.175194    , 0.108212    },6,"hu","D",9004,1,0.12165698},

 //       {{ 0.139452    ,-0.511179    , 0.830374    , 0.027272    , 0.139112    , 0.098139    },6}, // Dutch     (-XC) (better)
        {{ 0.024488    ,-0.495293    , 0.826486    , 0.149246    , 0.024359    , 0.219404    },6,"nl","D",40924,1,0.10994085},
        {{ 0.1069193225,-0.6193738021, 0.7649522686,-0.0790680200, 0.1081691280, 0.0428959952},6,"nl","S",40674,1,0.11003574},
        {{ 0.016243    ,-0.614317    , 0.760298    , 0.018110    , 0.021227    , 0.208614    },6,"nl","C",41089,1,0.11060505},
        {{-0.011390    ,-0.594109    , 0.774997    , 0.030348    ,-0.001950    , 0.212978    },6,"nl","E",41185,1,0.11088971},

        {{-0.021046    ,-0.659149    , 0.741112    ,-0.011371    ,-0.005483    , 0.125192    },6,"ro","E",50221,1,0.10066098},
        {{ 0.129777    ,-0.588471    , 0.779661    , 0.038950    , 0.129759    , 0.103123    },6,"ro","C",46469,1,0.11261428},
        {{-0.0016509028,-0.7015134891, 0.7114569708,-0.0355891928,-0.0009517376, 0.0209184461},6,"ro","S",46385,1,0.11366444},

        {{ 0.021767    ,-0.480966    , 0.778614    , 0.266345    , 0.021757    , 0.300908    },6,"sk","C",70835,0,0.06273769},
        {{ 0.004002    ,-0.540399    , 0.809001    , 0.131287    , 0.005075    , 0.190287    },6,"sk","E",74039,1,0.06525689},

        {{ 0.016221    ,-0.408998    , 0.765411    , 0.330531    , 0.016203    , 0.370256    },6,"sl","C",30333,0,0.13036071},
        {{ 0.0242105152,-0.6728127440, 0.7368843437, 0.0312092389, 0.0243645571, 0.0465905375},6,"sl","S",30425,0,0.13132361},
        {{-0.005696    ,-0.563334    , 0.773168    , 0.190243    , 0.002000    , 0.220540    },6,"sl","E",31467,1,0.13487890},

        {{ 0.008000    ,-0.238492    , 0.876932    , 0.278726    , 0.008168    , 0.310311    },6,"uk","C",32925,1,0.14237366},
        {{ 0.038750    ,-0.441348    , 0.888120    , 0.065188    , 0.043326    , 0.093945    },6,"uk","E",33696,1,0.14296821},
    };
static size_t goodParmsIndex = 0;
static double InitialDelta = 0.5;

void brown()
    {
    static int it = 0;
    if(it++ < 2)
        {
        normalise(parms.Matrix);
        return;
        }

    if(goodParmsIndex < sizeof(goodParms)/sizeof(goodParms[0]))
        {
        parms = goodParms[goodParmsIndex++];
        normalise(parms.Matrix);
        return;
        }

    double vector[6];
    size_t i;
    double tangens = InitialDelta*pow(0.995,it);
    /*
    Compute a randomly directed vector in 6 dimensional space.
    */
#if ZIGGURAT
    for(i = 0;i < sizeof(vector)/sizeof(vector[0]);++i)
        {
        vector[i] = RNOR; // Use Ziggurat algorithm (rnorrexp.c)
        }
#else
    double radius2 = 0.0;
    do
        {
        int i;
        radius2 = 0.0;
        for(i = 0;i < parms.ROWPARMS;++i)
            {
            vector[i] = rand() - (RAND_MAX/2);
            radius2 += vector[i]*vector[i]; 
            }
        }
    while(radius2 > ((double)RAND_MAX/2.0)*((double)RAND_MAX/2.0));
#endif
    normalise(vector);
    double inproduct = inner(parms.Matrix,vector); // Only first row. Ignore the rest.
    struct rotation diff = parms;
    times(diff.Matrix,-inproduct);
    plus(vector,diff.Matrix,parms.ROWPARMS);
    normalise(vector);

    /*
    We require that the delta is not only perpendicular to the currently best
    vector, but also to the previous delta. This could be extended to comprise
    even earlier delta's, up to three. (Though there are six values to define,
    there are only five degrees of freedom, because we look at a vector
    sweeping the unit sphere in six dimensional space. So we can require delta
    to be perpendicular to the current vector and to three previous delta's
    and still have one degree of freedom.)
    */
    inproduct = inner(previous,vector);
    times(previous,-inproduct);
    plus(vector,previous,parms.ROWPARMS);
    normalise(vector);

    times(vector,tangens);

    for(i = 0;i < sizeof(vector)/sizeof(vector[0]);++i)
        {
        previous[i] = vector[i];
        }

    plus(parms.Matrix,vector,parms.ROWPARMS);
    normalise(parms.Matrix);
    }

bool init(optionStruct * options)
    {
#if ZIGGURAT
    zigset(86947731);
#endif
    parms.init(options);
    double MinMaxInnerProduct = 1.0;
/*
    double MinInnerProduct = 1.0;
    size_t outlier = 1000000;
    size_t neighbour = 1000000;
    size_t furthest1 = 1000000;
    size_t furthest2 = 1000000;
*/
    for(size_t i = 0;i < sizeof(goodParms)/sizeof(goodParms[0]);++i)
        {
        double MaxInnerProduct = -1.0;
        for(size_t j = 0;j < sizeof(goodParms)/sizeof(goodParms[0]);++j)
            {
            if(i != j)
                {
                double InnerProduct = inner(goodParms[i].Matrix,goodParms[j].Matrix);
                if(InnerProduct > MaxInnerProduct) // Find closest meighbour of i
                    {
                    MaxInnerProduct = InnerProduct;
//                    neighbour = j;
                    }
/*
                if(InnerProduct < MinInnerProduct)
                    {
                    furthest1 = i;
                    furthest2 = j;
                    MinInnerProduct = InnerProduct;
                    }
*/
                }
            }
        if(MaxInnerProduct < MinMaxInnerProduct)
            {
            MinMaxInnerProduct = MaxInnerProduct;
//            outlier = i;
            }
        }
    /* Take the distance between the most outlying vector and its closest
       neighbour as the initial headroom for changing a vector. */
    InitialDelta = sqrt(1.0 - MinMaxInnerProduct * MinMaxInnerProduct); 
/*
    printf("InitialDelta %f outlier %d closest to %d\n",InitialDelta,outlier,neighbour);
    printf("Furthest distance %f between %d and %d\n",sqrt(1.0 - MinInnerProduct * MinInnerProduct),furthest1,furthest2);
*/
    return true;
    }

void printparms(int Nnodes,double weight,optionStruct * options)
    {
    int i;
    FILE * f = fopen(options->currentParms(),"a");
    ++openfiles;
    assert(f);
    fprintf(f
           ,"/*#nodes in tree: %d weight (%s used): %.*e , Suffix only: %s */\n"
           ,Nnodes
           ,   (  options->getWeightFunction() == esupport 
               || options->getWeightFunction() == eentropy
               )                                      ? "more support is better" 
             : options->getWeightFunction() == edepth ? "fewer non-wildcard characters is better" 
             : options->getWeightFunction() == esize  ? "fewer characters is better" 
             :                                          "not used"
           ,DBL_DIG+2
           ,weight
           ,options->suffixOnly() ? "yes" : "no"
           );
    fprintf(f,"        {\n        ");
    for(i = 0;i < NPARMS;++i)
        {
        fprintf(f,"%.*e", DBL_DIG+2,parms.Matrix[i]);
        if(((i+1) % parms.ROWPARMS) == 0)
            {
            if(i == NPARMS - 1)
                fprintf(f,"\n        ");
            else
                fprintf(f,",\n        ");
            }
        else
            fprintf(f,",");
        }
    fprintf(f,"}\n");
    --openfiles;
    fclose(f);
    }

#define RR 0
#define WR 1
#define RW 2
#define WW 3
#define RN 4
#define WN 5

void computeWeight(vertex * a)
    {
    double A = parms.Matrix[RR] * a->R__R + parms.Matrix[WR] * a->W__R + parms.Matrix[RW] * a->R__W + parms.Matrix[WW] * a->W__W;
#if _NA
    if (parms.ROWPARMS == 6)
        {
        A += parms.Matrix[RN] * a->R__NA + parms.Matrix[WN] * a->W__NA;
        }
#endif
    a->wght = A;
    }

// Notice that we want to sort from low penalty to high penalty, therefore sgn(a - b)
int comp_parms(const vertex * a,const vertex * b)
    {
    if (a->wght != b->wght)
        return a->wght > b->wght ? -1 : 1;
    return 0;
    }

static int nparms = 0;

static int comp_parms0_off(const vertex * a,const vertex * b)
    {   
    double A = parms.Matrix[RR]*a->R__R + parms.Matrix[WR]*a->W__R + parms.Matrix[RW]*a->R__W + parms.Matrix[WW]*a->W__W;
    double B = parms.Matrix[RR]*b->R__R + parms.Matrix[WR]*b->W__R + parms.Matrix[RW]*b->R__W + parms.Matrix[WW]*b->W__W;
#if _NA
    if(parms.ROWPARMS == 6)
        {
        A += parms.Matrix[RN]*a->R__NA + parms.Matrix[WN]*a->W__NA;
        B += parms.Matrix[RN]*b->R__NA + parms.Matrix[WN]*b->W__NA;
        }
#endif
    if(A != B)
        {
        return A > B ? -1 : 1;
        }
    return 0;
    }

struct funcstruct
    {
    const char * number;
    const char * name;
    int (*comp)(const vertex * a,const vertex * b);
    };

static struct funcstruct funcstructs[] =
    {
#if _NA
        {"1","fairly_good",comp_fairly_good},
        {"2","even_better",comp_even_better},
        {"3","affiksFEW3",comp_affiksFEW3},
        {"4","affiksFEW",comp_affiksFEW},
        {"5","affiksFEW2",comp_affiksFEW2},
        {"6","fixNA",comp_fixNA},
        {"7","fruit",comp_fruit},
        {"8","ice",comp_ice},
        {"9","pisang",comp_pisang},
        {"10","kiwi",comp_kiwi},
        {"11","carrot",comp_carrot},
        {"12","peen",comp_peen},
        {"13","beet",comp_beet},
        {"14","sugar",comp_sugar},
        {"15","affiksFEW2org",comp_affiksFEW2org},
#endif
        {"16","honey",comp_honey},
        {"17","koud",comp_koud},
        {"18","parms0",comp_parms0_off},
        {"19","parmsoff",comp_parms0_off},
        {0,0,0}
    };


void setCompetitionFunction(optionStruct * options)
    {
    size_t langlength = strlen(options->extra());
    const char * underscore = strchr(options->extra(),'_');
    if(underscore != NULL)
        langlength = (size_t)(underscore - options->extra());

    if(options->compfunc())
        {
        for(int i = 0;funcstructs[i].number;++i)
            if(!strcmp(options->compfunc(),funcstructs[i].number) || !strcmp(options->compfunc(),funcstructs[i].name))
                {
                comp = funcstructs[i].comp;
                if(comp == comp_parms0_off)
                    {
                    for(unsigned int j = 0;j < sizeof(bests)/sizeof(bests[0]);++j)
                        {
                        if(  bests[j].suffixonly == options->suffixOnly() 
                          && (langlength == strlen(bests[j].langbase)) // 20130125
                          && !strncmp(bests[j].langbase,options->extra(),strlen(bests[j].langbase))
                          )
                            {
                            printf("bests[%d].suffixonly == [%s] bests[%d].langbase == [%s]\n",j,bests[j].suffixonly ? "true" : "false",j,bests[j].langbase);
                            printf("comp = comp_parms0_off\n");
                            comp = comp_parms0_off;
                            printf("bests[%d].rows == [%d]\n",j,bests[j].rowss);
                            nparms = bests[j].rowss * parms.ROWPARMS;
                            if(nparms > NPARMS)
                                {
                                fprintf(stderr,"Too many rows of parameters in bestParms struct for %s (%d, max allowed %d)\n",options->extra(),nparms,NPARMS);
                                exit(-1);
                                }
                            parms = bests[j].val;
                            if(options->currentParms())
                                {
                                FILE * f = fopen(options->currentParms(),"w");
                                ++openfiles;
                                if(f)
                                    {
                                    fprintf(f,"bests[%d].suffixonly == [%s]\nbests[%d].langbase == [%s]\n",j,bests[j].suffixonly ? "true" : "false",j,bests[j].langbase);
                                    fprintf(f,"comp = comp_parms0_off\n");
                                    fprintf(f,"bests[%d].rows == [%d]\n",j,bests[j].rowss);
                                    fprintf(f,"  R->R     W->R     R->W     W->W\n");
                                    for(int k = 0;k < nparms;++k)
                                        {
                                        if(k % parms.ROWPARMS == 0)
                                            fprintf(f,"\n");
                                        fprintf(f,"%6f ",parms.Matrix[k]);
                                        }
                                    fprintf(f,"\n");
                                    --openfiles;
                                    fclose(f);
                                    }
                                }
                            break;
                            }
                        }
                    if(nparms == 0)
                        {
                        fprintf(stderr,"No parameters defined for \"%s\"\nChoose one of:\n",options->extra());
                        for(unsigned int j = 0;j < sizeof(bests)/sizeof(bests[0]);++j)
                            {
                            fprintf(stderr,"\t%s %s\n",bests[j].langbase,bests[j].suffixonly ? "suffix":"affix");
                            }
                        fprintf(stderr,"Or find optimal parameters for %s and put these in comp.cpp.\n",options->extra());
                        getchar();
                        exit(-1);
                        }
                    if(options->verbose())
                        {
                        printf("comp_parms0_off\n");
                        }
                    }
                return;
                }
        }
    if(options->numberOfParms() == 4 || options->numberOfParms() == 6)
        {
        comp = comp_parms0_off;
        return;
        }
    fprintf(stderr,"Error: Unknown competition function %s\n",options->compfunc());
    getchar();
    exit(2);
    }

