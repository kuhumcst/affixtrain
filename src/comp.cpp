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
#include "graph.h"

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
    // Danish 90.942165	+/- 0.589437 at 0.9856 of dataset, 5 iterations, 32327.400000 = 5.925881% rules, cutoff = 1
    // German 90.266461 +/-	0.509202 at 0.9856 of dataset, 7 iterations, 21539.428571 = 6.930653% rules, cutoff = 1
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
    // German 90.266461	+/-	0.509202 at 0.9856 of dataset, 7 iterations, 21539.428571 = 6.930653% rules, cutoff = 1
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
    //English 87.626771	+/- 0.060148 at 0.4928 (!) of dataset, 3 iterations, 933.000000 = 2.480262% rules
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
    // Slovene 86.133458	+/- 0.549185 at 0.9856 of dataset, 9 iterations, 40898.777778 = 20.904266% rules
    // English 87.803261	+/- 0.106156 at 0.4928 (!) of dataset, 3 iterations, 889.333333 = 2.364179% rules
    // Dutch 89.837692 +/- 0.412795 at 0.9856 of dataset, 7 iterations, 56640.285714 = 19.035104% rules, cutoff = 0
    // ALMOST BEST German 91.288892 +/- 0.670828 at 0.9856 of dataset, 7 iterations, 50584.857143 = 16.276480% rules, cutoff = 0
    // Swedish 91.873698 +/- 0.367967 at 0.9856 of dataset, 6 iterations, 9066.166667 = 1.924725% rules, cutoff = 2
    // ALMOST BEST Norwegian 87.535644 +/- 0.344659 at 0.9856 of dataset, 6 iterations, 48468 = 10.268492% rules, cutoff = 1
    // ALMOST BEST Greek 90.414875+/- 0.385254 at 0.9856 of dataset, 5 iterations, 120691.4 = 22.303999% rules, cutoff = 0
    // BEST Danish 92.796387 +/- 0.214267 at 0.9856 of dataset, 5 iterations, 67807 = 12.429587% rules, cutoff = 0
    // ALMOST BEST Russian 80.484806 +/- 0.409391 at 0.9856 of dataset, 6 iterations, 54630 = 14.022614% rules, cutoff = 1
    // BEST Polish 93.880103 +/- 0.077021 at 0.9856 of dataset, 2 iterations, 344944.5	= 10.165818% rules, cutoff = 0
//AMBI:
    // French ok 84.993631 ambi1 1.388535 ambi2 1.085987 diff 12.531847 rules 7318.375000 2.724738% cutoff 2
/*
0	0.985600	2	3490123.000000	0.000000	415069.500000	0.000000	11.892690	0.000000	0.000000	0.000000	0.000000	0.000000	45811.000000	299.500000	168.000000	0.000000	4714.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.044374	0.026347	0.033280	0.000000	0.015253	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.837821	0.587336	0.329457	0.000000	9.245387	1.282529	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.013865	6.709744	0.268664	92.007727	0.116994	0.038499
1	0.985600	2	3490123.000000	0.000000	198203.500000	0.000000	5.678983	0.000000	0.000000	0.000000	0.000000	0.000000	46176.500000	283.000000	241.000000	0.000000	4292.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.092907	0.063787	0.016640	0.000000	0.012480	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.554586	0.554978	0.472614	0.000000	8.417822	1.257035	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.784421	6.505795	0.472614	92.237170	0.231508	0.067725
2	0.985600	2	3490123.000000	0.000000	57342.000000	0.000000	1.642979	0.000000	0.000000	0.000000	0.000000	0.000000	46504.500000	192.000000	245.000000	0.000000	4051.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.134507	0.030507	0.044374	0.000000	0.059627	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	91.197811	0.376522	0.480458	0.000000	7.945208	0.999157	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.543212	6.522464	0.455945	92.478379	0.295613	0.065337
3	0.985600	2	3490123.000000	0.000000	34167.000000	0.000000	0.978963	0.000000	0.000000	0.000000	0.000000	0.000000	46470.500000	178.500000	210.000000	0.000000	4134.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.051307	0.029120	0.030507	0.000000	0.008320	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	91.131136	0.350048	0.411821	0.000000	8.106995	0.890318	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.498108	6.586198	0.392211	92.523484	0.282486	0.056203
4	0.985600	2	3490123.000000	0.000000	24896.500000	0.000000	0.713342	0.000000	0.000000	0.000000	0.000000	0.000000	46392.500000	166.500000	180.000000	0.000000	4254.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.059627	0.029120	0.044374	0.000000	0.013867	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.978173	0.326515	0.352990	0.000000	8.342321	0.809915	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.455945	6.624439	0.353970	92.565646	0.279628	0.050724
5	0.985600	2	3490123.000000	0.000000	19778.500000	0.000000	0.566699	0.000000	0.000000	0.000000	0.000000	0.000000	46335.500000	151.500000	180.500000	0.000000	4325.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.065174	0.006933	0.009707	0.000000	0.048534	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.866393	0.297100	0.353970	0.000000	8.482537	0.759908	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.429471	6.647971	0.330438	92.592120	0.277824	0.047351

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

int (*comp)(const vertex * a,const vertex * b) = comp_koud;
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
    double Matrix[ROWPARMS*ROWPARMS];
    } rotation;

static struct rotation parms    =
   /* R_R  W_R  R_W  W_W */   
    {{ 0.0, 1.0, 0.0, 0.0
    ,-1.0, 0.0, 0.0, 0.0
    , 0.0, 0.0, 0.0,-1.2
    , 0.0, 0.0,-1.2, 0.0
    }};

static void printmatrix(const char * msg,struct rotation arr)
    {
    printf("%s\n",msg);
    for(int i = 0;i < ROWPARMS;++i)
        {
        for(int j = 0;j < ROWPARMS;++j)
            {
            printf("%f ",arr.Matrix[ROWPARMS * i+j]);
            }
        putchar('\n');
        }
    putchar('\n');
    }

static void normalise(double * ROW)
    {
    double modulus = 0.0;
    for(int i = 0;i < ROWPARMS;++i)
        modulus += ROW[i] * ROW[i];
    modulus = sqrt(modulus);
    for(int i = 0;i < ROWPARMS;++i)
        ROW[i] /= modulus;
    }

static void orthogonalise(const char * Msg)
    {
    if(VERBOSE)
        {
        printf("%s: ",Msg);
        printmatrix("before orthogonalisation",parms);
        }
    for(int i = 1;i < sizeof(parms.Matrix)/(sizeof(parms.Matrix[0]) * ROWPARMS);++i)
        {
        int I = i * ROWPARMS;
        for(int j = 0; j < i;++j)
            {
            int J = j * ROWPARMS;
            double inner = 0.0;
            double Jmodulus = 0.0;
            for(int k = 0;k < ROWPARMS;++k)
                {
                inner += parms.Matrix[I+k]*parms.Matrix[J+k];
                Jmodulus += parms.Matrix[J+k]*parms.Matrix[J+k];
                }
            inner /= Jmodulus;
            for(int k = 0;k < ROWPARMS;++k)
                {
                parms.Matrix[I+k] -= inner * parms.Matrix[J+k];
                }
            }
        }
    for(int row = 0;row < ROWPARMS;++row)
        {
        normalise(parms.Matrix+row*ROWPARMS);
        }
    if(VERBOSE)
        {
        printf("%s: ",Msg);
        printmatrix("after orthogonalisation",parms);
        }
    }


struct rotation makeRotationMatrix(int r, double pone, double delta) // 0 <= r < 6
    {
    assert(-1.0 < delta && delta < 1.0);
    assert(0 <= r && r < 6);
    assert(pone == -1.0 || pone == 1.0);
    double contradelta = sqrt(1.0 - delta*delta);
    struct rotation Rotation;
    for(int i = 0;i < sizeof(Rotation.Matrix)/sizeof(Rotation.Matrix[0]);++i)
        {
        int col = i % ROWPARMS;
        int row = i / ROWPARMS;
        if(col == row)
            { // diagonal. 1 or something between 0 and 1
            switch(r)
                {
            case 0:
                if(col == 0 || col == 1)
                    Rotation.Matrix[i] = contradelta; 
                else
                    Rotation.Matrix[i] = 1.0; 
                break;
            case 1:
                if(col == 1 || col == 2)
                    Rotation.Matrix[i] = contradelta; 
                else
                    Rotation.Matrix[i] = 1.0; 
                break;
            case 2:
                if(col == 2 || col == 3)
                    Rotation.Matrix[i] = contradelta; 
                else
                    Rotation.Matrix[i] = 1.0; 
                break;
            case 3:
                if(col == 0 || col == 2)
                    Rotation.Matrix[i] = contradelta; 
                else
                    Rotation.Matrix[i] = 1.0; 
                break;
            case 4:
                if(col == 1 || col == 3)
                    Rotation.Matrix[i] = contradelta; 
                else
                    Rotation.Matrix[i] = 1.0; 
                break;
            case 5:
                if(col == 0 || col == 3)
                    Rotation.Matrix[i] = contradelta; 
                else
                    Rotation.Matrix[i] = 1.0; 
                break;
                }
            }
        else
            {
            switch(r)
                {
            case 0:
                if(row == 0 && col == 1)
                    Rotation.Matrix[i] = pone*delta; 
                else if(row == 1 && col == 0)
                    Rotation.Matrix[i] = -pone*delta;
                else
                    Rotation.Matrix[i] = 0.0; 
                break;
            case 1:
                if(row == 1 && col == 2)
                    Rotation.Matrix[i] = pone*delta; 
                else if(row == 2 && col == 1)
                    Rotation.Matrix[i] = -pone*delta;
                else
                    Rotation.Matrix[i] = 0.0; 
                break;
            case 2:
                if(row == 2 && col == 3)
                    Rotation.Matrix[i] = pone*delta; 
                else if(row == 3 && col == 2)
                    Rotation.Matrix[i] = -pone*delta;
                else
                    Rotation.Matrix[i] = 0.0; 
                break;
            case 3:
                if(row == 0 && col == 2)
                    Rotation.Matrix[i] = pone*delta; 
                else if(row == 2 && col == 0)
                    Rotation.Matrix[i] = -pone*delta;
                else
                    Rotation.Matrix[i] = 0.0; 
                break;
            case 4:
                if(row == 1 && col == 3)
                    Rotation.Matrix[i] = pone*delta; 
                else if(row == 3 && col == 1)
                    Rotation.Matrix[i] = -pone*delta;
                else
                    Rotation.Matrix[i] = 0.0; 
                break;
            case 5:
                if(row == 0 && col == 3)
                    Rotation.Matrix[i] = pone*delta; 
                else if(row == 3 && col == 0)
                    Rotation.Matrix[i] = -pone*delta;
                else
                    Rotation.Matrix[i] = 0.0; 
                break;
                }
            }
        }
    return Rotation;
    }

struct rotation rotate(struct rotation parms,struct rotation Rotation)
    {
    struct rotation newparms;
    for(int row = 0;row < ROWPARMS;++row)
        {
        for(int col = 0;col < ROWPARMS;++col)
            {
            newparms.Matrix[row * ROWPARMS + col] = 0.0;
            for(int i = 0;i < ROWPARMS;++i)
                newparms.Matrix[row * ROWPARMS + col] += parms.Matrix[row * ROWPARMS + i] * Rotation.Matrix[i * ROWPARMS + col];
            }
        }
    return newparms;
    }

struct bestParms
    {
    bool suffixonly;
    const char * langbase;
    int rows;
    double val[NPARMS];
    // Each row:
    // R__R W__R R__W W__W
    // Generally, good that Wrongs change to Rights (W__R > 0) and that Rights don't change to Wrongs (R__W < 0)
    // But what about rules that don't improve lemmatisation? (R__R > 0 or W__W > 0)
    // Intuitively difficult to decide!
    };

static bestParms  best_pl = 
    {
    false,
    "pl",
    3,
//iteration:15.14 count
/* 96974 90763.946221 */
/*
0	0.985600	2	3490123.000000	0.000000	289887.000000	0.000000	8.305925	0.000000	0.000000	0.000000	0.000000	0.000000	45478.500000	484.000000	292.500000	0.000000	4738.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.018027	0.069334	0.001387	0.000000	0.085974	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.185771	0.949150	0.573608	0.000000	9.291471	2.009099	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.727688	6.696998	0.281411	91.293903	0.075308	0.040326
1	0.985600	2	3490123.000000	0.000000	154622.000000	0.000000	4.430274	0.000000	0.000000	0.000000	0.000000	0.000000	45870.000000	495.500000	434.500000	0.000000	4193.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.019413	0.081814	0.029120	0.000000	0.033280	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.953523	0.971702	0.852078	0.000000	8.222697	2.234620	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.482556	6.226345	0.752064	91.539035	0.202321	0.107770
2	0.985600	2	3490123.000000	0.000000	45487.000000	0.000000	1.303307	0.000000	0.000000	0.000000	0.000000	0.000000	46324.500000	385.500000	337.500000	0.000000	3945.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.040214	0.056854	0.037440	0.000000	0.020800	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.844822	0.755986	0.661856	0.000000	7.737336	1.649246	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.040339	6.369502	0.608907	91.981252	0.226394	0.087256
3	0.985600	2	3490123.000000	0.000000	27160.500000	0.000000	0.778210	0.000000	0.000000	0.000000	0.000000	0.000000	46347.500000	326.500000	295.500000	0.000000	4023.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.070720	0.040214	0.037440	0.000000	0.067947	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.889926	0.640284	0.579491	0.000000	7.890299	1.432550	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.872669	6.418528	0.559881	92.148922	0.242875	0.080230
4	0.985600	2	3490123.000000	0.000000	20285.000000	0.000000	0.581212	0.000000	0.000000	0.000000	0.000000	0.000000	46308.000000	289.000000	254.000000	0.000000	4142.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.016640	0.013867	0.036054	0.000000	0.038827	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.812464	0.566744	0.498108	0.000000	8.122684	1.262918	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.741278	6.456769	0.521640	92.280313	0.260274	0.074751
5	0.985600	2	3490123.000000	0.000000	16391.500000	0.000000	0.469654	0.000000	0.000000	0.000000	0.000000	0.000000	46242.500000	273.000000	244.500000	0.000000	4233.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.009707	0.019413	0.045760	0.000000	0.036054	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.684015	0.535368	0.479478	0.000000	8.301139	1.196243	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.709901	6.492068	0.486341	92.311690	0.255144	0.069692

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      2
trainlines      3490123.000000
rules           45487.000000
rules%         1.303307
same%stdev     0.040214
ambi1%stdev    0.056854
ambi2%stdev    0.037440
ambi3%stdev    0.000000
diff%stdev     0.020800
same%          90.844822
ambi1%         0.755986
ambi2%         0.661856
ambi3%         0.000000
diff%          7.737336
amb.rules%     1.649246
false_amb%     1.040339
false_not_amb% 6.369502
true_amb%      0.608907
true_not_amb%  91.981252
precision       0.226394
recall          0.087256

bests[19].suffixonly == [false]
bests[19].langbase == [pl]
comp = comp_parms0_off
bests[19].rows == [3]
  R->R  W->R  R->W  W->W

     0    55   -76     3
     4   181  -177   113
    68   195  -151    60
*/
        {
            0,   55,  -76,    3, //6146838
            4,  181, -177,  113, //712742
           68,  195, -151,   60  //5628
                                 //        0
        }
//OnlyZeros 12 
//suffix only no 
    };

static bestParms  best_sl = 
    {
    false,
    "sl",
    3,
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          0
fraction        9856.000000
iterations      9
trainlines      204030.000000
rules           29635.333333
rules%         14.524988
same%stdev     0.941088
ambi1%stdev    0.196612
ambi2%stdev    0.591892
ambi3%stdev    0.000000
diff%stdev     1.043102
same%          80.332522
ambi1%         3.612770
ambi2%         3.022930
ambi3%         0.000000
diff%          13.031778
amb.rules%     7.605250
false_amb%     3.962988
false_not_amb% 5.076311
true_amb%      3.642262
true_not_amb%  87.318440
precision       0.314850
recall          0.417759
*/
//OnlyZeros 12 
//suffix only no 
//iteration:19.58
/* 21083 */
        {
            0,  505, -663,   20, //905780
          -77,  856,-1263,  297, //161407
         -768,  640,-1024,  768 //1542
//          512, -512,-1536, 1024  //0
                                 //        0
        }
    };

static bestParms best_sl_suffix =
    {
    true,
    "sl",
    3,
//OnlyZeros 12 
//suffix only yes 
//iteration:19.55
/* 26531 */
        {
          763,  537,  -90, -258, //444000
          540, 1508,-1252,   38, //38754
        -1488, -704,-1792,  272 //5586
        //-3072,    0,    0,-1024  //0
                                 //        0
        }
    };


static bestParms  best_fr = 
    {
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          1
fraction        9856.000000
iterations      8
trainlines      268628.750000
rules           16521.750000
rules%         6.150403
same%stdev     0.405192
ambi1%stdev    0.214864
ambi2%stdev    0.273387
ambi3%stdev    0.000000
diff%stdev     0.627035
same%          85.834674
ambi1%         2.640720
ambi2%         1.949180
ambi3%         0.000000
diff%          9.575426
amb.rules%     5.062721
false_amb%     3.351560
false_not_amb% 6.510132
true_amb%      1.711161
true_not_amb%  88.427147
precision       0.203364
recall          0.208138
*/
    false,
    "fr",
    3,
//OnlyZeros 12 
//suffix only no 
//iteration:19.63
/* 19238 */
        {
            0,  262, -304,   47, //728984
            3,  433, -833,   19, //139658
           16,  288,-1440,  144 //550
        //-1664, -256,-2304,  384  //0
                                 //        0
        }
    };

static bestParms best_fr_suffix =
    {
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          0
fraction        9856.000000
iterations      8
trainlines      268628.750000
rules           39414.250000
rules%         14.672387
same%stdev     0.355492
ambi1%stdev    0.324656
ambi2%stdev    0.255550
ambi3%stdev    0.000000
diff%stdev     0.644838
same%          87.204889
ambi1%         2.132518
ambi2%         1.247990
ambi3%         0.000000
diff%          9.414603
amb.rules%     3.618527
false_amb%     2.196848
false_not_amb% 6.799614
true_amb%      1.421679
true_not_amb%  89.581859
precision       0.244469
recall          0.172926
*/
    true,
    "fr",
    3,
//OnlyZeros 12 
//suffix only yes 
//iteration:19.60
/* 22512 */
        {
            0,  470, -537,   28, //308198
          219,  255,  -43,    2, //74208
           -4,   89, -112,   27 //6879
//         -224,  576,-1664,  480  //0
                                 //        0
        }
    };


static bestParms  best_sv = 
    {
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      6
trainlines      471038.000000
rules           6864.166667
rules%         1.457243
same%stdev     0.374487
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     0.374487
same%          92.256284
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          7.743716
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan
*/
    false,
    "sv",
    3,
        {
        //OnlyZeros 12 
        //iteration:19.62
        /* 33474 */
           0, 691, -793, 120,//971415
           2, 586,-1540, 440,//143569
        -432, 512,-1024,  16 //1224
        //-448,1024, -512,1024 //0
                             //0
        }
    };

static bestParms best_sv_suffix =
    {
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          0
fraction        9856.000000
iterations      6
trainlines      471038.000000
rules           51760.833333
rules%         10.988675
same%stdev     0.400911
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     0.400911
same%          91.917284
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          8.082716
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan
*/
    true,
    "sv",
    3,
        {
        //OnlyZeros 12 
        //iteration:19.58
        /* 41253 */
           0, 608, -608,  1,//285302
          -1, 600, -584, 80,//91856
         -82, 450,-1022,334 //17269
        //-448,-512,-2560,512 //0
                            //0
        }
    };

static bestParms best_is_suffix =
    {
    true,
    "is",
    3,
        {
// BEST Icelandic 75.010529 +/- 1.424540 at 0.9856 of dataset, 17 iterations, 19041.705882 = 33.208416% rules, cutoff = 0
// (old)Icelandic 73.620665 +/- 1.337148 at 0.9856 of dataset, 17 iterations, 18363.235294 = 32.025175% rules, cutoff = 0
/*
     -7,  41, -24,  26,
    -80,   0, -56, -24
*/

//iteration:331
/* 18952 */
// 74.877159 +/- 1.394266 at 0.9856 of dataset, 17 iterations, 18736.588235 = 32.676296% rules, cutoff = 0
    /*
0,       0,       0,       0,//            0
0,       0,       0,       0,//            0
0,       6,       -4,      2,//       229741
-10,     -6,      -11,     1 //        18595
*/
//OnlyZeros 12 
//iteration:19.58
/* 16365 */
// weightedcount: ret = ret + 1000.0/(1000.0 + rcount*rcount); (not sure!)
    -15,269, -174,126,//143775
    416,589, -136, 51,//41148
    -30,528,-2044,215,//5557
    //0,768,-2304,512 //0
    /*
cutoff          0
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           19312.470588
rules%         33.680625
same%stdev     1.273580
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.273580
same%          75.396603  (ACL, old algo: 73.2±1.4)
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          24.603397
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan
*/				    //0
        }
    };

#if 1
static bestParms best_is =
    {
    false,
    "is",
    3,
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


#elif 1
static bestParms best_is =
    {
    false,
    "is",
    3,
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57340.000000	0.000000	17314.882353	0.000000	30.196865	0.000000	0.000000	0.000000	0.000000	0.000000	596.352941	0.000000	0.000000	0.000000	241.647059	0.000000	0.000000	0.000000	0.000000	0.000000	1.502313	0.000000	0.000000	0.000000	1.502313	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.163835	0.000000	0.000000	0.000000	28.836165	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
1	0.985600	17	57340.000000	0.000000	9400.823529	0.000000	16.394879	0.000000	0.000000	0.000000	0.000000	0.000000	595.000000	0.000000	0.000000	0.000000	243.000000	0.000000	0.000000	0.000000	0.000000	0.000000	1.633472	0.000000	0.000000	0.000000	1.633472	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.002387	0.000000	0.000000	0.000000	28.997613	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
2	0.985600	17	57340.000000	0.000000	3618.000000	0.000000	6.309731	0.000000	0.000000	0.000000	0.000000	0.000000	581.529412	0.000000	0.000000	0.000000	256.470588	0.000000	0.000000	0.000000	0.000000	0.000000	1.819133	0.000000	0.000000	0.000000	1.819133	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	69.394918	0.000000	0.000000	0.000000	30.605082	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
3	0.985600	17	57340.000000	0.000000	2230.235294	0.000000	3.889493	0.000000	0.000000	0.000000	0.000000	0.000000	570.058824	0.000000	0.000000	0.000000	267.941176	0.000000	0.000000	0.000000	0.000000	0.000000	1.881354	0.000000	0.000000	0.000000	1.881354	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	68.026113	0.000000	0.000000	0.000000	31.973887	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
4	0.985600	17	57340.000000	0.000000	1616.352941	0.000000	2.818892	0.000000	0.000000	0.000000	0.000000	0.000000	560.235294	0.000000	0.000000	0.000000	277.764706	0.000000	0.000000	0.000000	0.000000	0.000000	2.021841	0.000000	0.000000	0.000000	2.021841	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	66.853854	0.000000	0.000000	0.000000	33.146146	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
5	0.985600	17	57340.000000	0.000000	1242.235294	0.000000	2.166438	0.000000	0.000000	0.000000	0.000000	0.000000	548.764706	0.000000	0.000000	0.000000	289.235294	0.000000	0.000000	0.000000	0.000000	0.000000	2.108042	0.000000	0.000000	0.000000	2.108042	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	65.485048	0.000000	0.000000	0.000000	34.514952	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          0
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           17314.882353
rules%         30.196865
same%stdev     1.502313
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.502313
same%          71.163835
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          28.836165
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan

bests[9].suffixonly == [false]
bests[9].langbase == [is]
comp = comp_parms0_off
bests[9].rows == [3]
  R->R  W->R  R->W  W->W

     2   771  -640   213
   -72   496 -1152   376
     0   512 -1792  -256

*/
//iteration:19.63 count
/* 14458 14190.158755 */
        {
            2,  771, -640,  213, //531302
          -72,  496,-1152,  376, //10780
            0,  512,-1792, -256 //927
          //256,-1280,-2048, 1024  //0
                                 //        0
        }
//OnlyZeros 12 
//suffix only no 
    };

#elif 0
static bestParms best_is =
    {
    false,
    "is",
    3,
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57340.000000	0.000000	17457.117647	0.000000	30.444921	0.000000	0.000000	0.000000	0.000000	0.000000	594.823529	0.000000	0.000000	0.000000	243.176471	0.000000	0.000000	0.000000	0.000000	0.000000	1.574792	0.000000	0.000000	0.000000	1.574792	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	70.981328	0.000000	0.000000	0.000000	29.018672	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
1	0.985600	17	57340.000000	0.000000	9578.352941	0.000000	16.704487	0.000000	0.000000	0.000000	0.000000	0.000000	595.470588	0.000000	0.000000	0.000000	242.529412	0.000000	0.000000	0.000000	0.000000	0.000000	1.699230	0.000000	0.000000	0.000000	1.699230	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.058543	0.000000	0.000000	0.000000	28.941457	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
2	0.985600	17	57340.000000	0.000000	3521.235294	0.000000	6.140975	0.000000	0.000000	0.000000	0.000000	0.000000	585.000000	0.000000	0.000000	0.000000	253.000000	0.000000	0.000000	0.000000	0.000000	0.000000	1.907909	0.000000	0.000000	0.000000	1.907909	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	69.809069	0.000000	0.000000	0.000000	30.190931	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
3	0.985600	17	57340.000000	0.000000	2128.058824	0.000000	3.711299	0.000000	0.000000	0.000000	0.000000	0.000000	571.529412	0.000000	0.000000	0.000000	266.470588	0.000000	0.000000	0.000000	0.000000	0.000000	1.795495	0.000000	0.000000	0.000000	1.795495	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	68.201600	0.000000	0.000000	0.000000	31.798400	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
4	0.985600	17	57340.000000	0.000000	1554.294118	0.000000	2.710663	0.000000	0.000000	0.000000	0.000000	0.000000	561.588235	0.000000	0.000000	0.000000	276.411765	0.000000	0.000000	0.000000	0.000000	0.000000	2.037858	0.000000	0.000000	0.000000	2.037858	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	67.015303	0.000000	0.000000	0.000000	32.984697	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
5	0.985600	17	57340.000000	0.000000	1213.058824	0.000000	2.115554	0.000000	0.000000	0.000000	0.000000	0.000000	549.411765	0.000000	0.000000	0.000000	288.588235	0.000000	0.000000	0.000000	0.000000	0.000000	2.148824	0.000000	0.000000	0.000000	2.148824	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	65.562263	0.000000	0.000000	0.000000	34.437737	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          1
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           9578.352941
rules%         16.704487
same%stdev     1.699230
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.699230
same%          71.058543
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          28.941457
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan

bests[9].suffixonly == [false]
bests[9].langbase == [is]
comp = comp_parms0_off
bests[9].rows == [3]
  R->R  W->R  R->W  W->W

    -4   746  -735   187
  -352   716  -976   600
   -64   192 -1664  -128

*/
//iteration:19.58 weights ret = ret + 5000.0/(5000.0 + rcount*rcount);
/* 14603 13595.621226 */
        {
           -4,  746, -735,  187, //535008
         -352,  716, -976,  600, //10333
          -64,  192,-1664, -128, //1303
            0,-1024,-2048,  512  //0
                                 //        0
        }
//OnlyZeros 12 
//suffix only no 
    };
#elif 0
static bestParms best_is =
    {
    false,
    "is",
    3,
//iteration:18.42 weights
/* 12033 11297.595310 */
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57340.000000	0.000000	17453.352941	0.000000	30.438355	0.000000	0.000000	0.000000	0.000000	0.000000	595.117647	0.000000	0.000000	0.000000	242.882353	0.000000	0.000000	0.000000	0.000000	0.000000	1.581923	0.000000	0.000000	0.000000	1.581923	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.016426	0.000000	0.000000	0.000000	28.983574	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
1	0.985600	17	57340.000000	0.000000	9576.058824	0.000000	16.700486	0.000000	0.000000	0.000000	0.000000	0.000000	595.941176	0.000000	0.000000	0.000000	242.058824	0.000000	0.000000	0.000000	0.000000	0.000000	1.677275	0.000000	0.000000	0.000000	1.677275	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.114699	0.000000	0.000000	0.000000	28.885301	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
2	0.985600	17	57340.000000	0.000000	3520.117647	0.000000	6.139026	0.000000	0.000000	0.000000	0.000000	0.000000	585.411765	0.000000	0.000000	0.000000	252.588235	0.000000	0.000000	0.000000	0.000000	0.000000	1.895769	0.000000	0.000000	0.000000	1.895769	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	69.858206	0.000000	0.000000	0.000000	30.141794	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
3	0.985600	17	57340.000000	0.000000	2126.411765	0.000000	3.708427	0.000000	0.000000	0.000000	0.000000	0.000000	571.764706	0.000000	0.000000	0.000000	266.235294	0.000000	0.000000	0.000000	0.000000	0.000000	1.793716	0.000000	0.000000	0.000000	1.793716	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	68.229679	0.000000	0.000000	0.000000	31.770321	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
4	0.985600	17	57340.000000	0.000000	1554.470588	0.000000	2.710971	0.000000	0.000000	0.000000	0.000000	0.000000	561.764706	0.000000	0.000000	0.000000	276.235294	0.000000	0.000000	0.000000	0.000000	0.000000	2.038936	0.000000	0.000000	0.000000	2.038936	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	67.036361	0.000000	0.000000	0.000000	32.963639	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
5	0.985600	17	57340.000000	0.000000	1211.705882	0.000000	2.113195	0.000000	0.000000	0.000000	0.000000	0.000000	549.352941	0.000000	0.000000	0.000000	288.647059	0.000000	0.000000	0.000000	0.000000	0.000000	2.140891	0.000000	0.000000	0.000000	2.140891	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	65.555244	0.000000	0.000000	0.000000	34.444756	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          1
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           9576.058824
rules%         16.700486
same%stdev     1.677275
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.677275
same%          71.114699
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          28.885301
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan

bests[9].suffixonly == [false]
bests[9].langbase == [is]
comp = comp_parms0_off
bests[9].rows == [3]
  R->R  W->R  R->W  W->W

    -4   746  -738   186
  -345   718  -976   601
   -57   197 -1661  -127
*/
        {
           -4,  746, -738,  186, //436098
         -345,  718, -976,  601, //8301
          -57,  197,-1661, -127, //1102
            7,-1019,-2045,  513  //0
                                 //        0
        }
//OnlyZeros 12 
//suffix only no 
    };
#elif 0

static bestParms best_is =
    {
    false,
    "is",
    3,
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57340.000000	0.000000	21632.823529	0.000000	37.727282	0.000000	0.000000	0.000000	0.000000	0.000000	566.588235	0.000000	0.000000	0.000000	271.411765	0.000000	0.000000	0.000000	0.000000	0.000000	1.599924	0.000000	0.000000	0.000000	1.599924	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	67.611961	0.000000	0.000000	0.000000	32.388039	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
1	0.985600	17	57340.000000	0.000000	13170.058824	0.000000	22.968362	0.000000	0.000000	0.000000	0.000000	0.000000	554.352941	0.000000	0.000000	0.000000	283.647059	0.000000	0.000000	0.000000	0.000000	0.000000	1.690210	0.000000	0.000000	0.000000	1.690210	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	66.151902	0.000000	0.000000	0.000000	33.848098	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
2	0.985600	17	57340.000000	0.000000	6774.647059	0.000000	11.814871	0.000000	0.000000	0.000000	0.000000	0.000000	541.058824	0.000000	0.000000	0.000000	296.941176	0.000000	0.000000	0.000000	0.000000	0.000000	1.712977	0.000000	0.000000	0.000000	1.712977	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	64.565492	0.000000	0.000000	0.000000	35.434508	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
3	0.985600	17	57340.000000	0.000000	4461.117647	0.000000	7.780114	0.000000	0.000000	0.000000	0.000000	0.000000	525.764706	0.000000	0.000000	0.000000	312.235294	0.000000	0.000000	0.000000	0.000000	0.000000	1.814436	0.000000	0.000000	0.000000	1.814436	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	62.740418	0.000000	0.000000	0.000000	37.259582	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
4	0.985600	17	57340.000000	0.000000	3216.000000	0.000000	5.608650	0.000000	0.000000	0.000000	0.000000	0.000000	508.352941	0.000000	0.000000	0.000000	329.647059	0.000000	0.000000	0.000000	0.000000	0.000000	1.669014	0.000000	0.000000	0.000000	1.669014	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	60.662642	0.000000	0.000000	0.000000	39.337358	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
5	0.985600	17	57340.000000	0.000000	2494.470588	0.000000	4.350315	0.000000	0.000000	0.000000	0.000000	0.000000	494.294118	0.000000	0.000000	0.000000	343.705882	0.000000	0.000000	0.000000	0.000000	0.000000	1.509891	0.000000	0.000000	0.000000	1.509891	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	58.984978	0.000000	0.000000	0.000000	41.015022	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          0
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           21632.823529
rules%         37.727282
same%stdev     1.599924
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.599924
same%          67.611961
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          32.388039
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan

bests[9].suffixonly == [false]
bests[9].langbase == [is]
comp = comp_parms0_off
bests[9].rows == [3]
  R->R  W->R  R->W  W->W

  -761   253    -3    -3
    30   902  -852   132
 -1536   768 -1536 -1024

*/
//iteration:19.54 weights
/* 18027 13013.328212 */
        {
         -761,  253,   -3,   -3, //569584
           30,  902, -852,  132, //16366
        -1536,  768,-1536,-1024, //616
         -256, 1792,-2048,-1024  //0
                                 //        0
        }
//OnlyZeros 12 
//suffix only no 
    };


#elif 0

static bestParms best_is =
    {
    false,
    "is",
    3,
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57340.000000	0.000000	18061.588235	0.000000	31.499107	0.000000	0.000000	0.000000	0.000000	0.000000	591.470588	0.000000	0.000000	0.000000	246.529412	0.000000	0.000000	0.000000	0.000000	0.000000	1.603291	0.000000	0.000000	0.000000	1.603291	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	70.581216	0.000000	0.000000	0.000000	29.418784	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
1	0.985600	17	57340.000000	0.000000	10025.352941	0.000000	17.484048	0.000000	0.000000	0.000000	0.000000	0.000000	592.235294	0.000000	0.000000	0.000000	245.764706	0.000000	0.000000	0.000000	0.000000	0.000000	1.725432	0.000000	0.000000	0.000000	1.725432	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	70.672469	0.000000	0.000000	0.000000	29.327531	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
2	0.985600	17	57340.000000	0.000000	3621.000000	0.000000	6.314963	0.000000	0.000000	0.000000	0.000000	0.000000	583.529412	0.000000	0.000000	0.000000	254.470588	0.000000	0.000000	0.000000	0.000000	0.000000	1.861211	0.000000	0.000000	0.000000	1.861211	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	69.633581	0.000000	0.000000	0.000000	30.366419	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
3	0.985600	17	57340.000000	0.000000	2173.176471	0.000000	3.789983	0.000000	0.000000	0.000000	0.000000	0.000000	570.411765	0.000000	0.000000	0.000000	267.588235	0.000000	0.000000	0.000000	0.000000	0.000000	1.756883	0.000000	0.000000	0.000000	1.756883	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	68.068230	0.000000	0.000000	0.000000	31.931770	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
4	0.985600	17	57340.000000	0.000000	1544.529412	0.000000	2.693633	0.000000	0.000000	0.000000	0.000000	0.000000	561.058824	0.000000	0.000000	0.000000	276.941176	0.000000	0.000000	0.000000	0.000000	0.000000	2.059757	0.000000	0.000000	0.000000	2.059757	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	66.952127	0.000000	0.000000	0.000000	33.047873	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
5	0.985600	17	57340.000000	0.000000	1212.470588	0.000000	2.114528	0.000000	0.000000	0.000000	0.000000	0.000000	548.411765	0.000000	0.000000	0.000000	289.588235	0.000000	0.000000	0.000000	0.000000	0.000000	2.151308	0.000000	0.000000	0.000000	2.151308	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	65.442931	0.000000	0.000000	0.000000	34.557069	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          1
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           10025.352941
rules%         17.484048
same%stdev     1.725432
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.725432
same%          70.672469
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          29.327531
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan

bests[9].suffixonly == [false]
bests[9].langbase == [is]
comp = comp_parms0_off
bests[9].rows == [3]
  R->R  W->R  R->W  W->W

    -3   187  -183    48
   -85   185  -241   149
    -9    53  -413   -31

Rule weighting: ret = ret + 50.0/(50.0 + rcount*rcount);
*/
//iteration:15.2 weights
/* 6773 6375.070938 */
        {
           -3,  187, -183,   48, //239208
          -85,  185, -241,  149, //4699
           -9,   53, -413,  -31, //590
            6, -252, -510,  129  //0
                                 //        0
        }
//OnlyZeros 12 
//suffix only no 
    };




#elif 0
static bestParms best_is =
    {
    false,
    "is",
    3,
        {
// Icelandic 74.631476 +/- 1.476088 at 0.9856 of dataset, 17 iterations, 20398.823529 = 35.575207% rules, cutoff = 0 SUFFIXONLY == 1 during production, SUFFIXONLY == 0 during parameter setting
/*
     0,  34, -30,   7,
    55, -11,  19,  65
*/
//iteration:283
/* 17487 */
    /*
0,       0,       0,       0,//       0
0,       0,       0,       0,//       0
0,       9,       -7,      2,//       636373
7,       4,       -12,     9,//       42035
*/
//OnlyZeros 12 
//iteration:16.35
/* 8192 */
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57340.000000	0.000000	17365.882353	0.000000	30.285808	0.000000	0.000000	0.000000	0.000000	0.000000	595.941176	0.000000	0.000000	0.000000	242.058824	0.000000	0.000000	0.000000	0.000000	0.000000	1.443415	0.000000	0.000000	0.000000	1.443415	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.114699	0.000000	0.000000	0.000000	28.885301	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
1	0.985600	17	57340.000000	0.000000	9567.470588	0.000000	16.685509	0.000000	0.000000	0.000000	0.000000	0.000000	596.411765	0.000000	0.000000	0.000000	241.588235	0.000000	0.000000	0.000000	0.000000	0.000000	1.594909	0.000000	0.000000	0.000000	1.594909	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	71.170855	0.000000	0.000000	0.000000	28.829145	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
2	0.985600	17	57340.000000	0.000000	3538.941176	0.000000	6.171854	0.000000	0.000000	0.000000	0.000000	0.000000	585.117647	0.000000	0.000000	0.000000	252.882353	0.000000	0.000000	0.000000	0.000000	0.000000	1.901780	0.000000	0.000000	0.000000	1.901780	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	69.823108	0.000000	0.000000	0.000000	30.176892	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
3	0.985600	17	57340.000000	0.000000	2142.647059	0.000000	3.736741	0.000000	0.000000	0.000000	0.000000	0.000000	572.294118	0.000000	0.000000	0.000000	265.705882	0.000000	0.000000	0.000000	0.000000	0.000000	1.828720	0.000000	0.000000	0.000000	1.828720	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	68.292854	0.000000	0.000000	0.000000	31.707146	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
4	0.985600	17	57340.000000	0.000000	1555.294118	0.000000	2.712407	0.000000	0.000000	0.000000	0.000000	0.000000	560.823529	0.000000	0.000000	0.000000	277.176471	0.000000	0.000000	0.000000	0.000000	0.000000	2.105806	0.000000	0.000000	0.000000	2.105806	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	66.924049	0.000000	0.000000	0.000000	33.075951	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan
5	0.985600	17	57340.000000	0.000000	1201.411765	0.000000	2.095242	0.000000	0.000000	0.000000	0.000000	0.000000	548.882353	0.000000	0.000000	0.000000	289.117647	0.000000	0.000000	0.000000	0.000000	0.000000	2.147923	0.000000	0.000000	0.000000	2.147923	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	65.499087	0.000000	0.000000	0.000000	34.500913	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.000000	0.000000	0.000000	100.000000	0.000000	nan

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          1
fraction        9856.000000
iterations      17
trainlines      57340.000000
rules           9567.470588
rules%         16.685509
same%stdev     1.594909
ambi1%stdev    0.000000
ambi2%stdev    0.000000
ambi3%stdev    0.000000
diff%stdev     1.594909
same%          71.170855
ambi1%         0.000000
ambi2%         0.000000
ambi3%         0.000000
diff%          28.829145
amb.rules%     0.000000
false_amb%     0.000000
false_not_amb% 0.000000
true_amb%      0.000000
true_not_amb%  100.000000
precision       0.000000
recall          nan

bests[9].suffixonly == [false]
bests[9].langbase == [is]
comp = comp_parms0_off
bests[9].rows == [3]
  R->R  W->R  R->W  W->W

     0     0     0     0
     0   267  -266    63
     7   261  -125   449
*/
    0  ,  0,   0,  0,//0
    0  ,267,-266, 63,//287175
    7  ,261,-125,449,//14213
    //569,955,-579, 63 //0
        }
    };
#else
#error best_is not declared
#endif

#if 0
// Based on list with 398930 full forms from CTS, University of Leeds (not publicly available) 
static bestParms best_ru =
    {
    false,
    "ru",
    3,
        {
    /*
     0, 10, -10,  1,
     0,  0,   0,  8,
    10,  0,  12,  0
    */
    // After introduction of parmsoff:
    //iteration:1021
    /* 33328 */
    /*
      0,256,-280, 26,//1185047
      1,277,-412, 24,//116109
    -32,176,-272,128,//1709
    */

//iteration:1215
/* 62652 */
       0,1023,-1117,107,//2139379
       9,1115,-1645, 95,//229605
    -121, 709,-1085,513,//3330
//263,1285,-1021,1537//0
        }
    };

static bestParms best_ru_suffix =
    {
    true,
    "ru",
    3,
        {
    // After introduction of parmsoff:
    //iteration:1085
    /* 45283 */
/*
     0,  53, -59,   2,//530032
     0, 146,-181,   0,//150055
     0, 363,-512, 116,//35222
    71,-379,-1021,449 //38711
*/
//iteration:1201
/* 69498 */
/*
cutoff          0
fraction        0.985600
iterations      6
trainlines      393185.000000 
rules           88034.333333 
rules%          22.390054
same%stdev      0.400184 
ambi1%stdev     0.084321 
ambi2%stdev     0.020349
ambi3%stdev     0.000000
diff%stdev      0.347374
same%           81.256165
ambi1%          0.272701
ambi2%          0.223383
ambi3%          0.000000 
diff%           18.247752
amb.rules%      0.696258 
false_amb%      0.646939 
false_not_amb%  1.784160 
true_amb%       0.049318 
true_not_amb%   97.519582
precision       0.036717
recall          0.026899
*/
    0,103, -116,  4,//787020
    0,290, -364,  1,//244334
    5,721,-1028,230,//91786
//135,-763,-2045,897//0
        }
/*
0	0.985600	2	3001378.000000	0.000000	345657.500000	0.000000	11.516627	0.000000	0.000000	0.000000	0.000000	0.000000	39459.000000	171.000000	132.500000	0.000000	4089.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.038700	0.022575	0.037087	0.000000	0.024187	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.982213	0.389948	0.302153	0.000000	9.325686	0.903037	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.775335	2.709112	0.127702	96.387850	0.076087	0.045016
1	0.985600	2	3001378.000000	0.000000	174555.500000	0.000000	5.815845	0.000000	0.000000	0.000000	0.000000	0.000000	39503.500000	219.500000	212.500000	0.000000	3916.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.098362	0.020962	0.037087	0.000000	0.082237	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.083691	0.500547	0.484585	0.000000	8.931178	1.164143	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.851729	2.524400	0.312414	96.311457	0.154977	0.110129
2	0.985600	2	3001378.000000	0.000000	53513.500000	0.000000	1.782964	0.000000	0.000000	0.000000	0.000000	0.000000	39780.000000	157.000000	164.000000	0.000000	3751.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.019350	0.003225	0.038700	0.000000	0.022575	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.714221	0.358022	0.373985	0.000000	8.553772	0.895056	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.701222	2.642981	0.193834	96.461963	0.121429	0.068328
3	0.985600	2	3001378.000000	0.000000	35667.500000	0.000000	1.188371	0.000000	0.000000	0.000000	0.000000	0.000000	39674.000000	141.500000	136.500000	0.000000	3900.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.022575	0.040312	0.030637	0.000000	0.012900	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.472498	0.322676	0.311274	0.000000	8.893551	0.851729	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.677278	2.662364	0.174450	96.485907	0.114094	0.061495
4	0.985600	2	3001378.000000	0.000000	27982.000000	0.000000	0.932305	0.000000	0.000000	0.000000	0.000000	0.000000	39543.000000	132.500000	132.000000	0.000000	4044.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.012900	0.014512	0.025800	0.000000	0.053212	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.173766	0.302153	0.301012	0.000000	9.223069	0.827784	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.663596	2.672626	0.164189	96.499590	0.110092	0.057878
5	0.985600	2	3001378.000000	0.000000	23424.500000	0.000000	0.780458	0.000000	0.000000	0.000000	0.000000	0.000000	39438.500000	123.000000	126.500000	0.000000	4164.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.053212	0.000000	0.004837	0.000000	0.058049	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.935465	0.280489	0.288470	0.000000	9.495576	0.811822	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.654474	2.679467	0.157347	96.508711	0.107309	0.055466

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          2
fraction        9856.000000
iterations      2
trainlines      3001378.000000
rules           53513.500000
rules%         1.782964
same%stdev     0.019350
ambi1%stdev    0.003225
ambi2%stdev    0.038700
ambi3%stdev    0.000000
diff%stdev     0.022575
same%          90.714221
ambi1%         0.358022
ambi2%         0.373985
ambi3%         0.000000
diff%          8.553772
amb.rules%     0.895056
false_amb%     0.701222
false_not_amb% 2.642981
true_amb%      0.193834
true_not_amb%  96.461963
precision       0.121429
recall          0.068328

bests[16].suffixonly == [true]
bests[16].langbase == [ru]
comp = comp_parms0_off
bests[16].rows == [3]
  R->R  W->R  R->W  W->W

     0   103  -116     4
     0   290  -364     1
     5   721 -1028   230
*/
    };
#else
// Based on much bigger data set (www.aot.ru) (https://github.com/bachan/libturglem/tree/master/share/turglem)
// Over 3 million full forms
static bestParms best_ru =
    {
    false,
    "ru",
    3,
#if 1
//iteration:16.27 count
/* 127813 119362.513238 */
        {
            0,  196, -203,    9, //8327849
            1,  147, -195,   27, //820691
           -1,  291, -743,  201  //3689
/*
0	0.985600	2	3001378.000000	0.000000	263646.000000	0.000000	8.784165	0.000000	0.000000	0.000000	0.000000	0.000000	39104.000000	544.500000	362.500000	0.000000	3841.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.025800	0.062887	0.043537	0.000000	0.045150	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.172672	1.241677	0.826644	0.000000	8.759008	2.410380	0.000000	0.000000	0.000000	0.000000	0.000000	nan	2.284958	2.711393	0.125422	94.878227	0.026712	0.044212
1	0.985600	2	3001378.000000	0.000000	141204.000000	0.000000	4.704639	0.000000	0.000000	0.000000	0.000000	0.000000	39449.000000	518.500000	411.500000	0.000000	3473.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.045150	0.082237	0.014512	0.000000	0.051600	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.959409	1.182386	0.938384	0.000000	7.919821	2.388717	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.955441	2.403539	0.433276	95.207744	0.099738	0.152733
2	0.985600	2	3001378.000000	0.000000	46031.500000	0.000000	1.533679	0.000000	0.000000	0.000000	0.000000	0.000000	39805.000000	384.500000	332.000000	0.000000	3330.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.022575	0.011287	0.016125	0.000000	0.027412	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.771231	0.876813	0.757092	0.000000	7.594865	1.831159	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.558652	2.564307	0.272508	95.604533	0.080390	0.096061
3	0.985600	2	3001378.000000	0.000000	28870.000000	0.000000	0.961892	0.000000	0.000000	0.000000	0.000000	0.000000	39782.000000	353.000000	289.500000	0.000000	3427.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.012900	0.029025	0.014512	0.000000	0.027412	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.718781	0.804980	0.660175	0.000000	7.816063	1.664690	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.432090	2.604214	0.232601	95.731096	0.075110	0.081994
4	0.985600	2	3001378.000000	0.000000	21713.000000	0.000000	0.723434	0.000000	0.000000	0.000000	0.000000	0.000000	39701.500000	330.000000	265.000000	0.000000	3555.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.033862	0.025800	0.041925	0.000000	0.017737	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.535209	0.752531	0.604305	0.000000	8.107954	1.547250	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.347715	2.637280	0.199535	95.815470	0.068925	0.070338
5	0.985600	2	3001378.000000	0.000000	17517.000000	0.000000	0.583632	0.000000	0.000000	0.000000	0.000000	0.000000	39631.000000	323.500000	239.000000	0.000000	3658.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.009675	0.004837	0.016125	0.000000	0.001612	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.374441	0.737709	0.545015	0.000000	8.342835	1.470857	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.300967	2.666925	0.169890	95.862218	0.061292	0.059887

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      2
trainlines      3001378.000000
rules           46031.500000
rules%         1.533679
same%stdev     0.022575
ambi1%stdev    0.011287
ambi2%stdev    0.016125
ambi3%stdev    0.000000
diff%stdev     0.027412
same%          90.771231
ambi1%         0.876813
ambi2%         0.757092
ambi3%         0.000000
diff%          7.594865
amb.rules%     1.831159
false_amb%     1.558652
false_not_amb% 2.564307
true_amb%      0.272508
true_not_amb%  95.604533
precision       0.080390
recall          0.096061

bests[15].suffixonly == [false]
bests[15].langbase == [ru]
comp = comp_parms0_off
bests[15].rows == [3]
  R->R  W->R  R->W  W->W

     0   196  -203     9
     1   147  -195    27
    -1   291  -743   201
*/
/* ==
            0,  392, -406,   18, //8327849
            2,  294, -390,   54, //820691
           -1,  291, -743,  201, //3689
         -153,  645, -509,  769  //0
                                 //        0
*/
        }
#else
//iteration:16.59 count
/* 126597 118267.248161 */
        {
            0,  394, -405,   19, //8225073
            3,  299, -387,   55, //818665
           -1,  293, -741,  201  //3650
         //-153,  645, -509,  769  //0
                                 //        0
        }
/*
0	0.985600	2	3001378.000000	0.000000	264098.500000	0.000000	8.799242	0.000000	0.000000	0.000000	0.000000	0.000000	39101.500000	534.500000	360.500000	0.000000	3855.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.017737	0.069337	0.046762	0.000000	0.040312	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.166971	1.218873	0.822083	0.000000	8.792073	2.377315	0.000000	0.000000	0.000000	0.000000	0.000000	nan	2.251893	2.711393	0.125422	94.911293	0.027094	0.044212
1	0.985600	2	3001378.000000	0.000000	141411.000000	0.000000	4.711536	0.000000	0.000000	0.000000	0.000000	0.000000	39444.000000	517.500000	412.500000	0.000000	3478.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.035475	0.082237	0.011287	0.000000	0.058049	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.948007	1.180106	0.940664	0.000000	7.931223	2.384156	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.946319	2.398978	0.437836	95.216866	0.101106	0.154341
2	0.985600	2	3001378.000000	0.000000	46099.000000	0.000000	1.535928	0.000000	0.000000	0.000000	0.000000	0.000000	39797.000000	379.000000	328.500000	0.000000	3347.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.006450	0.006450	0.020962	0.000000	0.008062	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.752987	0.864271	0.749111	0.000000	7.633631	1.807215	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.526726	2.556326	0.280489	95.636459	0.084131	0.098875
3	0.985600	2	3001378.000000	0.000000	28940.000000	0.000000	0.964224	0.000000	0.000000	0.000000	0.000000	0.000000	39773.000000	350.500000	284.500000	0.000000	3444.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.006450	0.053212	0.014512	0.000000	0.045150	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.698258	0.799279	0.648773	0.000000	7.853690	1.639606	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.402445	2.599653	0.237161	95.760741	0.077961	0.083601
4	0.985600	2	3001378.000000	0.000000	21759.500000	0.000000	0.724984	0.000000	0.000000	0.000000	0.000000	0.000000	39693.000000	329.500000	260.500000	0.000000	3569.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.022575	0.043537	0.030637	0.000000	0.035475	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.515826	0.751391	0.594044	0.000000	8.138739	1.527866	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.326051	2.635000	0.201815	95.837134	0.070715	0.071141
5	0.985600	2	3001378.000000	0.000000	17536.500000	0.000000	0.584282	0.000000	0.000000	0.000000	0.000000	0.000000	39625.500000	321.500000	235.000000	0.000000	3670.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.001612	0.011287	0.009675	0.000000	0.003225	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.361899	0.733148	0.535893	0.000000	8.369060	1.453754	0.000000	0.000000	0.000000	0.000000	0.000000	nan	1.283864	2.666925	0.169890	95.879321	0.062057	0.059887

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      2
trainlines      3001378.000000
rules           46099.000000
rules%         1.535928
same%stdev     0.006450
ambi1%stdev    0.006450
ambi2%stdev    0.020962
ambi3%stdev    0.000000
diff%stdev     0.008062
same%          90.752987
ambi1%         0.864271
ambi2%         0.749111
ambi3%         0.000000
diff%          7.633631
amb.rules%     1.807215
false_amb%     1.526726
false_not_amb% 2.556326
true_amb%      0.280489
true_not_amb%  95.636459
precision       0.084131
recall          0.098875

bests[15].suffixonly == [false]
bests[15].langbase == [ru]
comp = comp_parms0_off
bests[15].rows == [3]
  R->R  W->R  R->W  W->W

     0   394  -405    19
     3   299  -387    55
    -1   293  -741   201
*/
#endif
    };


static bestParms best_ru_suffix =
    {
    true,
    "ru",
    3,
/*
0	0.985600	2	3001378.000000	0.000000	389300.500000	0.000000	12.970725	0.000000	0.000000	0.000000	0.000000	0.000000	39464.500000	181.000000	135.000000	0.000000	4071.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.088687	0.003225	0.067724	0.000000	0.024187	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.994755	0.412752	0.307854	0.000000	9.284639	0.933823	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.808401	2.711393	0.125422	96.354784	0.071990	0.044212
1	0.985600	2	3001378.000000	0.000000	148398.500000	0.000000	4.944346	0.000000	0.000000	0.000000	0.000000	0.000000	39598.000000	212.000000	234.000000	0.000000	3808.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.148349	0.025800	0.054824	0.000000	0.067724	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.299188	0.483444	0.533613	0.000000	8.683754	1.188087	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.888215	2.536942	0.299872	96.274970	0.144426	0.105707
2	0.985600	2	3001378.000000	0.000000	50766.500000	0.000000	1.691440	0.000000	0.000000	0.000000	0.000000	0.000000	39761.000000	166.500000	179.500000	0.000000	3745.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.070949	0.030637	0.066112	0.000000	0.025800	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.670893	0.379686	0.409331	0.000000	8.540089	0.973730	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.773055	2.636140	0.200675	96.390130	0.114883	0.070740
3	0.985600	2	3001378.000000	0.000000	34596.500000	0.000000	1.152687	0.000000	0.000000	0.000000	0.000000	0.000000	39655.500000	150.500000	142.500000	0.000000	3903.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.027412	0.008062	0.072562	0.000000	0.053212	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.430311	0.343200	0.324957	0.000000	8.901532	0.901897	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.728587	2.663505	0.173310	96.434598	0.106294	0.061093
4	0.985600	2	3001378.000000	0.000000	27055.000000	0.000000	0.901419	0.000000	0.000000	0.000000	0.000000	0.000000	39534.500000	133.500000	141.500000	0.000000	4042.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.001612	0.027412	0.056437	0.000000	0.082237	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	90.154383	0.304433	0.322676	0.000000	9.218508	0.879093	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.716045	2.673766	0.163048	96.447140	0.102216	0.057476
5	0.985600	2	3001378.000000	0.000000	22910.500000	0.000000	0.763333	0.000000	0.000000	0.000000	0.000000	0.000000	39413.000000	124.000000	132.500000	0.000000	4182.500000	0.000000	0.000000	0.000000	0.000000	0.000000	0.029025	0.016125	0.066112	0.000000	0.079012	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	89.877315	0.282769	0.302153	0.000000	9.537763	0.868832	0.000000	0.000000	0.000000	0.000000	0.000000	nan	0.706923	2.674907	0.161908	96.456262	0.102750	0.057074

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          2
fraction        9856.000000
iterations      2
trainlines      3001378.000000
rules           50766.500000
rules%         1.691440
same%stdev     0.070949
ambi1%stdev    0.030637
ambi2%stdev    0.066112
ambi3%stdev    0.000000
diff%stdev     0.025800
same%          90.670893
ambi1%         0.379686
ambi2%         0.409331
ambi3%         0.000000
diff%          8.540089
amb.rules%     0.973730
false_amb%     0.773055
false_not_amb% 2.636140
true_amb%      0.200675
true_not_amb%  96.390130
precision       0.114883
recall          0.070740

bests[16].suffixonly == [true]
bests[16].langbase == [ru]
comp = comp_parms0_off
bests[16].rows == [3]
  R->R  W->R  R->W  W->W

     0    70   -78    11
     0   121  -162    20
    32  -192  -511   161
*/
//iteration:15.14 count
/* 108865 103648.985569 */
        {
//            0,    0,    0,    0, //0
            0,   70,  -78,   11, //1129078
            0,  121, -162,   20, //44737
           32, -192, -511,  161  //183867
                                 //        0
        }
//OnlyZeros 16 
//suffix only yes 
    };
#endif


static bestParms best_nl =
    {
    false,
    "nl",
    3,
        {
  /*
    0,   9,  -8,   2,
    7,  11,   9,   9
  */
// After introduction of parmsoff:
//iteration:1267
/* 12787 */
    // 90.294388 +/- 0.606321 at 0.9856 of dataset, 7 iterations, 18187.000000 = 6.112106% rules, cutoff = 1
     0,     738,     -926,   83,    //  463375
    -1,      16,     -255,    6,    //   65235
     8,     544,    -1144,   80     //     726
        }
    };


static bestParms best_nl_suffix =
    {
    true,
    "nl",
    3,
        {
// After introduction of parmsoff:
//iteration:578
/* 16759 */
    0,       5,      -5,    1,   //289215
    0,      31,     -31,    1,   // 21341
    0,      24,     -47,    1,   // 10435
        }
    };

static bestParms best_da_suffix =
    {
    true,
    "da",
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          0
fraction        9856.000000
iterations      5
trainlines      585173.000000
rules           61171.600000
rules%         10.453592
same%stdev     0.420111
ambi1%stdev    0.074450
ambi2%stdev    0.108916
ambi3%stdev    0.000000
diff%stdev     0.413541
same%          92.129861
ambi1%         0.480438
ambi2%         0.399572
ambi3%         0.000000
diff%          6.990130
amb.rules%     0.963254
false_amb%     0.582709
false_not_amb% 11.532881
true_amb%      0.380545
true_not_amb%  87.503865
precision       0.246154
recall          0.031943
/var/csttools/bin/makeaffixrules -i STOposUTF8.ph -c 0 -o flexrules.suffix.notags -e da -n 123 -f parms0 -s
*/
//OnlyZeros 12 
//suffix only yes 
//iteration:11.33
/* 11268 */
    3,
        {
            7,   31,  -24,    1, //94608
            1,   45,  -83,    4, //20473
           -3,   26,  -36,   14, //4091
          //-24,   48, -112,   32  //0
                                 //        0
        }

    };

static bestParms best_da =
    {
    false,
    "da",
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          1
fraction        9856.000000
iterations      5
trainlines      585173.000000
rules           33936.600000
rules%         5.799413
same%stdev     0.375131
ambi1%stdev    0.238159
ambi2%stdev    0.961303
ambi3%stdev    0.000000
diff%stdev     1.218753
same%          88.336306
ambi1%         1.198716
ambi2%         1.648234
ambi3%         0.000000
diff%          8.816744
amb.rules%     3.082412
false_amb%     2.321322
false_not_amb% 11.152337
true_amb%      0.761089
true_not_amb%  85.765252
precision       0.140845
recall          0.063885
*/
//OnlyZeros 12 
//suffix only no 
//iteration:11.61
/* 10481 */
     3,
        {
            0,   14,  -16,    2, //407715
           -2,   19,  -42,   12, //66486
            1,   45,  -50,   19, //317
//            8,    0,  -96,   88  //0
                                 //        0
        }

    };

/* transform singular to plural (only nouns) */
static bestParms best_daSP =
    {
    false,
    "daSP",
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57735.000000	6806.470588	4023.588235	11.789158	6.969063	758.411765	3.235294	2.176471	0.176471	80.000000	747.058824	7.235294	2.588235	0.000000	87.117647	0.819871	0.263374	0.146528	0.046558	0.896490	0.800827	0.344238	0.233456	0.000000	0.891467	89.859214	0.383329	0.257876	0.020909	9.478673	1.470588	88.514079	0.857262	0.306663	0.000000	10.321996	2.188458	0.689992	4.816002	0.780597	93.713410	0.361290	0.139477	1.282409	4.690549	0.906050	93.120992	0.261044	0.161893
1	0.985600	17	57735.000000	1690.647059	2244.117647	2.928288	3.886928	769.058824	0.000000	0.000000	0.000000	74.941176	749.470588	11.705882	10.588235	0.000000	72.235294	0.921081	0.000000	0.000000	0.000000	0.921081	0.821003	0.440873	0.393096	0.000000	0.910599	91.120714	0.000000	0.000000	0.000000	8.879286	0.000000	88.799833	1.386953	1.254530	0.000000	8.558684	3.589350	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	1.345135	3.352384	2.244215	93.058266	0.454802	0.400996
2	0.985600	17	57735.000000	1236.294118	1105.705882	2.141325	1.915140	765.117647	0.000000	0.000000	0.000000	78.882353	750.058824	15.294118	16.588235	0.000000	62.058824	0.985874	0.000000	0.000000	0.000000	0.985874	0.724919	0.466027	0.445415	0.000000	0.823514	90.653750	0.000000	0.000000	0.000000	9.346250	0.000000	88.869529	1.812099	1.965431	0.000000	7.352941	4.606914	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	1.080290	2.069975	3.526624	93.323111	0.620098	0.630137
3	0.985600	17	57735.000000	1031.117647	834.294118	1.785949	1.445040	759.588235	0.000000	0.000000	0.000000	84.411765	747.529412	16.529412	15.882353	0.000000	64.058824	0.946996	0.000000	0.000000	0.000000	0.946996	0.796044	0.591109	0.399349	0.000000	0.892046	89.998606	0.000000	0.000000	0.000000	10.001394	0.000000	88.569836	1.958461	1.881795	0.000000	7.589908	4.523279	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	1.136047	2.209367	3.387232	93.267354	0.598522	0.605230
4	0.985600	17	57735.000000	903.117647	693.058824	1.564246	1.200414	753.882353	0.000000	0.000000	0.000000	90.117647	744.647059	18.000000	14.352941	0.000000	67.000000	0.985874	0.000000	0.000000	0.000000	0.985874	0.937575	0.569768	0.507790	0.000000	1.071276	89.322554	0.000000	0.000000	0.000000	10.677446	0.000000	88.228325	2.132701	1.700585	0.000000	7.938389	4.495400	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	1.247561	2.348759	3.247839	93.155841	0.565534	0.580324
5	0.985600	17	57735.000000	801.647059	614.705882	1.388494	1.064702	748.764706	0.000000	0.000000	0.000000	95.235294	742.823529	17.529412	13.764706	0.000000	69.882353	1.046023	0.000000	0.000000	0.000000	1.046023	1.038197	0.571487	0.484062	0.000000	1.083397	88.716197	0.000000	0.000000	0.000000	11.283803	0.000000	88.012267	2.076945	1.630889	0.000000	8.279900	4.314190	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	1.240591	2.523000	3.073599	93.162810	0.553325	0.549191

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      17
trainlines      57735.000000
rules           1105.705882
rules%         1.915140
same%stdev     0.724919
ambi1%stdev    0.466027
ambi2%stdev    0.445415
ambi3%stdev    0.000000
diff%stdev     0.823514
same%          88.869529
ambi1%         1.812099
ambi2%         1.965431
ambi3%         0.000000
diff%          7.352941
amb.rules%     4.606914
false_amb%     1.080290
false_not_amb% 2.069975
true_amb%      3.526624
true_not_amb%  93.323111
precision       0.620098
recall          0.630137
*/
//OnlyZeros 12
//suffix only no
//iteration:19.45
/* 3237 */
     2,
        {
            1, 1068,-1263,  309, //172958
        //-32,   64, -544, 1024
           -1,    2,  -17,   32  //48
                                 //        0
        }
    };



static bestParms best_daSP_suffix =
    {
    true,
    "daSP",
/*
cutoff	fraction	iterations	trainlines	suffixrules	affixrules	suffix%	affix%	s-same	s-ambiguous			s-different	a-same	a-ambiguous			a-different	s-same-stddev%	s-ambiguous-stddev%			s-different-stddev%	a-same-stddev%	a-ambiguous-stddev%			a-different-stddev%	s-same%	s-ambiguous%			s-different%	s-amb.rules%	a-same%	a-ambiguous%			a-different%	a-amb.rules%	s_false_amb	s_false_not_amb	s_true_amb	s_true_not_amb	s_precision	s_recall	a_false_amb	a_false_not_amb	a_true_amb	a_true_not_amb	a_precision	a_recall
0	0.985600	17	57735.000000	6806.470588	4350.176471	11.789158	7.534730	758.411765	3.235294	2.176471	0.176471	80.000000	757.235294	4.411765	1.529412	0.000000	80.823529	0.819871	0.263374	0.146528	0.046558	0.896490	0.820374	0.274880	0.192502	0.000000	0.867881	89.859214	0.383329	0.257876	0.020909	9.478673	1.470588	89.719822	0.522721	0.181210	0.000000	9.576248	1.582102	0.689992	4.816002	0.780597	93.713410	0.361290	0.139477	0.731809	4.746306	0.850293	93.671592	0.367470	0.151930
1	0.985600	17	57735.000000	1690.647059	2643.470588	2.928288	4.578628	769.058824	0.000000	0.000000	0.000000	74.941176	754.235294	7.058824	8.941176	0.000000	73.764706	0.921081	0.000000	0.000000	0.000000	0.921081	0.900912	0.285563	0.328433	0.000000	0.835213	91.120714	0.000000	0.000000	0.000000	8.879286	0.000000	89.364371	0.836353	1.059381	0.000000	8.739894	2.752997	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	0.836353	3.679955	1.916643	93.567048	0.533981	0.342466
2	0.985600	17	57735.000000	1236.294118	1213.941176	2.141325	2.102609	765.117647	0.000000	0.000000	0.000000	78.882353	754.941176	8.647059	20.352941	0.000000	60.058824	0.985874	0.000000	0.000000	0.000000	0.985874	0.982833	0.433197	0.592330	0.000000	0.842473	90.653750	0.000000	0.000000	0.000000	9.346250	0.000000	89.448007	1.024533	2.411486	0.000000	7.115974	4.369947	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	0.738779	1.965431	3.631168	93.664622	0.710778	0.648817
3	0.985600	17	57735.000000	1031.117647	929.529412	1.785949	1.609993	759.588235	0.000000	0.000000	0.000000	84.411765	752.470588	9.647059	18.294118	0.000000	63.588235	0.946996	0.000000	0.000000	0.000000	0.946996	0.889728	0.429127	0.526160	0.000000	0.864664	89.998606	0.000000	0.000000	0.000000	10.001394	0.000000	89.155283	1.143016	2.167549	0.000000	7.534151	4.112071	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	0.808475	2.293003	3.303596	93.594926	0.671388	0.590286
4	0.985600	17	57735.000000	903.117647	772.000000	1.564246	1.337144	753.882353	0.000000	0.000000	0.000000	90.117647	749.764706	9.647059	17.411765	0.000000	67.176471	0.985874	0.000000	0.000000	0.000000	0.985874	0.944649	0.495548	0.479563	0.000000	1.014258	89.322554	0.000000	0.000000	0.000000	10.677446	0.000000	88.834681	1.143016	2.063005	0.000000	7.959297	4.000558	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	0.899080	2.495121	3.101478	93.504321	0.633001	0.554172
5	0.985600	17	57735.000000	801.647059	675.176471	1.388494	1.169440	748.764706	0.000000	0.000000	0.000000	95.235294	747.352941	9.705882	16.235294	0.000000	70.705882	1.046023	0.000000	0.000000	0.000000	1.046023	1.031164	0.454590	0.473062	0.000000	1.128523	88.716197	0.000000	0.000000	0.000000	11.283803	0.000000	88.548927	1.149986	1.923613	0.000000	8.377474	3.798439	0.000000	5.596599	0.000000	94.403401	0.000000	0.000000	0.899080	2.697240	2.899359	93.504321	0.617211	0.518057

New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          2
fraction        9856.000000
iterations      17
trainlines      57735.000000
rules           1213.941176
rules%         2.102609
same%stdev     0.982833
ambi1%stdev    0.433197
ambi2%stdev    0.592330
ambi3%stdev    0.000000
diff%stdev     0.842473
same%          89.448007
ambi1%         1.024533
ambi2%         2.411486
ambi3%         0.000000
diff%          7.115974
amb.rules%     4.369947
false_amb%     0.738779
false_not_amb% 1.965431
true_amb%      3.631168
true_not_amb%  93.664622
precision       0.710778
recall          0.648817

*/
//OnlyZeros 12 
//suffix only yes 
//iteration:19.63
/* 3439 */
    4,
        {
            1,  587, -939,   21, //37164
          -64,  527, -831,  193, //4158
       //-304,  768,-1088,  464, //967
          -19,   48,  -68,   29 //967
       //-960,-1008,-2048,  496  //0
          //-60,  -63, -128,   31  //0
                                 //        0
        }
    };


/* Based on 10000 words
int best_en[]       =
    {
    0, 9,-5, 4,
    0,-4,-5, 0,
    1, 0, 0,-1
    };
*/
static bestParms best_en =
    {
    false,
    "en",
    3,
        {
/* 9986 */
    /* Slightly worse than the values below!
    0,  8, -9,  2,
    1,-10, 11,  0,
    0,  4,-13,  0
    */
/* 10251 */
// 88.763636 +- 1.328464, cutoff 2, 1477.666667 rules, 15 iterations
/*
    0, 4,-5, 1,
    0,-5, 6, 0,
    0, 1, 0, 0
*/
/* 10397 20091021, using 100% of available training data, after 211 iterations */
// 88.484848 +- 1.471358, cutoff 2, 1672.200000 rules, 15 iterations
/*    
    0, 4,-3, 1,
    5, 4,-9,-4,*/
/* 10381 */
// iteration:395
/*
    0,31,-25,8,     //443260
    16,14,-36,-10,  //52623
*/
/* 9859 */
/*  0,35,-36, 2,
    1, 0,  5, 2,
    0,31,-40, 0
*/
// After introduction of parmsoff:
//iteration:94
/* 10398 */
/* 10398 */
/*
    0, 4,-3, 1,// 437782
    2, 1,-4,-2,// 58708
    7, 5, 3, 1,// 4
*/
//iteration:811
/* 10374 */
// 88.466667 +- 1.510243, cutoff 2, 1671.600000 = 2.221838% rules, 15 iterations
    0,  0,   0, 0,  //0
    0,123,-101,30,  //355497
    7, 77, -63,47   //38334
        }
    };

static bestParms best_ena = // English, ambiguous training pairs in training set derived from CELEX
    {
    false,
    "ena",
    3,
//OnlyZeros 12 
//suffix only no 
//iteration:19.63
/* 8960 8539.535721 */
        {
            1,  511, -348,  180, //325539
           -7,  696, -328,  573, //5597
           64,  512,-1536, -576  //592
                                 //        0
        }
    };




static bestParms best_en_suffix =
    {
    true,
    "en",
    3,
        {
    /* 14144, after 88 iterations*/
// 87.993939 +- 1.039632, cutoff 3, 1156.666667 rules, 15 iterations 
    /*
     1, 3, 0, 1,//421555
    -1, 2,-2, 1,//11288
     1, 6, 3, 0,//19
     */
// After introduction of parmsoff:
//iteration:258
/* 11917 */
/*
    0,  0,  0,  0,// 0
    0,  0,  0,  0,// 0
    0,  4,  0,  4,// 349160
    0,  7,-10, -2 // 22726
*/
//iteration:350
/* 11888 */
// 88.315152 +- 1.084324, cutoff 2, 1650.466667 = 2.193748% rules, 15 iterations    
//     0, 1,  0, 1, // 105490
     0, 1,  0, 1, //  13535
     0, 9, -7, 3, //  24911
    -7,-8,-12, 4  //   9828
        }
    };

static bestParms best_ge =
    {
    false,
    "ge",
    3,
        {
//based on 77566 pairs (deutsch_affix_2464_3.UTF8.txt)
//iteration:83
/* 11434 */
/*
    0, 5,-5, 0,//425724
    0, 3,-1, 2,//142447
    7, 5, 2, 1//51604
*/
// After introduction of parmsoff:
//iteration:322
/* 10970 */       //         
/*
     0, 0, 0, 0,  //0        Apparently, no parameter settings 
                  //         dedicated to this group are needed.
     0, 6,-5, 1,  //495542
     0, 8,-4, 4,  //1517
    11,19, 1,11   //53248
*/
//iteration:939
/* 10867 */
0,       0,       0,       0,  //          0
0,       189,     -165,    28, //     401859
7,       165,     -213,    121//      31601
//263,     421,     -93 ,    353,//          0
        }
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          0
fraction        9856.000000
iterations      7
trainlines      313549.000000
rules           31345.285714
rules%         9.996934
same%stdev     0.349160
ambi1%stdev    0.053246
ambi2%stdev    0.077088
ambi3%stdev    0.000000
diff%stdev     0.333239
same%          89.773648
ambi1%         0.296190
ambi2%         0.346075
ambi3%         0.000000
diff%          9.584087
amb.rules%     0.813743
false_amb%     0.751387
false_not_amb% 1.730374
true_amb%      0.062356
true_not_amb%  97.455883
precision       0.039841
recall          0.034783

bests[7].suffixonly == [false]
bests[7].langbase == [ge]
comp = comp_parms0_off
bests[7].rows == [3]
  R->R  W->R  R->W  W->W

     0     0     0     0
     0   189  -165    28
     7   165  -213   121
*/
    };


//based on 77566 pairs (deutsch_affix_2464_3.UTF8.txt)
//iteration:83
/* 11434 */
 static bestParms best_ge_suffix =
     {
     true,
     "ge",
     4,
//based on 77566 pairs (deutsch_affix_2464_3.UTF8.txt)
//iteration:77
/* 14206 */
        {
/*
    0, 3,-2, 2,//	408111
    0, 3, 0,-2,//	15424
    7, 5, 3, 1//	15449
*/
// After introduction of parmsoff:
//iteration:207
/* 13549 */
/*
     0, 3,  0, 3,//              290921
     0, 1, -1, 0,//              36836
     0, 5, -3, 2,//              15198
    -8,-6,-11, 1 //              12888
*/
//iteration:539
/* 13346 */
0,       0,       0,       0, //      0
//0,       0,       0,       0, //      0
0,       21,      -19,     2, //      172352
-11,     25,      -25,     25//      11892
         }
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          0
fraction        9856.000000
iterations      7
trainlines      313549.000000
rules           40932.857143
rules%         13.054692
same%stdev     0.351426
ambi1%stdev    0.068355
ambi2%stdev    0.048801
ambi3%stdev    0.000000
diff%stdev     0.367948
same%          89.751824
ambi1%         0.149654
ambi2%         0.130947
ambi3%         0.000000
diff%          9.967575
amb.rules%     0.386606
false_amb%     0.324250
false_not_amb% 1.730374
true_amb%      0.062356
true_not_amb%  97.883021
precision       0.087719
recall          0.034783

bests[8].suffixonly == [true]
bests[8].langbase == [ge]
comp = comp_parms0_off
bests[8].rows == [4]
  R->R  W->R  R->W  W->W

     0     0     0     0
     0    21   -19     2
   -11    25   -25    25
     0     0     04273647
*/
    };

static bestParms best_la =
    {
    false,
    "la",
    3,
//OnlyZeros 12 
//suffix only no 
//iteration:19.50
/* 6429 */
        {
            1,  289, -273,   27, //210484
           16,  808, -408,  512, //6084
         -448,  384, -896,  320//, //615
//          256,  256,-1280, 1024  //0
                                 //        0
        }
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      31
trainlines      17954.000000
rules           978.612903
rules%         5.450668
same%stdev     2.828105
ambi1%stdev    0.547465
ambi2%stdev    0.472805
ambi3%stdev    0.000000
diff%stdev     2.688214
same%          50.361830  (old lemmatiser: 49.540047, cutoff 1)
ambi1%         0.699129
ambi2%         0.870845
ambi3%         0.000000
diff%          48.068196
amb.rules%     2.845578
false_amb%     2.551208
false_not_amb% 6.328959
true_amb%      0.294370
true_not_amb%  90.825463
precision       0.054545
recall          0.044444
*/
    };

static bestParms best_la_suffix =
    {
    true,
    "la",
    3,
//OnlyZeros 12 
//suffix only yes 
//iteration:19.61
/* 7929 */
        {
          368, 1051,  -20, -440, //305094
        -1416, -310,-1984,    8, //6613
         -832,  416,-1488,   16//, //2386
//          512,  512,-1280,  256  //0
                                 //        0
        }
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          0
fraction        9856.000000
iterations      31
trainlines      17954.000000
rules           9680.612903
rules%         53.918976
same%stdev     3.264852
ambi1%stdev    0.585340
ambi2%stdev    0.825352
ambi3%stdev    0.000000
diff%stdev     3.120582
same%          46.363302  (old lemmatiser: 49.540047, cutoff 1)
ambi1%         1.275604
ambi2%         1.373727
ambi3%         0.000000
diff%          50.987367
amb.rules%     6.181774
false_amb%     5.458114
false_not_amb% 5.899669
true_amb%      0.723660
true_not_amb%  87.918558
precision       0.062171
recall          0.109259
*/
    };

static bestParms best_el =
    {
    false,
    "el",
    3,
//iteration:19.59
/* 31634 */
        {
            0,  611, -655,  126, //1141579
            2,  356, -317,   75, //162502
           32,  576,-1024,   32//, //2011
//          256,  768, -896, 1024  //0
                                 //        0
        }
/*
Old 86.901131	0.000000	0.000000	0.000000	13.098869 cutoff 1
or  87.638957	0.000000	0.000000	0.000000	12.361043 cutoff 1
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     no
cutoff          2
fraction        9856.000000
iterations      5
trainlines      556568.000000
rules           14376.800000
rules%         2.583117
same%stdev     0.351921
ambi1%stdev    0.139939
ambi2%stdev    0.166670
ambi3%stdev    0.000000
diff%stdev     0.286208
same%          84.259715
ambi1%         1.328087
ambi2%         2.788982
ambi3%         0.000000
diff%          11.623217
amb.rules%     4.604033
false_amb%     1.389572
false_not_amb% 2.506149
true_amb%      3.214461
true_not_amb%  92.889818
precision       0.536315
recall          0.561909
*/
    };

static bestParms best_el_suffix =
    {
    true,
    "el",
    3,
#if 0
//iteration:19.60
/* 35707 */
        {
            0,   30,  -32,    1, //872036
            1,   15,   -5,   11, //90168
           -9,  848,-1057,  404, //12291
        -2944,-1536,-2560, -128  //0
                                 //        0
        }
//OnlyZeros 12
//suffix only yes
#endif
//iteration:19.33
/* 62522 */
        {
            0,    0,    0,    0, //0
            0,  988, -929,  243, //1334941
            1,  549, -531,   72//, //81200
//        -1024,  128, -672, 1312  //0
                                 //        0
        }
/*
New algorithm, least wrongly lemmatised (MIN(diff)).
Suffix only     yes
cutoff          2
fraction        9856.000000
iterations      5
trainlines      556568.000000
rules           13966.000000
rules%         2.509307
same%stdev     0.287473
ambi1%stdev    0.091859
ambi2%stdev    0.174558
ambi3%stdev    0.000000
diff%stdev     0.172729
same%          85.993606
ambi1%         0.499262
ambi2%         2.619282
ambi3%         0.000000
diff%          10.887850
amb.rules%     3.357108
false_amb%     0.558288
false_not_amb% 2.678308
true_amb%      2.798819
true_not_amb%  93.964584
precision       0.714824
recall          0.511001
*/
    };


static struct bestParms bests[] =    
    {best_da
    ,best_da_suffix
    ,best_daSP
    ,best_en
    ,best_en_suffix
    ,best_ena

    ,best_fr
    ,best_fr_suffix
    ,best_ge
    ,best_ge_suffix
    ,best_is
    ,best_is_suffix
    ,best_la
    ,best_la_suffix
    ,best_el
    ,best_el_suffix
    ,best_ru
    ,best_ru_suffix
    ,best_nl
    ,best_nl_suffix
	,best_pl
    ,best_sl
    ,best_sl_suffix
    ,best_sv
    ,best_sv_suffix
    };

#if 0
static double best[NPARMS]    = 
  {// Before a chain of sibling rules is created, there are N candidate rules.
   // For each sibling created, the number n of remaining candidates is
   // decremented.
  0,  0,  0,  0, // n >= exp(N,0.75)
  0,  0,  0,  0, // n >= exp(N,0.5)
//0,  0,  0,  0, // n >= exp(N,0.25)
  0,  0,  0,  0  // n <  exp(N,0.25)
  };
#else
static double best[NPARMS]    =
   ///* R_R  W_R  R_W  W_W */   
   /* W_R  R_W  W_W  R_R */   
    { 1.0, 0.0, 0.0, 0.0
    , 0.0, 1.0, 0.0, 0.0
    , 0.0, 0.0, 1.0, 0.0
    , 0.0, 0.0, 0.0, 1.0
    };
#endif
//iteration:19.63
/* 387 381.726058 */
/*
0.000000,	10.000000,	-10.000000,	2.000000,
-13.581532,	-5.623026,	-4.283725,	6.696503,
-2.889162,	1.349916,	0.464059,	-4.429286,
-0.305579,	0.715326,	0.815889,	0.502817
*/
//iteration:19.63
/* 391 387.069724 */
/*
0.000000,	6.000000,	-6.000000,	2.000000,
-6.337474,	0.399136,	-1.608900,	-6.024108,
-3.947893,	0.883151,	2.100111,	3.650880,
0.463499,	1.732265,	1.476537,	-0.767184
*/

// Orthogonalised:

//iteration:19.39
/* number of nodes: 386, nodes/line: 0.387550 weight: 381.768147 blobs 1 lines 1254 * fraction 0.794328 = 996 lines*/
/*
        {         	         	         	          // # decisions
        0.000000,	4.100000,	-6.000000,	1.000000, //9886
        -4.167700,	0.043205,	-0.140436,	-1.019756, //784
        0.040775,	0.183468,	0.096672,	-0.172187, //17
        0.001182,	-0.003050,	-0.002846,	-0.004568  //0
        }         	         	         	          //(0 unresolved comparisons)
*/



//int D[NPARMS]        = {0};

static double  R[] = // 64 elements
    {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
    1,-1, 0, 0,
    1, 1, 0, 0,
    1, 0,-1, 0,
    1, 0, 1, 0,
    1, 0, 0, 1,
    1, 0, 0,-1,
    0, 1,-1, 0,
    0, 1, 1, 0,
    0, 1, 0,-1,
    0, 1, 0, 1,
    0, 0, 1,-1,
    0, 0, 1, 1
    };

#define NUMBER_OF_ROWS_IN_R (sizeof(R)/(ROWPARMS*sizeof(R[0])))
#if NPARMS % ROWPARMS != 0
#error NPARMS must be a multiple of ROWPARMS
#endif

static void plus(double * dest, double * term,int cols)
    {
    for(int col = 0;col < cols;++col)
        dest[col] += term[col];
    }

static void copy(double * dest,double * source,int cols)
    {
    for(int col = 0;col < cols;++col)
        dest[col] = source[col];
    }

static int pcnt[(NPARMS >> 2)+1] = {0,0,0,0,0};
static int improvements = 0;

void betterfound(int Nnodes,double weight,int swath,int iterations,const char * besttxt,int blobs,int lines,double fraction,int fraclines)
    {
    ++improvements;
    copy(best,parms.Matrix,NPARMS);
    FILE * f = fopen(besttxt,"a");
    if(f)
        {
        fprintf(f,"//iteration:%d.%d\n",swath,iterations);
        fprintf(f,"/* number of nodes: %d, nodes/line: %f weight: %f blobs %d lines %d * fraction %f = %d lines*/\n",
            Nnodes,(double)Nnodes/(double)fraclines,weight,blobs,lines,fraction,fraclines);
        fprintf(f,"        {         \t         \t         \t          // # decisions\n        ");
        int i = 0;
        for(;i < NPARMS;++i)
            {
            fprintf(f,"%f",parms.Matrix[i]);
            if(((i+1) % ROWPARMS) == 0)
                {
                if(i == NPARMS - 1)
                    fprintf(f,"  //%d\n        ",pcnt[i >> 2]);
                else
                    fprintf(f,", //%d\n        ",pcnt[i >> 2]);
                }
            else
                fprintf(f,",\t");
            }
        fprintf(f,"}         \t         \t         \t          //(%d unresolved comparisons)\n\n",pcnt[NPARMS >> 2]);
        fclose(f);
        }
    }

void worsefound()
    {
    copy(parms.Matrix,best,NPARMS);
    }
static int minparmsoff = 0;

void copybest()
    {
    copy(parms.Matrix,best,NPARMS); // go on with best result so far.
    }

//static const char * besttxt;
//static const char * parmstxt = NULL;
//static int OnlyZeros = NPARMS;
/*
static bool allZeros()
    {
    double x2 = 0;
    int i;
    for(i = 0;i < NPARMS;++i)
        {
        x2 += parms.Matrix[i]*parms.Matrix[i];
        }
    return x2 <= 0;
    }
*/


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

bool brown(/*const char * parmstxt*/)
    {
//    testAngle();
    static int it = 0;
    static double delta = 0.99;
    static double inc = (1.0 - delta) / 1365.0;
    int index = it;// - improvements; // to ensure that successful transformation is repeated.
    int r = (index / 2) % 6;
    double pone;
    if(index & 1)
        pone = 1.0;
    else
        pone = -1.0;
    setMinMaxIntegral(4);
    double ang = angle(delta,4);
//    printf("%d delta %1.4f angle %3.2f\n",it,delta,(180.0/pi)*ang);
    struct rotation Rotation;
    Rotation = makeRotationMatrix(r, pone, sin(ang));
    if((rand() % 3) == 2)
        {
        int G;
        do {G = rand() % 6;} while(G == r);
        double bigang = angle(0.5,4);
        struct rotation Forth = makeRotationMatrix(G, 1.0, sin(bigang));
        struct rotation Back = makeRotationMatrix(G, -1.0, sin(bigang));
        Rotation = rotate(Rotation,Forth);
        Rotation = rotate(Back,Rotation);
        }

    delta += inc;
    if(delta > 1.0 - inc)
        delta = 1.0 - inc;
    copy(parms.Matrix,best,NPARMS); // go on with best result so far.
    orthogonalise("brownStart");
    parms = rotate(parms,Rotation);
    ++it;
    return false;
    }

bool init()
    {
//    if(allZeros())
        {
        copy(parms.Matrix,best,NPARMS);
        }
    return true;
    }



void onlyZeros(const char * parmstxt,bool suffixonly)
    {
//    OnlyZeros = 0;
    for(unsigned int i = 0;i < sizeof(pcnt)/sizeof(pcnt[0]);++i)
        {
        if(pcnt[i] != 0)
            {
  //          OnlyZeros = (i+1) << 2;
            pcnt[i] = 0;
            }
        }
    if(parmstxt)
        {
        FILE * f = fopen(parmstxt,"a");
        assert(f);
    //    fprintf(f,"//OnlyZeros %d \n",OnlyZeros);
        fprintf(f,"//suffix only %s \n",suffixonly ? "yes" : "no");
        fclose(f);
        }
    }


void printparms(int Nnodes,double weight,const char * parmstxt)
    {
    int i;
    FILE * f = fopen(parmstxt,"a");
    assert(f);
    fprintf(f,"/* %d %f */\n",Nnodes,weight);
    fprintf(f,"        {\n        ");
    for(i = 0;i < NPARMS;++i)
        {
        fprintf(f,"%f",parms.Matrix[i]);
        if(((i+1) % ROWPARMS) == 0)
            {
            if(i == NPARMS - 1)
                fprintf(f,"  //%d\n        ",pcnt[i >> 2]);
            else
                fprintf(f,", //%d\n        ",pcnt[i >> 2]);
            }
        else
            fprintf(f,",");
        }
    fprintf(f,"                         //%9d\n",pcnt[i >> 2]);
    fprintf(f,"        }\n");
    fclose(f);
    }

#if 0
static int comp_parms(const vertex * a,const vertex * b)
    {
    //for(int o = 0;o < NPARMS;o += ROWPARMS)
    if(  a->R__R != b->R__R
      || a->W__R != b->W__R
      || a->R__W != b->R__W
      || a->W__W != b->W__W
      )
        {
        int off = minparmsoff;
        if(off < parmsoff)
            off = parmsoff;
        double A = 0.0;
        double B = 0.0;
        double N = a->R__R + a->W__R + a->R__W + a->W__W;
        //double N = sqrt((double)(a->R__R*a->R__R + a->W__R*a->W__R + a->R__W*a->R__W + a->W__W*a->W__W));
        double D = 20.0/N;
        double e = 4.0;
        for(int o = 0;o < NPARMS;o += ROWPARMS)
            {
            double x = parms.Matrix[o]*a->R__R + parms.Matrix[o+1]*a->W__R + parms.Matrix[o+2]*a->R__W + parms.Matrix[o+3]*a->W__W;
            double y = parms.Matrix[o]*b->R__R + parms.Matrix[o+1]*b->W__R + parms.Matrix[o+2]*b->R__W + parms.Matrix[o+3]*b->W__W;

            if(x < 0.0)
                A -= pow(D*-x,e);
            else
                A += pow(D*x,e);

            if(y < 0.0)
                B -= pow(D*-y,e);
            else
                B += pow(D*y,e);

            e -= 1.0;
            }
        ++pcnt[NPARMS >> 2];
        return A > B ? -1 : 1;
        }
    return 0;
    }
#else
static int comp_parms(const vertex * a,const vertex * b)
    {
    //for(int o = 0;o < NPARMS;o += ROWPARMS)
    if(  a->R__R != b->R__R
      || a->W__R != b->W__R
      || a->R__W != b->R__W
      || a->W__W != b->W__W
      )
        {
        int off = minparmsoff;
        if(off < parmsoff)
            off = parmsoff;
        for(int o = off;o < NPARMS;o += ROWPARMS)
            {
            //double A = parms.Matrix[o]*a->R__R + parms.Matrix[o+1]*a->W__R + parms.Matrix[o+2]*a->R__W + parms.Matrix[o+3]*a->W__W;
            //double B = parms.Matrix[o]*b->R__R + parms.Matrix[o+1]*b->W__R + parms.Matrix[o+2]*b->R__W + parms.Matrix[o+3]*b->W__W;
            double A = parms.Matrix[o+3]*a->R__R + parms.Matrix[o+0]*a->W__R + parms.Matrix[o+1]*a->R__W + parms.Matrix[o+2]*a->W__W;
            double B = parms.Matrix[o+3]*b->R__R + parms.Matrix[o+0]*b->W__R + parms.Matrix[o+1]*b->R__W + parms.Matrix[o+2]*b->W__W;
            if(A != B)
                {
                ++pcnt[o >> 2]; // For counting the number of times the first, second, third or fourth condition has been used.
                // (Hypothesis: with parms as doubles the first condition is used and only in very special cases the second.
                //  Addendum: This hypothesis holds.)
                return A > B ? -1 : 1;
                }
            }
        ++pcnt[NPARMS >> 2];
        }
    return 0;
    }
#endif

static int nparms = 0;

static int comp_parms0_off(const vertex * a,const vertex * b)
    {

    int off = minparmsoff;
    if(off < parmsoff)
        off = parmsoff;
    if(parmsoff >= nparms)
        {
        fprintf(stderr,"parmsoff is set to high. There are not enough rows of parameters. Fix in graph.cpp\n");
        exit(-1);
        }
    for(int o = off;o < nparms;o += ROWPARMS)
        {
            //double A = parms.Matrix[o]*a->R__R + parms.Matrix[o+1]*a->W__R + parms.Matrix[o+2]*a->R__W + parms.Matrix[o+3]*a->W__W;
            //double B = parms.Matrix[o]*b->R__R + parms.Matrix[o+1]*b->W__R + parms.Matrix[o+2]*b->R__W + parms.Matrix[o+3]*b->W__W;
            double A = parms.Matrix[o+3]*a->R__R + parms.Matrix[o+0]*a->W__R + parms.Matrix[o+1]*a->R__W + parms.Matrix[o+2]*a->W__W;
            double B = parms.Matrix[o+3]*b->R__R + parms.Matrix[o+0]*b->W__R + parms.Matrix[o+1]*b->R__W + parms.Matrix[o+2]*b->W__W;
        if(A != B)
            {
            return A > B ? -1 : 1;
            }
        }
    return 0;
    }

struct funcstruct
    {
    bool compute_parms;
    const char * number;
    const char * name;
    int (*comp)(const vertex * a,const vertex * b);
    };

static struct funcstruct funcstructs[] =
    {
        {true,"0","parms",comp_parms}, //makeaffix.exe mydata.txt 0 affixrules XX 123 parms
#if _NA
        {false,"1","fairly_good",comp_fairly_good},
        {false,"2","even_better",comp_even_better},
        {false,"3","affiksFEW3",comp_affiksFEW3},
        {false,"4","affiksFEW",comp_affiksFEW},
        {false,"5","affiksFEW2",comp_affiksFEW2},
        {false,"6","fixNA",comp_fixNA},
        {false,"7","fruit",comp_fruit},
        {false,"8","ice",comp_ice},
        {false,"9","pisang",comp_pisang},
        {false,"10","kiwi",comp_kiwi},
        {false,"11","carrot",comp_carrot},
        {false,"12","peen",comp_peen},
        {false,"13","beet",comp_beet},
        {false,"14","sugar",comp_sugar},
        {false,"15","affiksFEW2org",comp_affiksFEW2org},
#endif
        {false,"16","honey",comp_honey},
        {false,"17","koud",comp_koud},
        {false,"18","parms0",comp_parms0_off},
        {false,"19","parmsoff",comp_parms0_off},
        {false,0,0,0}
    };


bool setCompetitionFunction(const char * functionname,const char * extra,bool suffixonly,bool & compute_parms,const char * parmstxt)
    {
    int i;
    if(VERBOSE)
        printf("setCompetitionFunction(functionname:%s,extra:%s,suffixonly:%s,parmstxt:%s)"
	          , functionname
	          , extra
	          , suffixonly ? "true" : "false"
              , parmstxt ? parmstxt : "Not defined"
	          );

    size_t langlength = strlen(extra);
    const char * underscore = strchr(extra,'_');
    if(underscore != NULL)
        langlength = underscore - extra;

    for(i = 0;funcstructs[i].number;++i)
        if(!strcmp(functionname,funcstructs[i].number) || !strcmp(functionname,funcstructs[i].name))
            {
            comp = funcstructs[i].comp;
            compute_parms = funcstructs[i].compute_parms;
            if(comp == comp_parms0_off)
                {
                for(unsigned int j = 0;j < sizeof(bests)/sizeof(bests[0]);++j)
                    {
                    if(  bests[j].suffixonly == suffixonly 
                      && (langlength == strlen(bests[j].langbase)) // 20130125
                      && !strncmp(bests[j].langbase,extra,strlen(bests[j].langbase))
                      )
                        {
                        printf("bests[%d].suffixonly == [%s] bests[%d].langbase == [%s]\n",j,bests[j].suffixonly ? "true" : "false",j,bests[j].langbase);
                        printf("comp = comp_parms0_off\n");
                        comp = comp_parms0_off;
                        printf("bests[%d].rows == [%d]\n",j,bests[j].rows);
                        nparms = bests[j].rows * ROWPARMS;
                        if(nparms > NPARMS)
                            {
                            fprintf(stderr,"Too many rows of parameters in bestParms struct for %s (%d, max allowed %d)\n",extra,nparms,NPARMS);
                            exit(-1);
                            }
                        for(int k = 0;k < nparms;++k)
                            parms.Matrix[k] = bests[j].val[k];
                        if(parmstxt)
                            {
                            FILE * f = fopen(parmstxt,"w");
                            if(f)
                                {
                                fprintf(f,"bests[%d].suffixonly == [%s]\nbests[%d].langbase == [%s]\n",j,bests[j].suffixonly ? "true" : "false",j,bests[j].langbase);
                                fprintf(f,"comp = comp_parms0_off\n");
                                fprintf(f,"bests[%d].rows == [%d]\n",j,bests[j].rows);
                                fprintf(f,"  R->R  W->R  R->W  W->W\n");
                                for(int k = 0;k < nparms;++k)
                                    {
                                    if(k % ROWPARMS == 0)
                                        fprintf(f,"\n");
                                    fprintf(f,"%6d",parms.Matrix[k]);
                                    }
                                fprintf(f,"\n");
                                fclose(f);
                                }
                            }
                        break;
                        }
                    }
                if(nparms == 0)
                    {
                    fprintf(stderr,"No parameters defined for \"%s\"\nChoose one of:\n",extra);
                    for(unsigned int j = 0;j < sizeof(bests)/sizeof(bests[0]);++j)
                        {
                        fprintf(stderr,"\t%s %s\n",bests[j].langbase,bests[j].suffixonly ? "suffix":"affix");
                        }
                    fprintf(stderr,"Or find optimal parameters for %s and put these in comp.cpp.\n",extra);
                    getchar();
                    exit(-1);
                    }
                if(VERBOSE)
                    {
                    printf("comp_parms0_off\n");
                    }
                }
            for(i = 0;i < nparms; ++i)
                {
                if(parms.Matrix[i])
                    {
                    minparmsoff = i / ROWPARMS;
                    minparmsoff *= ROWPARMS;
                    break;
                    }
                }
            if(VERBOSE)
                printf("minparmsoff = %d \n",minparmsoff);
            return true;
            }
    return false;
    }

