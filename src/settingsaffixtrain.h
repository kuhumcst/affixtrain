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

#ifndef SETTINGSAFFIXTRAIN_H
#define SETTINGSAFFIXTRAIN_H

// turn off warnings against dangerous standard functions:
#define _CRT_SECURE_NO_DEPRECATE

// OS-specific stuff:
#ifdef WIN32
#define COPY "copy "
#define MKDIR "mkdir "
#define DIRSEP '\\'
#define ZD "%Iu"
// %zd not working in Microsoft C, so use %Id instead. 
// See https://msdn.microsoft.com/en-us/library/tcxf1dw6.aspx

#else
#define COPY "cp -f "
#define MKDIR "mkdir "
#define DIRSEP '/'
#define ZD "%zd"
#endif

#define DOIMPEDANCE 0

// Full form - lemma pairs that require unique rules are discarded if
// PRUNETRAININGPAIRS 1
#define PRUNETRAININGPAIRS 0

#define SMALLMEMORY 0
#if SMALLMEMORY
#define MAXPAIRS 8000000 //1000000 //500000 //100000 //500000 
#endif
// Set this number as high as possible. A high value entails very heavy use 
// of memory!

#define MINLINES 50000.0 
// Start determination of optimal parameters with minimally MINLINES lines (double)

#define AMBIGUOUS 1
// AMBIGUOUS 1
// Handle ambiguous data by redoing training a second time without the
// homographs that were handled in the first round.
// AMBIGUOUS 0
// Disambiguate the training data, either by lemma frequency or heuristically by taking the first

#define PESSIMISTIC 0

#define RULESASTEXTINDENTED 0
/* RULESASTEXTINDENTED 1
produce human readable tree:

threshold 0
|    {^*$    ^*$}
        R+    Buffalo    Buffalo,semi-metro    semi-metro,pornofoto    pornofoto, ...
 |    {^*n$    ^*n$}
        R+    mandarijn    mandarijn,waslijn    waslijn, ...
  |    {^*en$    ^*en$}
        R+    overhebben    overhebben
   |    {^*gen$    ^*g$}
        R+    boterbergen    boterberg
    |    {^*ngen$    ^*ng$}
        R+    lammergangen    lammergang
     |    {^*ingen$    ^*ing$}
        R+    plaatspanningen    plaatspanning,stukkenrekeningen    stukkenrekening, ...
      |    {^*elingen$    ^*eling$}
        R+    omwisselingen    omwisseling,wildelingen    wildeling,loonregelingen    loonregeling
...
 |    {^k*u$    ^k*uen$}
        R    keu    keuen
threshold 0:3309 words 1553 nodes 1322 nodes with words
*/

#define BRACMATOUTPUT 1
/* BRACMATOUTPUT 1

affixtrain -b binaryFlexRules -t myprettyrules

produces output myprettyrules.bra as a Bracmat expression:

((((=?W: ?A ).(=!A)).1)
,((((=?W n&@(!W: ?A )).(=!A n)).2)
,((((=?W e&@(!W: ?A )).(=!A en)).3)
,((((=?W g&@(!W: ?A )).(=!A g)).4)
,((((=?W n&@(!W: ?A )).(=!A ng)).5)
,((((=?W i&@(!W: ?A )).(=!A ing)).6)
,((((=?W el&@(!W: ?A )).(=!A eling)).7)

...

(((=k ?W u&@(!W: ?A )).(=k !A uen)).1553)
)



Bracmat function to lemmatise a word, given the structure stored in 'brafile'
(see declaration after this comment)


  ( lemmatise
  =     apply tree woord pat rep subtree lemma nr rule
      , A B C D E F W prune
      , first second one none two ntwo sub
    .   !arg:(?tree.?woord)
      & ( apply
        =   pat,rep,woord
          .   !arg:(?woord.(=?pat).?rep)
            & @(!woord:!pat)
            & !rep
        )
      & ( prune
        =   AA AAA BBB CCC MM ZZ
          .     whl
              ' (   !arg
                  : ?AA (?AAA.?BBB) ?MM (!AAA.?CCC) ?ZZ
                & !AA (!AAA.!BBB !CCC) !MM !ZZ:?arg
                )
            & !arg
        )
      & ( !tree:&(.)
        |   !tree:((?rule.#?nr),?subtree) ?tree
          & (   apply$(!woord.!rule):(=?lemma)
              & (   lemmatise$(!subtree.!W):?sub
                  & prune$!sub
                | (str$!lemma.!nr)
                )
            | !tree:~&lemmatise$(!tree.!woord)
            | ~`
            )
        |   !tree:(?.?):(?first.~#?second)
          & lemmatise$(!first.!woord):?first
          & lemmatise$(!second.!woord):?second
          & prune$(!first !second)
        |   !tree:((?.#):(?rule.#?nr)) ?tree
          & (   apply$(!woord.!rule):(=?lemma)
              & (str$!lemma.!nr)
            | !tree:~&lemmatise$(!tree.!woord)
            )
        )
  )
& get$"dict4ambi.flexrules.pass4.cutoff0.accumulated.txt.bra":?tree
& ( tree
  =   (((=?W:?A).(=!A)).1)
    , ( ( (   ( ( (=?W k&@(!W:?A))
                . (=!A kken)
                )
              . 2
              )
              ( ( (=?W n&@(!W:?A))
                . (=!A n)
                )
              . 3
              )
              ( ((=?W:?A o ?B).(=!A !B en))
              . 4
              )
              ( ( (=?W l&@(!W:?A))
                . (=!A l)
                )
              . 5
              )
              ( ( (=?W a&@(!W:))
                . (=a)
                )
              . 6
              )
              ( ( (=?W b&@(!W:?A))
                . (=!A b)
                )
              . 7
              )
              ( ( (=?W c&@(!W:?A))
                . (=!A c)
                )
              . 8
              )
              ( ( (=?W d&@(!W:?A))
                . (=!A)
                )
              . 9
              )
          .   ( ( (=?W k&@(!W:?A))
                . (=!A kken)
                )
              . 10
              )
              ( ( (=?W n&@(!W:?A))
                . (=!A n)
                )
              . 11
              )
              ( ((=?W:?A o ?B).(=!A !B en))
              . 12
              )
              ( ( (=a ?W&@(!W:?A))
                . (=a !A)
                )
              . 13
              )
              ( ( (=?W l&@(!W:?A))
                . (=!A l)
                )
              . 14
              )
              ( ( (=?W d&@(!W:?A))
                . (=!A)
                )
              . 15
              )
          )
        .   ( ( ( (=?W d&@(!W:?A))
                . (=!A)
                )
              . 16
              )
            ,   ( ( (=?W c&@(!W:?A))
                  . (=!A)
                  )
                . 17
                )
                ( ( (=?W o&@(!W:?A))
                  . (=!A den)
                  )
                . 18
                )
            )
            ( ( (=?W k&@(!W:?A))
              . (=!A kken)
              )
            . 19
            )
            ( ( (=?W n&@(!W:?A))
              . (=!A n)
              )
            . 20
            )
            ( ((=?W:?A o ?B).(=!A !B en))
            . 21
            )
        )
      .   ( ( ( (=?W d&@(!W:?A))
              . (=!A)
              )
            . 22
            )
          ,   ( ( (=?W o&@(!W:?A))
                . (=!A den)
                )
              . 23
              )
              ( ( (=?W bc&@(!W:?A))
                . (=!A)
                )
              . 24
              )
          )
          ( ( (=?W k&@(!W:?A))
            . (=!A kken)
            )
          . 25
          )
          ( ( (=?W n&@(!W:?A))
            . (=!A n)
            )
          . 26
          )
          ( ( (=a ?W&@(!W:?A))
            . (=a !A)
            )
          . 27
          )
          ( ((=?W:?A o ?B).(=!A !B en))
          . 28
          )
      )
  )
&     abcd
      abc
      ab
      a
      kald
      kal
      koop
      kopen
      bak
      bakken
      dood
      doden
  : ?testWords
&   whl
  ' ( !testWords:%?word ?testWords
    & out$(!word ":" lemmatise$(!tree.!word))
    )
& ;

    How to lemmatise the word 'acculadestationerne':

      get$"pretty.bra":?tree
    & lemmatise$(!tree.acculadestationerne):(?lemma.?rulenr) (|(?lemma2.?rulenr2))

    Go to https://github.com/BartJongejan/Bracmat.

*/

#define RULESASTEXT 0
/* RULESASTEXT 1
Produce output as text. This format is close to the binary format, only made 
"readable".

0            
6416        n    n


6260        e    en

500        g    g


132        n    ng

84        i    ing
56        el    eling

...

0k    k    u    uen
*/

#define WRITEINVERTED 0
/* WRITEINVERTED 1
Reverse words in training set and write them to file. You can base another
training on this set and see whether results get better.
Background: words are scanned from left to right. If an "infix" occurs more
than once in a word, the replacement may land in the wrong place. There may
be a bias indicating that "infixes" tend to be near the end, rather than near
the start, of a word. In that case, reversing words may improve results.
*/

#define _NA 0
/* NA 1
Allow weight functions (see comp.cpp) to use the numbers W__NA and N__NA,
the number of wrongly (correctly) lemmatized words that a candidate rule 
does not apply to. These numbers are not updated (=decremented) as siblings
are added, in contrast to W__W, W__R, R__W and R__R.
20150309: now they are updated!
*/

#define LEMMACLASS 0
/* LEMMACLASS 1
Read lemma class from file. No code in the training program currently uses
the lemma class.
*/

#define LEMMAINL 0
/* LEMMAINL 1
Read lemma frequency from file. No code in the training program currently uses
the lemma frequency.
Also enables word frequency.
*/

#define WORDCLASS 0
/* WORDCLASS 1
Read word class from file. 
*/

#if defined __BORLANDC__
#define CONSTSTRCHR
#define REFER(v) v;
#elif defined _MSC_VER
#define REFER(v) v;
#else
#define REFER(v)
#endif

#define ANY   4 //':'
#define START 5 //'^'
#define END   6 //'$'

#if defined _WIN64 || defined _WIN32
/*Microsoft*/
#define LONG    long long
#define LONGU   "%llu"
#define LONGD   "%lld"
#define LONG0D  "%0lld"
#define LONGX   "%llX"
#define STRTOUL _strtoui64
#define STRTOL  _strtoi64
#define FSEEK   _fseeki64
#define FTELL   _ftelli64
// \r in printf
#define STARTLINE 13
#else
#define LONG    long
#define LONGU   "%lu"
#define LONGD   "%ld"
#define LONG0D  "%0ld"
#define LONGX   "%lX"
#define STRTOUL strtoul
#define STRTOL  strtol
#define FSEEK   fseek
#define FTELL   ftell
// \n in printf
#define STARTLINE 10
#endif

typedef enum {failure,wrong,right} matchResult;

#endif
