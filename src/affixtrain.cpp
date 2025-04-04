/*
AFFIXTRAIN - supervised learning of affix rules for CSTLEMMA

Copyright (C) 2014  Center for Sprogteknologi, University of Copenhagen

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

#define VERSION "3.94"

#include "affixtrain.h"
#include "testrules.h"
#include "flexcombi.h"
#include "optionaff.h"
#include "applyaffrules.h"
#include "node.h"
#include "comp.h"
#include "trainingpair.h"
#include "ruletemplate.h"
#include "vertex.h"
#include "hashtable.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "utf8func.h"


static const char* SCUT = ".cutoff";
static const char* emptyRules = "\rV3\r\0\0\0\0\t\t\t\n";

#define CHECK(a)

int openfiles = 0;

static const char* tempFolder(const char* filename, optionStruct* options)
    {
    static char* fullname = NULL;
    size_t length;
    int written;
    CHECK("AglobTempDir");
    if(options->tempDir())
        {
        if(fullname)
            delete[] fullname;
        length = strlen(options->tempDir()) + strlen(filename) + 3;
        fullname = new char[length];
        written = sprintf(fullname, "%s%c%s", options->tempDir(), DIRSEP, filename);
        if((size_t)written >= length)
            {
            fprintf(stderr, "tempDir: Buffer overrun");
            exit(-1);
            }
        return fullname;
        }
    CHECK("BglobTempDir");
    return filename;
    }

union pointers
    {
    unsigned char* uchars;
    char* chars;
    const char* cchars;
    const short int* shorts;
    const long int* longs;
    };

class tagClass
    {
    public:
        char* name;
        tagClass* next;
        tagClass(const char* name, size_t length, tagClass* next)
            {
            this->name = new char[length + 1];
            strncpy(this->name, name, length);
            this->name[length] = '\0';
            this->next = next;
            }
        ~tagClass()
            {
            delete[] name;
            delete next;
            }
        bool has(const char* str, size_t length)
            {
            return  ((length == strlen(name))
                     && !strncmp(str, name, length)
                     )
                || (next && next->has(str, length));
            }
    };


class shortRulePair;
FILE* fopenOrExit(const char* name, const char* mode, const char* descriptionOfFile) // must not contain format specifiers
    {
    assert(name);
    assert(mode);
    assert(descriptionOfFile);
    FILE* fp = fopen(name, mode);
    ++openfiles;
    if(!fp)
        {
        fprintf
        (stderr
         , "Open files: %d\n%s: Cannot open file %s for %s (%s). Exiting.\n"
         , openfiles
         , descriptionOfFile
         , name
         , strchr(mode, 'r') ? "reading" : strchr(mode, 'a') ? "appending" : "writing"
         , mode
        );
        exit(-1);
        }
    return fp;
    }

class ruleTemplate;






const char* find(const char* string
                 , int c
                 , const char* end
)
    {
    while(string < end && *string != c)
        {
        ++string;
        }
    return (string < end) ? string : NULL;
    }

struct aFile
    {
    union pointers file;
    union pointers* Lines;
    const char* eob;
    char* fname;
    const char* columns;
    size_t size;
    size_t lines;
    int filesize;
    aFile(const char* fname, optionStruct* options, const char* columns) :Lines(NULL), eob(NULL), columns(columns), lines(0), filesize(0)
        {
        assert(fname);

        FILE* fp = fopenOrExit(fname, "rb", "Input file");
        options->info(5, "reading file %s ...", fname);
        CHECK("CglobTempDir");
        this->fname = new char[strlen(fname) + 1];
        strcpy(this->fname, fname);
        file.chars = NULL;
        size = 0;
        while(getc(fp) != EOF)
            {
            ++size;
            }
        file.chars = new char[size + 1];
        rewind(fp);
        size = 0;
        int kar;
        while((kar = getc(fp)) != EOF)
            {
            file.chars[size++] = (char)kar;
            }
        file.chars[size] = '\0';
        --openfiles;
        fclose(fp);

        eob = file.chars + size;

        const char* buf;
        size_t line;
        const char* nl;
        const char* tab;
        for(buf = file.chars, line = 0; buf < eob;)
            {
            for(;;)
                {
                nl = find(buf, '\n', eob);
                if(nl)
                    {
                    tab = find(buf, '\t', nl);
                    if(tab)
                        ++line;
                    else
                        {
                        buf = nl + 1;
                        continue;
                        }
                    }
                break;
                }
            if(nl)
                {
                buf = nl + 1;
                }
            else
                {
                tab = find(buf, '\t', eob);
                if(tab)
                    ++line;
                buf = eob;
                }
            }
        options->info(5, "line " ZD "\n", line);

        lines = line;
        Lines = new pointers[lines + 1];

        for(buf = file.cchars, line = 0; buf < eob;)
            {
            for(;;)
                {
                nl = find(buf, '\n', eob);
                if(nl)
                    {
                    tab = find(buf, '\t', nl);
                    if(tab)
                        Lines[line++].cchars = buf;
                    else
                        {
                        buf = nl + 1;
                        continue;
                        }
                    }
                break;
                }
            if(nl)
                {
                buf = nl + 1;
                }
            else
                {
                tab = find(buf, '\t', eob);
                if(tab)
                    Lines[line++].cchars = buf;
                buf = eob;
                }
            }
        options->info(5, "line " ZD " read\n%s:", line, fname);
        }
    ~aFile()
        {
        delete[] file.chars;
        delete[] Lines;
        delete[] this->fname;
        }
    };

static const int unequalBit = 1 << 0;
static const int equalBit = 1 << 1;
static const int transition = unequalBit | equalBit;



static bool isSpace(int a)
    {
    CHECK("RglobTempDir");
    switch(a)
        {
        case ' ':
        case '\t':
        case '\r':
        case '\n':
            return true;
        default:
            return false;
        }
    }

static tagClass* collectTags(optionStruct* options)
// "123456" means Word, Lemma, Wordfreq, Lemmafreq, Wordclass, Lemmaclass
    {
    struct aFile afile(options->wordLemmaList(), options, NULL);

    CHECK("SglobTempDir");
    tagClass* Tags = NULL;
    size_t line;
    pointers* Lines = afile.Lines;
    options->info(5, "Checking " ZD " lines in %s for tags %s\n", afile.lines, afile.fname, options->columns());
    for(line = 0; line < afile.lines; ++line)
        {
        const char* cols[18];
        size_t lengths[18];
        const char* limit = (line < afile.lines - 1) ? Lines[line + 1].cchars : afile.eob;

        const char* q = Lines[line].cchars;
        cols[0] = q;
        unsigned int ii = 0;
        while((q = find(q, '\t', limit)) != NULL
              && (ii < sizeof(cols) / sizeof(cols[0]) - 1)
              )
            {
            lengths[ii] = (size_t)(q - cols[ii]);
            cols[++ii] = ++q;
            }
        lengths[ii] = (size_t)(limit - cols[ii] - 1);
        unsigned int maxii = ++ii;
        cols[maxii] = q ? q : limit;
        const char* column = options->columns();
        if(column == 0)
            column = "FBT";
        for(ii = 0
            ; *column && (ii < maxii)
            ; ++column, ++ii
            )
            {
            switch(*column)
                {
                case '3':
                case 'T':
                case 't':
                    {
                    const char* C = cols[ii];
                    size_t L = lengths[ii];
                    // Trim CR, white space and other obnoxious stuff.
                    while(L > 0 && 0 < C[L - 1] && C[L - 1] <= ' ')
                        --L;
                    if(!Tags || !(Tags->has(C, L)))
                        {
                        Tags = new tagClass(C, L, Tags);
                        }
                    break;
                    }
                default:
                    ;
                }
            }
        }
    if(Tags == 0)
        {
        const char* cols = options->columns();
        if(cols == 0 || strchr(cols, '3') || strchr(cols, 't') || strchr(cols, 'T'))
            {
            options->setn("FB"); // No column with tags found.
            options->info(0, "No column with tags found in %s.\n"
                          "Changed the column format (-n option) to FB.\n",
                          options->wordLemmaList());
            }
        }
    else
        {
        if(options->columns() == 0)
            {
            options->setn("FBT"); // Column with tags found.
            options->info(0, "Column with tags found in %s.\n"
                          "I suppose the column format (-n option) is FBT.\n"
                          "If you want to ignore the tags, define -n FBO.\n"
                          "(Or whatever the actual column order is.)\n",
                          options->wordLemmaList());
            }
        }
    return Tags;
    }

static int compare(const void* arg1, const void* arg2)
    {
    CHECK("UglobTempDir");
    const trainingPair* A = *(const trainingPair* const*)arg1;
    const trainingPair* B = *(const trainingPair* const*)arg2;
    int ret = A->cmpWord(B);
    if(!ret)
        {
        ret = A->cmpLemma(B);
#if WORDCLASS
        if(!ret)
            {
            ret = A->cmpWordClass(B);
            }
#endif
#if LEMMACLASS
        if(!ret)
            {
            ret = A->cmpLemmaClass(B);
            }
#endif
        }
    return ret;
    }

static int markAmbiguous(size_t allPairs, trainingPair* TrainingPair, optionStruct* options)
    {
    options->info(2, "markAmbiguous\n");
    trainingPair** pTrainingPair = new trainingPair * [allPairs];
    size_t j;
    for(j = 0; j < allPairs; ++j)
        pTrainingPair[j] = TrainingPair + j;
    qsort((void*)pTrainingPair, allPairs, sizeof(trainingPair*), compare);

    for(j = 0; j < allPairs; ++j)
        // Mark words containing space for skipping
        {
        size_t nn = pTrainingPair[j]->itsWordlength();
        const char* s;
        // do not use words containing space
        for(s = pTrainingPair[j]->itsWord(); nn-- > 0;)
            {
            if(s[nn] == ' ')
                {
                pTrainingPair[j]->set(b_skip);
                }
            }
        }
    int n = 0;



    for(j = 0; j < allPairs; ++j)
        {
        if(!pTrainingPair[j]->isset(b_skip))
            {
            size_t k;
            for(k = j + 1
                ;    k < allPairs
                //&& !pTrainingPair[j]->isset(b_skip)
                && !pTrainingPair[j]->cmpWord(pTrainingPair[k])
#if WORDCLASS
                && !pTrainingPair[j]->cmpWordClass(pTrainingPair[k])
#endif
                ; ++k
                )
                {
                if(!pTrainingPair[j]->cmpLemma(pTrainingPair[k]) /*|| pTrainingPair[j]->cmpLemmaClass(pTrainingPair[k])*/) //TODO Tags!
                    pTrainingPair[k]->set(b_doublet | b_skip);
                }
            }
        }



    for(j = 0; j < allPairs;// - 1;
#if AMBIGUOUS
#else
        ++j
#endif
        )
        {
        if(pTrainingPair[j]->isset(b_skip))
            ++j;
        else
            {
            size_t k;
#if AMBIGUOUS
            trainingPair* pAlt = pTrainingPair[j];
#endif
            for(k = j + 1
                ;    k < allPairs
                //&& !pTrainingPair[j]->isset(b_skip)
                && !pTrainingPair[j]->cmpWord(pTrainingPair[k])
#if WORDCLASS
                && !pTrainingPair[j]->cmpWordClass(pTrainingPair[k])
#endif
                ; ++k
                )
                {
                if(pTrainingPair[j]->cmpLemma(pTrainingPair[k]) /*|| pTrainingPair[j]->cmpLemmaClass(pTrainingPair[k])*/) //TODO Tags!
                    {
#if AMBIGUOUS
                    if(!pTrainingPair[k]->isset(b_skip))
                        {
                        pTrainingPair[k]->set(b_ambiguous);
                        pAlt->Alt = pTrainingPair[k];
                        pAlt = pTrainingPair[k];
                        ++n;
                        }
#else
#if LEMMAINL
                    switch(pTrainingPair[j]->cmpFreq(pTrainingPair[k]))
                        {
                        case 1:
                            pTrainingPair[k]->set(b_skip);
                            break;
                        case -1:
                            pTrainingPair[j]->set(b_skip);
                            break;
                        default:
#else
                            {
#endif
                            /* Take the first one */
                            if(!pTrainingPair[k]->isset(b_ambiguous | b_skip))
                                {
                                pTrainingPair[k]->set(b_ambiguous);
                                if(famb)
                                    pTrainingPair[k]->fprint(famb);
                                ++n;
                                }
                            }
#endif
                    }
                else
                    pTrainingPair[k]->set(b_doublet | b_skip);
                }
#if AMBIGUOUS
            // close the ring of homographs
            if(pAlt != pTrainingPair[j])
                {
                pTrainingPair[j]->set(b_ambiguous); // ALL members in the ring are marked b_ambiguous!
                pAlt->Alt = pTrainingPair[j];
                }
            j = k;
#else
            //        if(allFile && j < allPairs)
            //          pTrainingPair[j]->fprintAll(allFile);
#endif
            }
        }

    delete[] pTrainingPair; // Bart 20081008
    options->info(2, "markAmbiguous DONE. Returns %d\n", n);
    return n;
    }

#if PESSIMISTIC
static int compare2(const void* arg1, const void* arg2)
    {
    const trainingPair* A = *(const trainingPair* const*)arg1;
    const trainingPair* B = *(const trainingPair* const*)arg2;
    int ret = A->cmpLemma(B);
    if(!ret)
        {
#if LEMMACLASS
        ret = A->cmpLemmaClass(B);
        if(!ret)
#endif
            {
            ret = A->cmpWord(B);
            }
        }
    return ret;
    }

static int markParadigms(size_t allPairs, trainingPair * TrainingPair, FILE * fparadigms)
    {
    //    FILE * allFile = fopenwb("allFile.txt");
    trainingPair** pTrainingPair = new trainingPair * [allPairs];
    size_t j;
    for(j = 0; j < allPairs; ++j)
        pTrainingPair[j] = TrainingPair + j;
    qsort((void*)pTrainingPair, allPairs, sizeof(trainingPair*), compare2);

    int n = 0;
    for(j = 0; j < allPairs/* - 1*/;)
        {
        size_t k;
        trainingPair* pAlt = pTrainingPair[j];
        if(fparadigms)
            {
            if(!pTrainingPair[j]->isset(b_skip))
                pTrainingPair[j]->fprint(fparadigms);
            }
        for(k = j + 1
            ;    k < allPairs
            && !pTrainingPair[j]->isset(b_skip)
            && !pTrainingPair[j]->cmpLemma(pTrainingPair[k])
            ; ++k
            )
            {
            if(!pTrainingPair[k]->isset(b_skip))
                {
                pAlt->AltLemma = pTrainingPair[k];
                pAlt = pTrainingPair[k];
                if(fparadigms)
                    {
                    fprintf(fparadigms, "PARADIGM:");
                    pTrainingPair[k]->fprint(fparadigms);
                    }
                }
            }
        if(pAlt != pTrainingPair[j])
            {
            pAlt->AltLemma = pTrainingPair[j];
            }
        j = k;
        }
    delete[] pTrainingPair; // Bart 20081008
    return n;
    }
#endif

static trainingPair* globTrainingPair;

static trainingPair* readTrainingPairs(aFile & afile, size_t & pairs, optionStruct * options)
    {
    options->info(2, "readTrainingPairs\n");
    const char* tag = options->POStag();
    if(afile.lines == 0)
        return 0;
    trainingPair* TrainingPair = new trainingPair[afile.lines];
    globTrainingPair = TrainingPair;
    CHECK("TglobTempDir");
    pairs = 0;
    size_t line;
//    size_t taglength = tag ? strlen(tag) : 0;
    pointers* Lines = afile.Lines;
    options->info(2, "readLines with tag %s\n", (tag && *tag) ? tag : "_ (UNDEFINED)");
    for(line = 0; line < afile.lines; ++line)
        {
        const char* cols[18];
        const char* Word = NULL;
        size_t lengths[18];
        size_t wordlength = 0;
#if WORDCLASS
        const char* WordClass = NULL; // unused
        size_t wordclasslength = 0; // unused
#endif
        const char* LemmaHead = NULL;
        size_t lemmalength = 0;
#if LEMMACLASS
        const char* LemmaClass = NULL;
        size_t lemmaclasslength = 0;
#endif
#if LEMMAINL
        long Inl = 0;
        long Lemma_Inl = 0;
#endif
        const char* limit = (line < afile.lines - 1) ? Lines[line + 1].cchars : afile.eob;

        const char* q = Lines[line].cchars;
        cols[0] = q;
        unsigned int ii = 0;
        while((q = find(q, '\t', limit)) != NULL
              && (ii < sizeof(cols) / sizeof(cols[0]) - 1)
              )
            {
            lengths[ii] = (size_t)(q - cols[ii]);
            cols[++ii] = ++q;
            }
        lengths[ii] = (size_t)(limit - cols[ii] - 1);
        unsigned int maxii = ++ii;
        cols[maxii] = q ? q : limit;
        const char* column;
        // "-n 123456" means Word, Lemma, Wordfreq, Lemmafreq, Wordclass, Lemmaclass
        for(column = afile.columns, ii = 0
            ;    *column
            && (ii < maxii)
            ; ++column, ++ii
            )
            {
            const char* C = cols[ii];
            size_t L = lengths[ii];
            // Trim CR, white space and other obnoxious stuff.
            while(L > 0 && 0 < C[L - 1] && C[L - 1] <= ' ')
                --L;
            switch(*column)
                {
                case '1':
                case 'F':
                case 'f':
                case 'W':
                case 'w':
                    Word = C;
                    wordlength = L;
                    break;
                case '2':
                case 'B':
                case 'b':
                case 'L':
                case 'l':
                    LemmaHead = C;
                    lemmalength = L;
                    break;
                case '3':
                case 'T':
                case 't':
                    /*
                    if(tag && (taglength == (size_t)L) && !strncmp(C, tag, taglength))
                        {
                        doUse = true;
                        }
                    break;
                    */
#if LEMMAINL
                case '3':
                    Inl = strtol(C, NULL, 10);
                    break;
                case '4':
                    Lemma_Inl = strtol(C, NULL, 10);
                    break;
#endif
#if WORDCLASS
                case '5':
                    WordClass = C;
                    wordclasslength = L;
                    break;
#endif
#if LEMMACLASS
                case '6':
                    LemmaClass = C;
                    lemmaclasslength = L;
                    break;
#endif
                default:
                    ;
                }
            }
        if(Word /*&& doUse*/)
            {
            //doUse = ((tag == NULL) || !*tag);
            while(wordlength > 0 && isSpace(Word[wordlength - 1]))
                --wordlength;
            while(lemmalength > 0 && isSpace(LemmaHead[lemmalength - 1]))
                --lemmalength;
#if WORDCLASS
            while(wordclasslength > 0 && isSpace(WordClass[wordclasslength - 1]))
                --wordclasslength;
#endif
#if LEMMACLASS
            while(lemmaclasslength > 0 && isSpace(LemmaClass[lemmaclasslength - 1]))
                --lemmaclasslength;
#endif
            if(wordlength > 0 && lemmalength > 0)
                {
                TrainingPair[pairs].init
                (Word
                 , wordlength
                 , LemmaHead
                 , lemmalength
#if WORDCLASS
                 , WordClass
                 , wordclasslength
#endif
#if LEMMACLASS
                 , LemmaClass
                 , lemmaclasslength
#endif
#if LEMMAINL
                 , Inl, Lemma_Inl
#endif
                );
                ++pairs;
                }
            else if(wordlength + lemmalength > 0)
                {
                fprintf(stderr, "Word [%.*s] and/or lemma [%.*s] is empty\n", (int)wordlength, Word, (int)lemmalength, LemmaHead);
                exit(-1);
                }
            }
        }
    if(tag && *tag)
        options->info(2, ZD " characters and " ZD " lines, " ZD " selected with tag %s       \n", afile.size, afile.lines, pairs, tag);
    else
        options->info(2, ZD " characters and " ZD " lines%c", afile.size, afile.lines, STARTLINE);
    options->info(2, "readTrainingPairs DONE\n");
    return TrainingPair;
    }

void FilterTagLines(const char* fname, optionStruct * options, const char* tag)
    {
    assert(fname);
    assert(tag);

    size_t tabcount;
    size_t tagtabcount;
    const char* column;

    for(column = options->columns(), tagtabcount = 0
        ;    *column
        ; ++column
        )
        {
        switch(*column)
            {
            case '3':
            case 'T':
            case 't':
                break;
            default:
                ++tagtabcount;
                continue;
            }
        break;
        }
    char* tagfilename = new char[strlen(fname) + strlen(tag) + 2];
    sprintf(tagfilename, "%s_%s", fname, tag);

    options->seti(tagfilename);
    options->setk(tag);
    options->setArgstring();
    options->seto(0);

    FILE* fp = fopenOrExit(fname, "rb", "Input file");
    FILE* fpt = fopenOrExit(tagfilename, "wb", "Single tag input file");

    char line[10000];
    int i = 0;
    int kar;
    bool empty = true;
    unsigned int lines = 0;
    while((kar = getc(fp)) != EOF)
        {
        line[i++] = (char)kar;
        if(kar == '\n')
            {
            while(--i >= 0 && 0 <= line[i] && line[i] <= ' ')
                {
                line[i] = 0;
                }
            if(i <= 0) // empty line
                {
                if(!empty)
                    {
                    fputc('\n', fpt);
                    empty = true;
                    }
                }
            tabcount = 0;
            for(i = 0; line[i] && tabcount < tagtabcount;)
                {
                if(line[i] == '\t')
                    {
                    ++tabcount;
                    }
                ++i;
                }
            if(tabcount == tagtabcount)
                {
                const char* p = tag;
                while(*p && line[i] == *p)
                    {
                    ++i;
                    ++p;
                    }
                if(*p == 0 && (line[i] <= ' '))
                    {
                    fwrite(line, strlen(line), 1, fpt);
                    fputc('\n', fpt);
                    empty = false;
                    ++lines;
                    }
                }
            i = 0;
            }
        }
    options->setBlobs(1);
    options->setLines(lines);
    fclose(fp);
    fclose(fpt);
    }

static void markTheAmbiguousPairs(trainingPair * TrainingPair, size_t pairs, optionStruct * options)
    {
    markAmbiguous(pairs, TrainingPair, options);

#if PESSIMISTIC
    char* filename = "paradigms.txt";
    FILE* fparadigms = fopenOrExit(filename, "wb", "fparadigms");
    markParadigms(pairs, TrainingPair, fparadigms);
    if(fparadigms)
        {
        --openfiles;
        fclose(fparadigms);
        }
#endif
    }

/* Create binary output */
static void rearrange(const char* filename
                      , FILE * folel
#if RULESASTEXT
                      , FILE * foleltxt
#endif
)
    {
    CHECK("fglobTempDir");
    FILE* fo = fopenOrExit(filename, "wb", "rearrange output");
    long end;
    end = ftell(folel);
    char* buf = new char[(size_t)end + 1]; // contents of folel, a textual file.
    rewind(folel);
    if(fread(buf, (size_t)end, 1, folel) != 1)
        return;
    buf[end] = '\0'; // 20140224 new
    unsigned int n = 0;
    int i;
    for(i = 0; i < end; ++i)
        if(buf[i] == '\n')
            ++n;
    char** pbuf = new char* [n + 1]; // lines
    unsigned int* size = new unsigned int[n + 1]; // binary size taken up by each rule
    unsigned int* cumsize = new unsigned int[n + 1]; // idem, cumulative
    unsigned int* ind = new unsigned int[n + 1]; // nesting levels, also used to store cumulative sizes
    pbuf[0] = buf + 0;
    n = 0;
    bool doind = true;  // first field on line is an indentation (or nesting
    // level) number
    ind[0] = 0;
    cumsize[0] = 0;
    size[0] = 0;
    for(i = 0; i < end; ++i)
        {
        if(buf[i] == '\n')
            {
            size[n] += sizeof(int) + sizeof(int); // room for index and '\0' byte
            size[n] >>= 2;
            size[n] <<= 2; // rounded to nearest word boundary
            pbuf[++n] = buf + i + 1; // start of next line
            cumsize[n] = size[n - 1];
            cumsize[n] += cumsize[n - 1];
            size[n] = 0;
            doind = true; // read first data on line as nesting level
            ind[n] = 0; // initialize nesting level to 0
            }
        else if(buf[i] == '\t')
            {
            if(doind)
                doind = false;
            else
                size[n]++;
            }
        else if(doind)
            {// read nesting level
            ind[n] *= 10;
            ind[n] += buf[i] - '0';
            }
        else
            size[n]++;
        }
    pbuf[n] = NULL;
    unsigned int lev[50];
    unsigned int j;
    for(j = 0; j < sizeof(lev) / sizeof(lev[0]); ++j)
        lev[j] = 0;
    for(j = 0; j < n; ++j)
        {
        unsigned int oldj = lev[ind[j]]; // find previous sibling
        if(oldj)
            ind[oldj] = cumsize[j]; // tell previous sibling 
                                    // where its next sibling is
        for(unsigned int k = ind[j] + 1; k < 50 && lev[k]; ++k)
            {
            lev[k] = 0; // forget about previous sibling'c children
            }
        lev[ind[j]] = j; // update who's the current node at this level.
        ind[j] = 0; // the level information is not needed anymore. Space can 
        // be reused to keep the position to the next sibling 
        // (if any).
        }
    char* p = buf;
    for(j = 0; j < n; ++j)
        {
        long pos = ftell(fo);
#ifdef _DEBUG
        long ppos = cumsize[j];
        assert(pos == ppos);
#endif
        if(ind[j] >= (unsigned int)pos)
            ind[j] -= (unsigned int)pos; // We only need to know how far to jump from here.
        fwrite(ind + j, sizeof(unsigned int), 1, fo);

#if RULESASTEXT
        fprintf(foleltxt, "%d", ind[j]);
#endif
        unsigned int written = sizeof(int);

        while(*p && *p != '\t')
            ++p;
        assert(*p);
        ++p;
        while(*p != '\n')
            { // write the pattern and the replacement (which are intertwined)
#if RULESASTEXT
            fputc(*p, foleltxt);
#endif
            fputc(*p++, fo);
            written++;
            }
#if RULESASTEXT
        fputc(*p, foleltxt);
#endif
        fputc(*p++, fo); // write the end-of-line marker
        written++;
        assert(written <= size[j]);
        assert(written > size[j] - sizeof(int));
        while(written < size[j]) // write until word boundary hit
            {
#if RULESASTEXT
            fputc('\n', foleltxt);
#endif
            fputc('\n', fo);
            written++;
            }
        }
    --openfiles;
    fclose(fo);
    delete[] buf;
    delete[] size;
    delete[] cumsize;
    delete[] ind;
    }

//int Nnodes = 0;

static void countNodes(node * tree, countAndWeight * Count, optionStruct * options)
    {
    Count->setNnodes(0);
    Count->setNnodes(tree->count());
    switch(options->getWeightFunction())
        {
        case econstant:
            break;
        case edepth:
            Count->setCountByDepth(0);
            Count->setCountByDepth(tree->countByDepth(1));
            break;
        case esupport:
            Count->setWeight(0.0);
            Count->setWeight(tree->weightedcount());
            break;
        case eentropy:
            Count->setWeight(0.0);
            Count->setWeight(tree->entropy(Count->getNnodes()) / log((double)(Count->getNnodes())));
            break;
        case esize:
            Count->setCountBySize(0);
            Count->setCountBySize(tree->countBySize());
            break;
        }
    }

static const char* flexRuleFileName(const char* ext, int threshold, const char* nflexrules)
    {
    if(nflexrules)
        return nflexrules;
    else
        {
        static char* filename = NULL;
        delete[] filename;
        size_t size = strlen("rules_.lem") + strlen(ext) + 20;
        filename = new char[size];
        if((int)size <= sprintf(filename, "rules_%d%s.lem", threshold, ext))
            {
            fprintf(stderr, "flexRuleFileName: filename too small");
            exit(-1);
            }
        return filename;
        }
    }

static bool writeRules(node * tree, const char* ext, int threshold, const char* nflexrules, optionStruct * options)
    {
    CHECK("gglobTempDir");
    char filename[1000];
    const char* FlexRuleFileName = flexRuleFileName(ext, threshold, nflexrules);
#if RULESASTEXTINDENTED
    sprintf(filename, "tree_%d%s.txt", threshold, ext);
    FILE* foo = fopenOrExit(filename, "wb", "indented rules");

    if(foo)
#endif
        {
#if RULESASTEXTINDENTED
        fprintf(foo, "threshold %d\n", threshold);
        int Nnodes = 0;
        int NnodesR = 0;
        int N = tree->print(foo, 0, Nnodes, NnodesR);
//        fprintf(foo,"threshold %d:%d words %d nodes %d weight %f nodes with words\n\n",threshold,N,Nnodes,NnodesR,weight);
        fprintf(foo, "threshold %d:%d words %d nodes %d nodes with words\n\n", threshold, N, Nnodes, NnodesR);
        --openfiles;
        fclose(foo);
        sprintf(filename, "rules_%d%s.txt", threshold, ext);
        foo = fopenOrExit(filename, "wb", "rules");
#endif
        if(256 <= sprintf(filename, "rules_%d%s.lel", threshold, ext))
            {
            fprintf(stderr, "writeRules: filename 2 small");
            exit(-1);
            }
        FILE* folel = fopenOrExit(tempFolder(filename, options), "wb+", "writeRules");
#if RULESASTEXT
        filename[strlen(filename) - 1] += 2; // change ".lel" to ".len"
        FILE* foleltxt = fopenOrExit(filename, "wb", "Text version");
        filename[strlen(filename) - 1] -= 2; // change ".len" back to ".lel"
#endif
        if(folel
#if RULESASTEXTINDENTED
           && foo
#endif
#if RULESASTEXT
           && foleltxt
#endif
           )
            {
#if RULESASTEXTINDENTED
            //fprintf(foo,"tree={%d %f}\n",Nnodes,weight); // "rules_%d%s.txt"
            fprintf(foo, "tree={%d}\n", Nnodes); // "rules_%d%s.txt"
#endif
            int nr = 0;
            strng L("");
            strng R("");
            options->info(6, "Calling printRules\n");
            printRules
            (tree
#if RULESASTEXTINDENTED
             , foo // "rules_%d%s.txt"
#endif
             , 0
             , folel // "rules_%d%s.lel"
             , 0
             , &L
             , &R
             , nr
             , options
            );
            options->info(6, "Calling printRules DONE\n");
#if RULESASTEXTINDENTED
            --openfiles;
            fclose(foo);
#endif
            options->info(6, "Calling rearrange\n");
            rearrange
            (FlexRuleFileName// the binary output of the training, 
            // third command line argument
             , folel     // .lel file, textual output with relative 
#if RULESASTEXT
             , foleltxt  // .len file, textual output with absolute 
             // positions, like in the binary output.
#endif
            );
            options->info(6, "Calling rearrange DONE\n");
            --openfiles;
            fclose(folel);
            if(remove(tempFolder(filename, options))) // del ".lel"
                {
                options->info(0, "cannot remove %s\n", filename);
                }
#if RULESASTEXT
            --openfiles;
            fclose(foleltxt);
#endif
            return true;
            }
        else
            {
#if RULESASTEXTINDENTED
            if(foo)
                {
                --openfiles;
                fclose(foo);
                }

#endif
            if(folel)
                {
                --openfiles;
                fclose(folel);
                }
#if RULESASTEXT
            if(foleltxt)
                {
                --openfiles;
                fclose(foleltxt);
                }
#endif
            }
        }
    return false;
    }

trainingPair* createChainOfPointersToTrainingPairs(trainingPair * TrainingPair, size_t allPairs)
    {
    trainingPair* train;
    trainingPair** ptrain = &train;
    for(size_t pairs = 0; pairs < allPairs; ++pairs)
        {
#if AMBIGUOUS
        if(!TrainingPair[pairs].isset(b_doublet | b_skip))
#else
        if(!TrainingPair[pairs].isset(b_ambiguous | b_doublet | b_skip))
#endif
            {
            *ptrain = TrainingPair + pairs;
            ptrain = &TrainingPair[pairs].Next;
            }
        }
    return train;
    }

trainingPair* saveChainOfPointersToTrainingPairsToFile(trainingPair * train, const char* fileName, const char* tagName)
    {
    FILE* fp;
    if(tagName && tagName[0])
        {
        char* name = new char[strlen(fileName) + strlen(tagName) + 2];
        sprintf(name, "%s_%s", fileName, tagName);
        fp = fopen(name, "w");
        delete[] name;
        }
    else
        {
        fp = fopen(fileName, "w");
        }
    trainingPair* ptrain = train;
    for(; ptrain; ptrain = ptrain->next())
        {
        ptrain->fprintTraining(fp);
        }
    fclose(fp);
    return train;
    }

static bool doTraining
(aFile & afile
 , const char* ext
 , int cutoff
 , const char* nflexrulesFormat
 , char* pairsToTrainInNextPassName
 , countAndWeight * Counts
 , optionStruct * options
)
    {
    size_t allPairs;
    trainingPair* TrainingPair = readTrainingPairs(afile, allPairs, options);
    if(TrainingPair == 0)
        return false;
    markTheAmbiguousPairs(TrainingPair, allPairs, options);
    trainingPair* train = createChainOfPointersToTrainingPairs(TrainingPair, allPairs);

    bool moreToDo = false;
    VertexPointerCount = 0;

    hashTable Hash(10);
    node* top;
    vertex ROOT(StartAnyEnd, StartAnyEnd);
    bool New;
    vertex* best = Hash.getVertex(&ROOT, New);

    options->info(2, "Going to build a decision tree\n");

    top = new node(best);
    options->info(2, "Top node is created\n");

    trainingPair* Right = NULL;
#if PRUNETRAININGPAIRS
    fpmourn = fopen("prunedTrainingPairs.txt", "w");
#endif
//    fprune = fopen("prunedTrainingPairs.txt","w");
//    train->printAll(fprune,"ALL\n",'\n');
//    fprintf(fprune,"____________\n");
    top->init(&Right, &train, 0, options);
    /* After top->init, the chain of pointers to training pairs, that starts
       with 'train', is broken and spread over the generated nodes in the
       decision tree. To traverse the training pairs, we have to use the
       table 'TrainingPair', as is done below. */
#if PRUNETRAININGPAIRS
    fclose(fpmourn);
#endif
    options->info(2, "Decision tree built\n");
    //    fclose(fprune);
    top = top->cleanup(NULL, 0);

    FILE* nexttrain = pairsToTrainInNextPassName ? fopenOrExit(tempFolder(pairsToTrainInNextPassName, options), "wb", "nexttrain") : NULL;
    if(nexttrain)
        {
        size_t pairs = 0;
        int donepairs = 0;
        for(pairs = 0; pairs < allPairs; ++pairs)
            {
#if AMBIGUOUS
            if(!TrainingPair[pairs].isset(b_doublet | b_skip))
#else
            if(!TrainingPair[pairs].isset(b_ambiguous | b_doublet | b_skip))
#endif
                {
                if(!TrainingPair[pairs].isset(b_ambiguous))
                    {
                    TrainingPair[pairs].printMore(nexttrain);
                    }
                else if(TrainingPair[pairs].getRes() != yes)
                    {
                    TrainingPair[pairs].printMore(nexttrain);
                    }
                else
                    {
                    donepairs++;
                    }
                }
            }

        if(donepairs)
            {
            moreToDo = true;
            }
        }

    if(nexttrain)
        {
        --openfiles;
        fclose(nexttrain); /* The training that still has to be done, containing all unambiguous pairs and all that remains of the ambiguous pairs. */
        }

    if(nflexrulesFormat)
        {
        assert(cutoff >= 0);
        char name[500];
        if(sizeof(name) <= (size_t)sprintf(name, nflexrulesFormat, 0))
            {
            fprintf(stderr, "doTraining: name 2 small");
            exit(-1);
            }
        countNodes(top, Counts + 0, options);
        options->info(6, "Calling writeRules %s %s\n", name, ext);
        writeRules(top, ext, 0, name, options);
        for(int thresh = 1; thresh <= cutoff; thresh++)
            {
#ifdef OLDPRUNE
            top->pruneAll(thresh);
#else
            assert(top->IfPatternFails == NULL);
            if(top->IfPatternSucceeds)
                top->IfPatternSucceeds = top->IfPatternSucceeds->Prune((unsigned long)thresh);
#endif
            top = top->cleanup(NULL, 0);
            if(sizeof(name) <= (size_t)sprintf(name, nflexrulesFormat, thresh))
                {
                fprintf(stderr, "doTraining: name 2 small");
                exit(-1);
                }
            countNodes(top, Counts + thresh, options);
            options->info(6, "Calling writeRules %s %s thresh %d\n", name, ext, thresh);
            writeRules(top, ext, thresh, name, options);
            }
        options->info(6, "Calling writeRules %s %s DONE\n", name, ext);
        }
    else
        {
        countNodes(top, Counts + 0, options);
        }
    delete top;
    building = false; // Signal to ~vertexPointer() to not access nodes.
    delete[] TrainingPair;
    building = true; // Signal to ~vertexPointer() to access nodes.
    return moreToDo;
    }

const unsigned int partOfFile(const char* fbuf, const double fraction, optionStruct * options)
    {
    if(options->currentParms() && !flog)
        {
        flog = fopenOrExit(options->currentParms(), "a", "log file");
        fprintf(flog, "%s: blobs=%d lines=%d fraction=%f most penalized=%d\n", options->wordLemmaList(), options->blobs(), options->lines(), fraction, node::mostPenalized);
        --openfiles;
        fclose(flog);
        flog = 0;
        }
    FILE* f2 = fopenOrExit(fbuf, "w", "computeParms");
    double bucket = fraction;
    int kar;
    FILE* f = fopen(options->wordLemmaList(), "rb");
    ++openfiles;
    bool tabSeen = false;
    if((double)options->blobs() * fraction > 1.0)
        {
        int bl = 1;
        int li = 0;
        int blbs = 0; // If there are no non-empty lines, there are no blobs
                     // either.
        int prevkar = 0;
        unsigned int lineProfile = 0;
        while((kar = fgetc(f)) != EOF)
            {
            if(kar == '\r')
                continue;
            if(kar == '\n')
                {
                if(lineProfile > 1)
                    { // Just finished non-empty line
                    lineProfile >>= 1; // push line before previous line out of
                                       // memory
                    // if lineProfile is 1, we may have a bunch of empty lines
                    // before a new non-empty line appears
                    }
                if(prevkar == '\n')
                    {
                    if(bucket >= 1.0)
                        {
                        ++bl;
                        bucket -= 1.0;
                        }
                    bucket += fraction;
                    }
                else if(prevkar != 0)
                    {
                    if(bucket >= 1.0)
                        {
                        ++li;
                        fputc(kar, f2);
                        }
                    }
                }
            else
                {
                if(lineProfile == 1)
                    {
                    // previous line was empty. Before that, there was a 
                    // non-empty line. Blob boundary detected!
                    ++blbs;
                    }
                lineProfile |= 4; // set third bit when making non-empty line
                if(bucket >= 1.0)
                    fputc(kar, f2);
                }
            prevkar = kar;
            }
        if(lineProfile & 4)
            { // last line was not finished off with newline
            ++li;
            }
        if(li > 0)
            ++blbs;
        options->info(0, "%d of %d lines, %d of %d clusters\n", li, options->lines(), bl, blbs);
        if(options->currentParms() && !flog)
            {
            flog = fopenOrExit(options->currentParms(), "a", "log file");
            ++openfiles;
            fprintf(flog, "Read %d blobs\n", bl);
            --openfiles;
            fclose(flog);
            flog = 0;
            }
        options->setReadBlobs(bl);
        }
    else if((double)options->lines() * fraction > 1.0)
        {
        while((kar = fgetc(f)) != EOF)
            {
            if(bucket >= 1.0)
                fputc(kar, f2);
            if(kar == '\n')
                {
                if(tabSeen)
                    {
                    if(bucket >= 1.0)
                        bucket -= 1.0;
                    bucket += fraction;
                    tabSeen = false;
                    }
                }
            else if(kar == '\t')
                tabSeen = true;
            }
        }
    else
        {
        while((kar = fgetc(f)) != EOF)
            {
            fputc(kar, f2);
            }
        }
    --openfiles;
    fclose(f);
    --openfiles;
    fclose(f2);

    f2 = fopenOrExit(fbuf, "r", "computeParms");
    unsigned int fraclines = 0;
    tabSeen = false;
    while((kar = fgetc(f2)) != EOF)
        {
        if(kar == '\n')
            {
            if(tabSeen)
                {
                tabSeen = false;
                ++fraclines;
                }
            }
        else if(kar == '\t')
            tabSeen = true;
        }
    assert(fraclines != 0);
    --openfiles;
    fclose(f2);
    options->setReadLines(fraclines);
    if(options->currentParms() && !flog)
        {
        flog = fopenOrExit(options->currentParms(), "a", "log file");
        fprintf(flog, "Use %d lines of %d\n", fraclines, options->lines());
        --openfiles;
        fclose(flog);
        flog = 0;
        }
    options->info(0, "reading %d lines\ntmpnam %s\n", fraclines, fbuf);
    return fraclines;
    }

void computeParms(optionStruct * options)
    {
    CHECK("iglobTempDir");
    int maxswath = options->swaths();//MAXSWATH;
    int currentNo = 0;
    int brownNo = 0;
    double currentweight = 0.0;
    double brownweight = 0.0;
    double fraction = 0.0; // 0.0 <= fraction <= 1.0
    double factor = 0.0;
    double iterationsfactor = 1;
    double miniterations = options->minIterations();//MINITERATIONS;
    double maxiterations = options->maxIterations();//MAXITERATIONS;
    if(maxiterations > options->lines())
        maxiterations = options->lines();
    node::mostPenalized = options->expectedCutoff() + 1; // parameter to weight function
    if(options->minfraction() > 0.0)
        {
        factor = pow(options->maxfraction() / options->minfraction(), 1.0 / (double)maxswath);
        // minfraction * factor^0 == minfraction
        // minfraction * factor^maxswath == maxfraction
        }

    if(miniterations > 0)
        {
        iterationsfactor = pow(maxiterations / miniterations, 1.0 / (double)maxswath);
        // maxiterations * iterationsfactor^0 == maxiterations
        // maxiterations * iterationsfactor^-maxswath == miniterations
        }
    else
        iterationsfactor = 1;
    const char* filename = options->wordLemmaList();
    const char* fbuf;
    fbuf = dup(tempFolder("trainFraction", options));

    char ext[100];
    ext[0] = '\0';
    if(sizeof(ext) <= (size_t)sprintf(ext, "%s", options->extra()))
        {
        fprintf(stderr, "computeParms: ext 2 small");
        exit(-1);
        }
    countAndWeight Count;
    int br1 = 0, br2 = 0;
    for(int swath = 0; swath <= maxswath; ++swath)
        {
        int blobs = 1;
        size_t lines = 0;
        size_t fraclines = 0;

        options->info(5, "Computing parameters: swath %d of %d\n", swath, maxswath);

        if(options->minfraction() > 0.0)
            {
            fraction = (swath == maxswath) ? options->maxfraction() : options->minfraction() * pow(factor, (double)swath);
            if(fraction == 1.0)
                filename = options->wordLemmaList();
            else
                {
                filename = fbuf;
                options->info(5, "Computing parameters: take a sample of training data\n");
                fraclines = partOfFile(fbuf, fraction, options);
                options->info(5, "Computing parameters: take a sample of training data DONE\n");
                }
            CHECK("D1globTempDir");
            brown(); // until not all parms are zero
            ++br1;
            if(swath == 0)
                init(options);
            else
                copybest(); // go on with best result so far.
            size_t filelines;

            options->info(5, "Computing parameters: initial training\n");

            aFile afile(filename, options, options->columns());
            filelines = afile.lines;

            options->info(6, "A Calling doTraining %s\n", ext);
            doTraining
            (/* aFile &                                     */  afile
             ,/* const char *                                */  ext
             ,/* int cutoff                                  */  0
             ,/* const char * nflexrulesFormat               */  0//options->flexrules()
             ,/* char * pairsToTrainInNextPassName           */  NULL
             ,/* countAndWeight *                            */  &Count
             , options
            ); // sets Nnodes
            options->info(6, "A Calling doTraining %s DONE\n", ext);
            if(lines == 0)
                fraclines = lines = filelines;

            brownNo = Count.getNnodes();
            currentNo = brownNo;
            switch(options->getWeightFunction())
                {
                case esupport:
                case eentropy:
                    brownweight = Count.getWeight();
                    break;
                case edepth:
                    brownweight = (double)Count.getCountByDepth();
                    break;
                case esize:
                    brownweight = (double)Count.getCountBySize();
                    break;
                case econstant:
                default:
                    brownweight = 0.0;
                }
            currentweight = brownweight;
            betterfound(currentNo, currentweight, swath, -1, blobs, lines, fraction, fraclines, false, options);
            printparms(Count.getNnodes(), brownweight, options);
            }
        int looplimit = (int)(maxiterations * pow(iterationsfactor, -swath));
        options->info(5, "Computing parameters: iterate and train, improve parameters by trial and error.\n");

        for(int iterations = 0; iterations < looplimit; ++iterations)
            {
            CHECK("D2aglobTempDir");
            brown();
            ++br2;
            CHECK("D2bglobTempDir");
            if(options->currentParms() && !flog)
                {
                flog = fopenOrExit(options->currentParms(), "a", "log file");
                CHECK("D2dglobTempDir");
                fprintf(flog, "//iteration:%d.%d %s\n", swath, iterations
                        , options->getWeightFunction() == esupport ? "weights"
                        : options->getWeightFunction() == eentropy ? "entropy"
                        : options->getWeightFunction() == edepth ? "depth"
                        : options->getWeightFunction() == esize ? "size"
                        : "count");
                CHECK("D2eglobTempDir");
                --openfiles;
                fclose(flog);
                flog = 0;
                }
            CHECK("D2fglobTempDir");
            options->info(5, "%d.%d ", swath, iterations);
            CHECK("D2globTempDir");

            aFile afile(filename, options, options->columns());
            size_t filelines = afile.lines;

            options->info(6, "B Calling doTraining %s\n", ext);
            doTraining
            (/* aFile &                                     */  afile
             ,/* const char *                                */  ext
             ,/* int cutoff                                  */  0
             ,/* const char * nflexrulesFormat               */  0//options->flexrules()
             ,/* char * pairsToTrainInNextPassName           */  NULL
             ,/* countAndWeight *                            */  &Count
             , options
            ); // sets Nnodes
            options->info(6, "B Calling doTraining %s DONE\n", ext);
            if(lines == 0)
                fraclines = lines = filelines;

            options->info(5, "%c%d %d %f %d %f           \n", STARTLINE, iterations, currentNo, currentweight, brownNo, brownweight);

            if(currentNo == 0)
                {
                currentNo = Count.getNnodes();
                switch(options->getWeightFunction())
                    {
                    case esupport:
                    case eentropy:
                        currentweight = Count.getWeight();
                        break;
                    case edepth:
                        currentweight = (double)Count.getCountByDepth();
                        break;
                    case esize:
                        currentweight = (double)Count.getCountBySize();
                        break;
                    case econstant:
                    default:
                        ;
                    }
                printparms(currentNo, currentweight, options);
                }
            else
                {
                brownNo = Count.getNnodes();
                switch(options->getWeightFunction())
                    {
                    case esupport:
                    case eentropy:
                        brownweight = Count.getWeight();
                        break;
                    case edepth:
                        brownweight = (double)Count.getCountByDepth();
                        break;
                    case esize:
                        brownweight = (double)Count.getCountBySize();
                        break;
                    case econstant:
                    default:
                        brownweight = 0.0;
                    }
                printparms(brownNo, brownweight, options);

                options->info(5, "swath %d brownNo %d currentNo %d\n", swath, brownNo, currentNo);

                if((brownNo <= currentNo
                    && options->getWeightFunction() == econstant
                    )
                   || (brownweight <= currentweight
                       && (options->getWeightFunction() == esupport
                           || options->getWeightFunction() == eentropy
                           || options->getWeightFunction() == edepth
                           || options->getWeightFunction() == esize
                           )
                       )
                   )
                    {
                    bool improvement = (((brownNo < currentNo)
                                         && options->getWeightFunction() == econstant
                                         )
                                        || ((brownweight < currentweight)
                                            && (options->getWeightFunction() == esupport
                                                || options->getWeightFunction() == eentropy
                                                || options->getWeightFunction() == edepth
                                                || options->getWeightFunction() == esize
                                                )
                                            )
                                        );
                    options->info(5, "%s\n", improvement ? "IMPROVEMENT" : "same");
                    currentNo = brownNo;
                    currentweight = brownweight;
                    betterfound(currentNo, currentweight, swath, iterations, blobs, lines, fraction, fraclines, improvement, options);
                    }
                else
                    {
                    worsefound();
                    }
                }
            }
        options->info(5, "br1 %d br2 %d\n", br1, br2);
        }
    remove(fbuf);
    delete[] fbuf;
    }

bool canwriteindir(const char* name)
    {
    size_t L = strlen(name);
    bool res = false;
    char* testfile = new char[L + 20];
    if(L > 0)
        sprintf(testfile, "%s%cT1E2S3T4F5I6L7E8", name, DIRSEP);
    else
        sprintf(testfile, "T1E2S3T4F5I6L7E8");
    FILE* fptest = fopen(testfile, "w");
    if(fptest)
        {
        fclose(fptest);
        remove(testfile);
        res = true;
        }
    delete[] testfile;
    return res;
    }

bool haswritabledir(const char* name)
    {
    size_t L = strlen(name);
    if(L > 0)
        {
        bool res = false;
        char* dirname = new char[L + 20];
        strcpy(dirname, name);
        if(dirname[L - 1] == '\\' || dirname[L - 1] == '/')
            dirname[L - 1] = 0;

        if(canwriteindir(dirname))
            res = true;
        else
            {
            char command[1024];
            sprintf(command, "%s %s", MKDIR, dirname);
            if(system(command)) return false;
            res = canwriteindir(dirname);
            }
        delete[] dirname;
        return res;
        }
    else
        return canwriteindir("");
    }

void trainRules(optionStruct * options, countAndWeight * Counts)
    {
    CHECK("jglobTempDir");
    options->info(2, "trainRules\n");
    assert(options->flexrules() != NULL);
    const char* tag = options->POStag();
    const char* nflexrules = options->flexrules();
    const char* fname = options->wordLemmaList();
    bool moreToDo = true;
    int passes = 0;
    char pairsToTrainInNextPassName[1024];
    char ingestedFractionOfAmbiguousPairsName[1024];
    char allPairsName[1024];
    char allIngestedPairsName[1024];
    char pairsToTrainInNextPassFormat[1024];
    char ingestedFractionOfAmbiguousPairsFormat[1024];
    char allPairsFormat[1024];
    char allIngestedPairsFormat[1024];
    char bestRulesFormat[1024];
    char numbersFormat[1024];
    //char command[1024];
    char FlexrulePassFormat[1024];
    char AccumulatedFlexrulePassFormat[1024];
    const char* displaytag = tag ? tag : "";
    sprintf(FlexrulePassFormat, "%s%c%%s.pass%%d.cutoff%%%%d", options->tempDir(), DIRSEP);
    sprintf(AccumulatedFlexrulePassFormat, "%s%c%%s.pass%%d.cutoff%%%%d.accumulated", options->tempDir(), DIRSEP);
    sprintf(pairsToTrainInNextPassFormat, "pairsToTrainInNextPass.%s%s.pass%%d", options->extra(), displaytag);
    sprintf(ingestedFractionOfAmbiguousPairsFormat, "ingestedFractionOfAmbiguousPairs.%s%s.pass%%d", options->extra(), displaytag);
    sprintf(allPairsFormat, "allPairs.%s%s.pass%%d", options->extra(), displaytag);
    sprintf(allIngestedPairsFormat, "allIngestedPairs.%s%s.pass%%d", options->extra(), displaytag);
    sprintf(numbersFormat, "numbers.%s%s.pass%%d", options->extra(), displaytag);
    char nflexrulesTag[1256];
    if(nflexrules)
        {
/*        if (tag && *tag)
            sprintf(nflexrulesTag, "%s.%s", nflexrules, tag);
        else*/
        strcpy(nflexrulesTag, nflexrules);
        nflexrules = nflexrulesTag;
        }

    const char* accumulatedFormat = FlexrulePassFormat;
    const char* accumulatedFormatPrev = NULL;

    do
        {
        ++passes;
        options->info(3, "do pass %d\n", passes);
        char flexrulesPass[1256];
        sprintf(flexrulesPass, FlexrulePassFormat, nflexrules, passes);
        char ext[100];
        ext[0] = '\0';
        if(options->extra())
            strcpy(ext, options->extra());

        sprintf(pairsToTrainInNextPassName, pairsToTrainInNextPassFormat, passes);
        sprintf(ingestedFractionOfAmbiguousPairsName, ingestedFractionOfAmbiguousPairsFormat, passes);
        sprintf(allPairsName, allPairsFormat, passes);
        sprintf(allIngestedPairsName, allIngestedPairsFormat, passes);
        char spass[20];
        if(sizeof(spass) <= (size_t)sprintf(spass, ".pass%d", passes))
            {
            fprintf(stderr, "trainRules: spass 2 small");
            exit(-1);
            }
        if(sizeof(ext) <= (size_t)sprintf(ext, "%s%s", options->extra(), spass))
            {
            fprintf(stderr, "trainRules: ext 2 small");
            exit(-1);
            }

        if(options->externalTrainer())
            {
            options->info(1, "options->externalTrainer() %s\n", options->externalTrainer());
            if(sizeof(bestRulesFormat) <= (size_t)sprintf(bestRulesFormat, accumulatedFormat, nflexrules, 1))
                {
                fprintf(stderr, "trainRules: bestRulesFormat 2 small");
                exit(-1);
                }
            for(int threshold = 0; threshold <= options->cutoff(); ++threshold)
                {
                options->info(2, "threshold %d\n", threshold);
                char* dest = new char[strlen(bestRulesFormat) + 10];
                sprintf(dest, bestRulesFormat, threshold);
                char* command = new char[strlen(options->externalTrainer()) + strlen(fname) + strlen(dest) + strlen("noofrules.txt") + 15];
                sprintf(command, "%s %s %s %s", options->externalTrainer(), fname, dest, "noofrules.txt");
                delete[] dest;
                options->info(2, "External trainer command: [%s]\n", command);
                if(!system(command))
                    {
                    FILE* noofrules = fopen("noofrules.txt", "r");
                    if(noofrules)
                        {
                        char buffer[100];
                        if(fread(buffer, 1, sizeof(buffer), noofrules) == 0)
                            {
                            fprintf(stderr, "Cannot read characters from a buffer open \"noofrules.txt\".\n");
                            Counts[threshold].setNnodes(1);
                            }
                        else
                            {
                            long rulecount = strtol(buffer, NULL, 10);
                            Counts[threshold].setNnodes(rulecount);
                            }
                        fclose(noofrules);
                        }
                    else
                        Counts[threshold].setNnodes(1);
                    }
                else
                    {
                    Counts[threshold].setNnodes(1);
                    }
                delete[] command;
                }
            moreToDo = false;
            }
        else
            {
            aFile afile(fname, options, passes > 1 ? "12" : options->columns());
            options->info(6, "C Calling doTraining %s, %s, %s\n", ext, flexrulesPass, pairsToTrainInNextPassName);
            moreToDo = doTraining
            (/* aFile &          */  afile
             ,/* const char *     */  ext
             ,/* int cutoff       */  options->cutoff()
             ,/* const char *     */  flexrulesPass
             ,/* char *           */  pairsToTrainInNextPassName
             ,/* countAndWeight * */  Counts
             ,/* optionStruct *   */  options
            );
            options->info(6, "C Calling doTraining %s, %s, %s DONE\n", ext, flexrulesPass, pairsToTrainInNextPassName);
            }
        /*
        Re-do the training, but only with those pairs that made it into
        the set of ingested pairs. The idea is to avoid the "noise" caused
        by the siblings of ambiguous pairs that didn't make it.
        The result of this excercise can be a better rule tree.
        */
        if(options->redo() && moreToDo)
            {
            options->info(5, "More training to do with file \"%s\"\n", pairsToTrainInNextPassName);

            aFile afile(tempFolder(allIngestedPairsName, options), options, options->columns());

            options->info(6, "D Calling doTraining %s, %s\n", ext, flexrulesPass);
            if(doTraining
            (/* aFile &                           */  afile
             ,/* const char *                      */  ext
             ,/* int                               */  options->cutoff()
             ,/* const char *                      */  flexrulesPass
             ,/* char * pairsToTrainInNextPassName */  NULL
             ,/* countAndWeight *                  */  Counts
             ,/* optionStruct *                    */  options
            )
               ) // sets Nnodes
                {
                options->info(0, "\"doTraining\" should return 0\n");
                }
            options->info(6, "D Calling doTraining %s, %s DONE\n", ext, flexrulesPass);
            }
        else if(moreToDo)
            {
            options->info(1, "No retraining done on ingested pairs, although ambiguous pairs were found and may have caused noise. (Faster) \n");
            }

        char AccumulatedFlexrulesPass[256];
        char NextAccumulatedFlexrulesPassFormat[256];
        if(sizeof(bestRulesFormat) <= (size_t)sprintf(bestRulesFormat, accumulatedFormat, nflexrules, passes))
            {
            fprintf(stderr, "trainRules: bestRulesFormat 2 small");
            exit(-1);
            }
        if(passes > 1)
            {
            assert(accumulatedFormatPrev != NULL);
            if(sizeof(AccumulatedFlexrulesPass) <= (size_t)sprintf(AccumulatedFlexrulesPass, accumulatedFormatPrev, nflexrules, passes - 1))
                {
                fprintf(stderr, "trainRules: AccumulatedFlexrulesPass 2 small");
                exit(-1);
                }
            }
        if(sizeof(NextAccumulatedFlexrulesPassFormat) <= (size_t)sprintf(NextAccumulatedFlexrulesPassFormat, AccumulatedFlexrulePassFormat, nflexrules, passes))
            {
            fprintf(stderr, "trainRules: NextAccumulatedFlexrulesPassFormat 2 small");
            exit(-1);
            }

        for(int cut = 0; cut <= options->cutoff(); ++cut)
            {
            char nextbestflexrules[1150];
            if(sizeof(nextbestflexrules) <= (size_t)sprintf(nextbestflexrules, flexrulesPass, cut))
                {
                fprintf(stderr, "trainRules: nextbestflexrules 2 small");
                exit(-1);
                }
            //prettyPrint(nextbestflexrules);
            if(passes > 1)
                {
                char bestflexrules[1150], newbestflexrules[1150];
                if(sizeof(bestflexrules) <= (size_t)sprintf(bestflexrules, AccumulatedFlexrulesPass, cut))
                    {
                    fprintf(stderr, "trainRules: bestflexrules 2 small");
                    exit(-1);
                    }
                if(sizeof(newbestflexrules) <= (size_t)sprintf(newbestflexrules, NextAccumulatedFlexrulesPassFormat, cut))
                    {
                    fprintf(stderr, "trainRules: newbestflexrules 2 small");
                    exit(-1);
                    }
                options->info(5, "flexcombi best %s + next best %s -> combined %s\n", bestflexrules, nextbestflexrules, bestflexrules);
                if(!flexcombi(bestflexrules, nextbestflexrules, newbestflexrules))
                    break;
                //prettyPrint(newbestflexrules);
                }
            else
                { /* Purpose: use flexcombi to insert version in binary flexrule file.
                     If the training set is unambiguous, only one pass is needed,
                     and the version number would be missing, causing cstlemma to
                     protest that the flexrules are of obsolete design.
                   */
                char newbestflexrules[1150];
                if(sizeof(newbestflexrules) <= (size_t)sprintf(newbestflexrules, NextAccumulatedFlexrulesPassFormat, cut))
                    {
                    fprintf(stderr, "trainRules: newbestflexrules 2 small");
                    exit(-1);
                    }
                options->info(5, "flexcombi best %s + next best %s -> combined %s\n", nextbestflexrules, nextbestflexrules, nextbestflexrules);
                if(!flexcombi(nextbestflexrules, nextbestflexrules, newbestflexrules))
                    break;
                remove(nextbestflexrules);
                rename(newbestflexrules, nextbestflexrules);
                }
            }
        fname = tempFolder(pairsToTrainInNextPassName, options);

        accumulatedFormatPrev = accumulatedFormat;
        accumulatedFormat = AccumulatedFlexrulePassFormat;
        }
    while(moreToDo && passes < options->maxPasses());

    options->info(5, "\n"); // to compensate for missing newline in readTrainingPairs
    if(options->cutoff() >= 0)
        {
        char dest[1000];
        dest[0] = '\0';
        for(int cut = 0; cut <= options->cutoff(); ++cut)
            {
            char dirname[500];
            const char* lastslash = strrchr(nflexrules, DIRSEP);
            const char* filename;
            if(lastslash)
                {
                filename = lastslash + 1;
                sprintf(dirname, "%.*s%c%d", (int)(lastslash - nflexrules), nflexrules, DIRSEP, cut);
                }
            else
                {
                filename = nflexrules;
                sprintf(dirname, "%d", cut);
                }
            options->info(3, "dirname %s\n", dirname);

            char testfile[512];
            sprintf(testfile, "%s%cTESTFILE", dirname, DIRSEP);
            options->info(3, "testfile %s\n", testfile);
            FILE* fptest = fopen(testfile, "w");
            ++openfiles;
            bool hasDir = false;
            if(fptest)
                {
                options->info(3, "testfile created\n");
                --openfiles;
                fclose(fptest);
                remove(testfile);
                options->info(3, "testfile deleted\n");
                hasDir = true;
                }
            else
                {
                options->info(3, "testfile not created\n");
                char command[1024];
                sprintf(command, "%s %s", MKDIR, dirname);
                options->info(3, "command %s\n", command);
                if(!system(command))
                    {
                    fptest = fopen(testfile, "w");
                    ++openfiles;
                    if(fptest)
                        {
                        options->info(3, "testfile created\n");
                        --openfiles;
                        fclose(fptest);
                        remove(testfile);
                        options->info(3, "testfile deleted\n");
                        hasDir = true;
                        }
                    }
                }
            if(hasDir)
                {
                char finalPath[2048];
                sprintf(finalPath, "%s%c%s", dirname, DIRSEP, filename);
                remove(finalPath);
                sprintf(dest, bestRulesFormat, cut);
                rename(dest, finalPath);
                char* lastdest = dest + strlen(dest);
                char* lastfinalPath = finalPath + strlen(finalPath);
                if(prettytxt)
                    {
                    strcpy(lastdest, prettytxt);
                    strcpy(lastfinalPath, prettytxt);
                    rename(dest, finalPath);
                    }
                if(prettybra)
                    {
                    strcpy(lastdest, prettybra);
                    strcpy(lastfinalPath, prettybra);
                    rename(dest, finalPath);
                    }
                }
            }
        }
    for(int cut = 0; cut <= options->cutoff(); ++cut)
        {
        char scut[100];
        sprintf(scut, "%s%d", SCUT, cut);
        for(int passs = 1;; ++passs)
            {
            char dest[1500];
            char spass[100];
            sprintf(spass, ".pass%d", passs);
            sprintf(dest, "%s%s%s", nflexrules, scut, spass);
            options->info(3, "Remove %s \n", dest);
            if(remove(dest))
                break;
            }
        }
    }


FILE* flog = NULL;

static void initOutput(const char* path)
    {
    if(!path)
        return;
    FILE* fp = fopenOrExit(path, "w", "initOutput");
    --openfiles;
    fclose(fp);
    }

void ParameterComputation(optionStruct * poptions/*, const char * tagName*/)
    {
    if(poptions->computeParms())
        {
        poptions->info(2, "Computing parameters\n");

        if(poptions->currentParms())
            initOutput(poptions->currentParms());

        if(poptions->bestParms())
            initOutput(poptions->bestParms());

        computeParms(poptions);

        poptions->info(2, "Computing parameters DONE\n");
        }
    }

void ComputeDelta(optionStruct * poptions)
    {
    poptions->info(2, "Going to compute delta.\n");

    init(poptions); // TODO Check that penalties are the best ones, not the last ones tried.

    poptions->info(2, "Computing delta DONE.\n");
    }

void TestRules(optionStruct * poptions)
    {
    poptions->info(1, "TestRules\n");
    if(poptions->test() || (poptions->trainTest() && !poptions->tenfoldCrossValidation()))
        {
        poptions->info(2, "Going to test the rules.\n");

        poptions->printArgFile();
        testrules(poptions);
        poptions->info(2, "Testing the rules DONE.\n");
        }
    }

void CreateFlexRules(optionStruct * poptions)
    {
    if(poptions->createFlexRules())
        {
        setCompetitionFunction(poptions);
        poptions->setReadLines(poptions->lines());
        poptions->info(2, "createFlexRules()\n");

        countAndWeight* Counts = new countAndWeight[1 + (size_t)poptions->cutoff()];
        trainRules(poptions, Counts);
        delete[]Counts;
        }
    }

static bool printNiceInput(optionStruct * options, const char* tagName)
    {
    aFile afile(options->wordLemmaList(), options, options->columns());
    size_t allPairs;
    trainingPair* TrainingPair = readTrainingPairs(afile, allPairs, options);
    if(TrainingPair)
        {
        markTheAmbiguousPairs(TrainingPair, allPairs, options);
        trainingPair* train = createChainOfPointersToTrainingPairs(TrainingPair, allPairs);
        saveChainOfPointersToTrainingPairsToFile(train, "clearedFile", tagName);
        return true;
        }
    return false;
    }

static bool AllProcessing(optionStruct * poptions, const char* tagName)
    {
    poptions->info(0, "AllProcessing %s\n", tagName);

    if(printNiceInput(poptions, tagName))
        {
        ParameterComputation(poptions/*,tagName*/);
        ComputeDelta(poptions);
        TestRules(poptions);
        CreateFlexRules(poptions);
        poptions->info(0, "AllProcessing %s DONE\n", tagName);
        return true;
        }
    return false;
    }


int main(int argc, char** argv)
    {
    if(argc < 2)
        {
        printf("affixtrain - supervised learning of affix rules for CSTLEMMA, version " VERSION "\n");
        printf("%s -h for usage\n", argv[0]);
        return 0;
        }
    optionStruct options;
    switch(options.readArgs(argc, argv))
        {
        case GoOn:
            break;
        case Leave:
            return 0;
        case Error:
            fprintf(stderr, "Error: Error in options. Exiting\n");
            return 1;
        }

    if(options.verbose())
        {
        options.print(stderr);
        options.printArgFile();
        }


    if(options.rawRules())
        {
        prettyPrint(options.rawRules());
        prettyPrintBracmat(options.rawRules());
        if(options.wordPerLineList())
            {
            ::lemmatiseFile(options.wordPerLineList(), options.rawRules(), options.lemmas());
            }
        }
    else
        {
        if(options.tempDir())
            {
            if(!haswritabledir(options.tempDir()))
                {
                fprintf(stderr, "Cannot create file %s. Did you specify an existing and writable temp directory? (option -j)\n", tempFolder("testFile", &options));
                return -1;
                }
            }
        else if(!canwriteindir("."))
            {
            fprintf(stderr, "Cannot write in the working directory.\n");
            return -1;
            }

        tagClass* Tags = NULL;
        Tags = collectTags(&options);
        tagClass* theTag = Tags;
        if(theTag)
            {
            options.info(0, "Doing Tags\n");
            optionStruct originalOptions(options);
            while(theTag)
                {
                optionStruct tagOptions(options);
                tagOptions.info(0, "Doing tag %s\n", theTag->name);
                FilterTagLines(originalOptions.wordLemmaList(), &tagOptions, theTag->name);
                if(!AllProcessing(&tagOptions, theTag->name))
                    break;
                theTag = theTag->next;
                }
            }
        else
            {
            options.info(0, "NOT doing Tags\n");
            if(!AllProcessing(&options, ""))
                {
                FILE* fp;
                char name[256];
                sprintf(name, "empty-flexrules-%s", options.argstring());
                fp = fopen(name, "wb");
                fwrite(emptyRules,1,12,fp);
                fclose(fp);
                }
            }


        options.info(0, "\nAffixTrain OK\n");
        }
    cleanUpOptions();
    return 0;
    }
