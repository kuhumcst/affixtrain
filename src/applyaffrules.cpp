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
#include "affixtrain.h"
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

typedef unsigned char typetype;

struct var
    {
    const char * s;
    const char * e;
    };

class buffer
    {
    public:
        char * buf;
        long buflen;
        buffer()
            {
            buf = 0;
            buflen = 0;
            }
        ~buffer()
            {
            delete [] buf;
            }
    };

static char * result = 0;

static bool readRules(FILE * flexrulefile, buffer * Buffer)
    {
    if (flexrulefile)
        {
        int start;
        if (fread(&start, sizeof(int), 1, flexrulefile) != 1)
            return 0;
        fseek(flexrulefile, 0, SEEK_END);
        Buffer->buflen = ftell(flexrulefile);
        if (start == *(int*)"\rV3\r")
            {
            fseek(flexrulefile, sizeof(int), SEEK_SET);
            Buffer->buflen -= sizeof(int);
            }
        else if (start == 0)
            {
            fseek(flexrulefile, 0, SEEK_SET);
            }
        else
            return false;
        Buffer->buf = new char[Buffer->buflen + 1];
        if (Buffer->buf && Buffer->buflen > 0)
            {
            if (fread(Buffer->buf, 1, Buffer->buflen, flexrulefile) != (size_t)Buffer->buflen)
                return 0;
            Buffer->buf[Buffer->buflen] = '\0';
            }
        return true;
        }
    return false;
    }

static bool substr(const char ** fields, int k, const char * w, const char * wend, var * vars, int vindex)
    {
    if (w == wend)
        return false;
    const char * f = fields[k];
    const char * e = fields[k + 1] - 1;
    const char * p = w;
    assert(f != e);
    const char * ff;
    const char * pp;
    do
        {
        while ((p < wend) && (*p != *f))
            {
            ++p;
            }
        if (p == wend)
            return false;
        pp = ++p;
        ff = f + 1;
        while (ff < e)
            {
            if (pp == wend)
                return false;
            if (*pp != *ff)
                break;
            ++pp;
            ++ff;
            }
        } while (ff != e);
        vars[vindex].e = p - 1;
        vars[vindex + 1].s = pp;
        return true;
    }

static const char * samestart(const char ** fields, const char * s, const char * we)
    {
    const char * f = fields[0];
    const char * e = fields[1] - 1;
    while ((f < e) && (s < we) && (*f == *s))
        {
        ++f;
        ++s;
        }
    // On success: return pointer to first unparsed character
    // On failure: return 0
    return f == e ? s : 0;
    }

static const char * sameend(const char ** fields, const char * s, const char * wordend)
    {
    const char * f = fields[2];
    const char * e = fields[3] - 1;
    const char * S = wordend - (e - f);
    if (S >= s)
        {
        s = S;
        while (f < e && *f == *S)
            {
            ++f;
            ++S;
            }
        }
    // On success: return pointer to successor of last unparsed character
    // On failure: return 0
    return f == e ? s : 0;
    }

static char * rewrite(const char *& word, const char *& wordend, const char * p)
    {
    var vars[20];
    const char * fields[44]; // 44 = (2*20 + 3) + 1
    // The even numbered fields contain patterns
    // The odd numbered fields contain replacements
    // The first two fields (0,1) refer to the prefix
    // The third and fourth (2,3) fields refer to the affix
    // The remaining fields (4,5,..,..) refer to infixes, from left to right
    // output=fields[1]+vars[0]+fields[5]+vars[1]+fields[7]+vars[2]+...+fields[2*n+3]+vars[n]+...+fields[3]
    const char * wend = wordend;
    fields[0] = p;
    int findex = 1;
    while (*p != '\n')
        {
        if (*p == '\t')
            fields[findex++] = ++p;
        else
            ++p;
        }
    fields[findex] = ++p;
    //printpat(fields, findex);
    // fields[findex] points to character after \n. 
    // When 1 is subtracted, it points to the character following the last replacement.
    // p is now within 3 bytes from the first Record of the subtree
    //        printpat(fields,findex);
    // check Lpat
    vars[0].s = samestart(fields, word, wend);
    if (vars[0].s)
        {
        // Lpat succeeded
        vars[0].e = wend;
        char * destination = NULL;
        int printed = 0;
        if (findex > 2) // there is more than just a prefix
            {
            const char * newend = sameend(fields, vars[0].s, wend);
            if (newend)
                wend = newend;
            else
                return 0; //suffix didn't match

            int k;
            const char * w = vars[0].s;
            int vindex = 0;
            for (k = 4; k < findex; k += 2)
                {
                if (!substr(fields, k, w, wend, vars, vindex))
                    break;
                ++vindex;
                w = vars[vindex].s;
                }
            if (k < findex)
                return 0;

            vars[vindex].e = newend;
            //                     length of prefix       length of first unmatched         length of suffix
            ptrdiff_t resultlength = (fields[2] - fields[1] - 1) + (vars[0].e - vars[0].s) + (fields[4] - fields[3] - 1);/*20120709 int -> ptrdiff_t*/
            int m;
            for (m = 1; 2 * m + 3 < findex; ++m)
                {
                int M = 2 * m + 3;
                //                    length of infix       length of unmatched after infix
                resultlength += (fields[M + 1] - fields[M] - 1) + (vars[m].e - vars[m].s);
                }
            destination = new char[resultlength + 1];
            printed = sprintf(destination, "%.*s%.*s", (int)(fields[2] - fields[1] - 1), fields[1], (int)(vars[0].e - vars[0].s), vars[0].s);
            for (m = 1; 2 * m + 3 < findex; ++m)
                {
                int M = 2 * m + 3;
                printed += sprintf(destination + printed, "%.*s%.*s", (int)(fields[M + 1] - fields[M] - 1), fields[M], (int)(vars[m].e - vars[m].s), vars[m].s);
                }
            printed += sprintf(destination + printed, "%.*s", (int)(fields[4] - fields[3] - 1), fields[3]);
            word = vars[0].s;
            wordend = newend;
            }
        else if (vars[0].e == vars[0].s) // whole-word match: everything matched by "prefix"
            {
            //++news;
            destination = new char[(fields[2] - fields[1] - 1) + 1];
            printed = sprintf(destination, "%.*s", (int)(fields[2] - fields[1] - 1), fields[1]);
            }
        else
            return 0; // something unmatched

        return destination;
        }
    else
        {
        // Lpat failed
        return 0; // prefix failed
        }
    }

static char ** addLemma(char ** lemmas, const char * lemma)
    {
    if (lemma)
        {
        if (lemmas)
            {
            int i;
            for (i = 0; lemmas[i]; ++i)
                {
                if (!strcmp(lemmas[i], lemma))
                    {
                    return lemmas;
                    }
                }
            char ** nlemmas = new char *[i + 2];
            for (i = 0; lemmas[i]; ++i)
                {
                nlemmas[i] = lemmas[i];
                }
            delete[] lemmas;
            lemmas = nlemmas;
            lemmas[i] = new char[strlen(lemma) + 1];
            strcpy(lemmas[i], lemma);
            lemmas[++i] = 0;
            }
        else
            {
            lemmas = new char *[2];
            lemmas[1] = 0;
            lemmas[0] = new char[strlen(lemma) + 1];
            strcpy(lemmas[0], lemma);
            }
        }
    return lemmas;
    }


static char ** lemmatiseerV3(const char * word, const char * wordend, const char * buf, const char * maxpos, const char * parentcandidate, char ** lemmas);

static char ** chainV3(const char * word, const char * wordend, const char * buf, const char * maxpos, const char * parentcandidate, char ** lemmas)
    {
    for (int next = *(int*)buf
         ;
         ; buf += next, next = *(int*)buf
         )
        {
        assert((next & 3) == 0);
        assert(next == -4 || next == 0 || next == 4 || next >= 12);
        if (next == -4 || next == 4)
            {
            // add parent candidate to lemmas.
            lemmas = addLemma(lemmas, parentcandidate);
            }
        else
            {
            char ** temp = lemmatiseerV3(word, wordend, buf + sizeof(int), next > 0 ? buf + next : maxpos, parentcandidate, lemmas);
            if (temp)
                {
                lemmas = temp;
                }
            else
                {
                lemmas = addLemma(lemmas, parentcandidate);
                }
            }
        if (next <= 0)
            break;
        }
    return lemmas;
    }

static char ** lemmatiseerV3(const char * word, const char * wordend, const char * buf, const char * maxpos, const char * parentcandidate, char ** lemmas)
    {
    if (maxpos <= buf)
        return 0;
    const char * cword = word;
    const char * cwordend = wordend;
    int pos = 0;
    pos = *(int*)buf;
    const char * until;
    char ** result;
    assert((pos & 3) == 0);
    assert(pos >= 0);
    if (pos == 0)
        until = maxpos;
    else
        until = buf + pos;
    typetype type;
    const char * p = buf + sizeof(int);
    type = *(typetype*)p;
    /*
    buf+4:   
    first bit  0: Fail branch is unambiguous, buf points to tree. (A)
    first bit  1: Fail branch is ambiguous, buf points to chain. (B)
    second bit 0: Success branch is unambiguous, buf+8 points to tree (C)
    second bit 1: Success branch is ambiguous, buf+8 points to chain (D)
    */
    if (type < 4)
        {
        ++p;
        }
    else
        {
        type = 0; // no ambiguity
        }
    char * candidate = rewrite(cword, cwordend, p);
    p = strchr(p, '\n');
    ptrdiff_t off = p - buf;
    off += sizeof(int);
    off /= sizeof(int);
    off *= sizeof(int);
    p = buf + off;
    if (candidate)
        {
        const char * defaultCandidate = candidate[0] ? candidate : parentcandidate;
        /* 20150806 A match resulting in a zero-length candidate is valid for
        descending, but if all descendants fail, the candidate is overruled by
        an ancestor that is not zero-length. (The top rule just copies the
        input, so there is a always a non-zero length ancestor.) */
        switch (type)
            {
            case 0:
            case 1:
                {
                /* Unambiguous children. If no child succeeds, take the
                candidate, otherwise take the succeeding child's result. */
                char ** childcandidates = lemmatiseerV3(cword, cwordend, p, until, defaultCandidate, lemmas);
                result = childcandidates ? childcandidates : addLemma(lemmas, defaultCandidate);
                delete[] candidate;
                break;
                }
            case 2:
            case 3:
                {
                /* Ambiguous children. If no child succeeds, take the
                candidate, otherwise take the succeeding children's result
                Some child may in fact refer to its parent, which is our
                current candidate. We pass the candidate so it can be put
                in the right position in the sequence of answers. */
                char ** childcandidates = chainV3(cword, cwordend, p, until, defaultCandidate, lemmas);
                result = childcandidates ? childcandidates : addLemma(lemmas, defaultCandidate);
                delete[] candidate;
                break;
                }
            default:
                result = lemmas;
            }
        }
    else
        {
        switch (type)
            {
            case 0:
            case 2:
                {
                /* Unambiguous siblings. If no sibling succeeds, take the
                parent's candidate. */
                char ** childcandidates = lemmatiseerV3(word, wordend, until, maxpos, parentcandidate, lemmas);
                result = childcandidates ? childcandidates : addLemma(lemmas, parentcandidate);
                break;
                }
            case 1:
            case 3:
                {
                /* Ambiguous siblings. If a sibling fails, the parent's
                candidate is taken. */
                char ** childcandidates = chainV3(word, wordend, until, maxpos, parentcandidate, lemmas);
                result = childcandidates ? childcandidates : addLemma(lemmas, parentcandidate);
                break;
                }
            default:
                result = lemmas;
            }
        }
    return result;
    }

static char * concat(char ** L)
    {
    if (L)
        {
        int lngth = 0;
        int i;
        for (i = 0; L[i]; ++i)
            {
            lngth += strlen(L[i]) + 1;
            }
        char * ret = new char[lngth];
        ret[0] = 0;
        if(L[0])
            {
            strcat(ret, L[0]);
            delete[] L[0];
            }
        for (i = 1; L[i]; ++i)
            {
            strcat(ret, "|");
            strcat(ret, L[i]);
            delete[] L[i];
            }
        delete[] L;
        ret[lngth-1] = 0;
        return ret;
        }
    else
        return 0;
    }

bool readRules(const char * filename,buffer * Buffer)
    {
    FILE * f = fopen(filename,"rb");
    ++openfiles;
    if(f)
        {
        bool result = readRules(f,Buffer);
        --openfiles;
        fclose(f);
        return result;
        }
    return false;
    }

void deleteRules()
    {
    delete [] result;
    result = NULL;
    }


const char * applyRules(const char * word,buffer * Buffer)
    {
    if(Buffer->buf)
        {
        int len = strlen(word);
        delete [] result;
        result = NULL;
        char ** lemmas = 0;
        result = concat(lemmatiseerV3(word,word+len,Buffer->buf,Buffer->buf+Buffer->buflen,0,lemmas));
        return result;
        }
    return NULL;
    }

bool lemmatiseFile(const char * OneWordPerLineFile,const char * rulefile,const char * resultFile)
    {
    buffer Buffer;
    if(!readRules(rulefile,&Buffer))
        return false;

    FILE * OWPLF = fopen(OneWordPerLineFile,"rb");
    ++openfiles;
    if(!OWPLF)
        {
        deleteRules();
        return false;
        }

    FILE * rf = fopen(resultFile,"wb");
    ++openfiles;
    if(!rf)
        {
        --openfiles;
        fclose(OWPLF);
        deleteRules();
        return false;
        }

    fseek(OWPLF, 0, SEEK_END);
    long wend = ftell(OWPLF);
    fseek(OWPLF, 0, SEEK_SET);
    char * wbuf = new char[wend + 1];
    if (wbuf && wend > 0)
        {
        if (fread(wbuf, 1, wend, OWPLF) != (size_t)wend)
            return 0;
        wbuf[wend] = '\0';
        }

    char * q;
    for(char * p = wbuf;(q = strchr(p,'\n')) != 0;p = q + 1)
        {
        *q = 0;
        for(char * r = q - 1;r >= p && (*r == '\r' || *r == ' ' || *r == '\t');--r)
            *r = 0;
        fprintf(rf,"%s\n",applyRules(p,&Buffer));
        }

    delete [] wbuf;
    deleteRules();
    --openfiles;
    fclose(rf);
    --openfiles;
    fclose(OWPLF);
    return true;
    }

