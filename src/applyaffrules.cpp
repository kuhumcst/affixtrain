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
#if DOTEST
#include <stdio.h>
#include <string.h>
#include <assert.h>

struct var
    {
    const char * s;
    const char * e;
    };

static char * buf = NULL;
static long end;
static char * result = NULL;
static bool NewStyle = false;
static bool newStyleRules()
    {
    return NewStyle;
    }

static const char * samestart(const char ** fields,const char * s,const char * we)
    {
    const char * f = fields[0];
    const char * e = fields[1]-1;
    while((f < e) && (s < we) && (*f == *s))
        {
        ++f;
        ++s;
        }
    return f == e ? s : NULL;
    }

static const char * sameend(const char ** fields,const char * s,const char * wordend)
    {
    const char * f = fields[2];
    const char * e = fields[3]-1;
    const char * S = wordend - (e-f);
    if(S >= s)
        {
        s = S;
        while(f < e && *f == *S)
            {
            ++f;
            ++S;
            }
        }
    return f == e ? s : NULL;
    }

static bool substr(const char ** fields,int k,const char * w,const char * wend,var * vars,int vindex)
    {
    if(w == wend)
        return false;
    const char * f = fields[k];
    const char * e = fields[k+1]-1;
    const char * p = w;
    assert(f != e);
    const char * ff;
    const char * pp;
    do
        {
        while((p < wend) && (*p != *f))
            {
            ++p;
            }
        if(p == wend)
            return false;
        pp = ++p;
        ff = f+1;
        while(ff < e)
            {
            if(pp == wend)
                return false;
            if(*pp != *ff)
                break;
            ++pp;
            ++ff;
            }
        }
    while(ff != e);
    vars[vindex].e = p-1;
    vars[vindex+1].s = pp;
    return true;
    }

/*
static void printpat(const char ** fields,int findex)
    {
    printf("findex %d {%.*s]",findex,fields[1] - fields[0] - 1,fields[0]);//,vars[0].e - vars[0].s,vars[0].s);
    for(int m = 1;2*m+3 < findex;++m)
        {
        int M = 2*m+3;
        printf("[%.*s]",fields[M] - fields[M-1] - 1,fields[M-1]);//,vars[m].e - vars[m].s,vars[m].s);
        }
    if(findex > 2)
        printf("[%.*s}\n",fields[3] - fields[2] - 1,fields[2]);
    else
        printf("\n");
    }
*/
static bool lemmatiseer(const char * word,const char * wordend,const char * buf,int maxpos)
    {
    int pos = 0;
    do
        {
        var vars[20];
        const char * fields[44];
        // output=fields[1]+vars[0]+fields[5]+vars[1]+fields[7]+vars[2]+...+fields[2*n+3]+vars[n]+...+fields[3]
        const char * wend = wordend;
        buf += pos;
        pos = *(int*)buf;
        const char * p = buf + sizeof(int);
        fields[0] = p;
        int findex = 1;
        while(*p != '\n')
            {
            if(*p == '\t')
                fields[findex++] = ++p;
            else
                ++p;
            }
        fields[findex] = ++p; // p is now within 3 bytes from the next record.
//        printpat(fields,findex);
        // check Lpat
        vars[0].s = samestart(fields,word,wend);
        if(vars[0].s)
            {
            // Lpat succeeded
            vars[0].e = wend;
            if(findex > 2)
                {
                const char * newend = sameend(fields,vars[0].s,wend);
                if(newend)
                    wend = newend;
                else
                    continue;

                int k;
                const char * w = vars[0].s;
                int vindex = 0;
                for(k = 4;k < findex;k += 2)
                    {
                    if(!substr(fields,k,w,wend,vars,vindex))
                        break;
                    ++vindex;
                    w = vars[vindex].s;
                    }
                if(k < findex)
                    continue;
                vars[vindex].e = newend;
                long nxt = p - buf;
                nxt += 3;
                nxt >>= 2;
                nxt <<= 2;
                if(pos)
                    maxpos = pos;
                if(nxt >= maxpos || !lemmatiseer(vars[0].s,wend,buf+nxt,maxpos))
                    {
                    int resultlength = (fields[2] - fields[1] - 1) + (vars[0].e - vars[0].s) + (fields[4] - fields[3] - 1);
                    int m;
                    for(m = 1;2*m+3 < findex;++m)
                        {
                        int M = 2*m+3;
                        resultlength += (fields[M+1] - fields[M] - 1) + (vars[m].e - vars[m].s);
                        }
                    result = new char[resultlength+1];
                    int printed = sprintf(result,"%.*s%.*s",fields[2] - fields[1] - 1,fields[1],vars[0].e - vars[0].s,vars[0].s);

                    for(m = 1;2*m+3 < findex;++m)
                        {
                        int M = 2*m+3;
                        printed += sprintf(result+printed,"%.*s%.*s",fields[M+1] - fields[M] - 1,fields[M],vars[m].e - vars[m].s,vars[m].s);
                        }
                    sprintf(result+printed,"%.*s",fields[4] - fields[3] - 1,fields[3]);
                    }
                }
            else if(vars[0].e == vars[0].s)
                { 
                result = new char[(fields[2] - fields[1] - 1)+1];
                sprintf(result,"%.*s",fields[2] - fields[1] - 1,fields[1]);
                }
            else
                continue;

    
            return true;
            }
        else
            {
            // Lpat failed
            continue;
            }
        }
    while(pos);
    return false;
    }


static bool readRules(FILE * flexrulefile)
    {
    fseek(flexrulefile,0,SEEK_END);
    end = ftell(flexrulefile);

    buf = new char[end+1];// 20140224 +1
    rewind(flexrulefile);
    fread(buf,1,end,flexrulefile);
    buf[end] = '\0';// 20140224 new
    fclose(flexrulefile);
    NewStyle = true;
    return true;
    }


bool readRules(const char * filename)
    {
    FILE * f = fopen(filename,"rb");
    if(f)
        {
        return readRules(f);
        }
    return false;
    }



void deleteRules()
    {
    delete [] buf;
    buf = NULL;
    delete [] result;
    result = NULL;
    NewStyle = false;
    }


const char * applyRules(const char * word)
    {
    if(buf)
        {
        int len = strlen(word);
        delete [] result;
        result = NULL;
        lemmatiseer(word,word+len,buf,end);
        return result;
        }
    return NULL;
    }

#endif
