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
#ifndef STRNG_H
#define STRNG_H

#include <stddef.h> // gcc ptrdiff_t
#include <stdarg.h>
#include <string.h>
#include <assert.h>

typedef enum {dif_smaller,dif_equal,dif_bigger,dif_incommensurable,dif_incompatible
/*
        dif_equalsuccessor, // rule is identical with a successor
        dif_addedsuccessor, 
        dif_identical       // used if both pattern and replacement are the same. In that case we have the same object in both hands.
        */
        } edif;  // rule has been added as a successor


char * dup(const char * buf);
char * dupn(const char * buf,size_t n);

extern int StrngCount;

class strng
    {
    private:
        char * Txt;
    public:
        char * itsTxt() const
            {
            return Txt;
            }
        const char * itsCTxt() const
            {
            return Txt;
            }
        void dup(const char * buf){Txt = ::dup(buf);}
        strng(const char * buf);
        strng(const char * buf,size_t length);
        strng(const strng * s)
            {
            Txt = ::dup(s->Txt);
            ++StrngCount;
            }
        strng(int kar)
            {
            Txt = new char[2];
            Txt[0] = (char)kar;
            Txt[1] = '\0';
            ++StrngCount;
            }
        strng():Txt(0)
            {
            ++StrngCount;
            }
        ~strng()
            {
            delete [] Txt;
            --StrngCount;
            }
        bool eq(const char * s)// returns Txt - s
            {
            return !strcmp(Txt, s);
            }
        bool eq(strng & s)// returns Txt - s
            {
            return !strcmp(Txt, s.Txt);
            }

        edif dif(strng * s); // returns this->Txt - s (1, 0 or -1)
        void checkIntegrity();
        bool hasWildCard();
        void cat(const strng * a,...)
            {
            va_list ap;
            va_start(ap,a);
            const strng * i = a;
            size_t len = strlen(this->Txt);
            while(i)
                {
                len += strlen(i->Txt);
                i = va_arg(ap,const strng *);
                }
            va_end(ap);
            char * buf = new char[len+1];
            char * p = buf;
            const char * j = this->Txt;
            while(*j)
                {
                *p++ = *j++;
                }
            va_list ap2;
            va_start(ap2,a);
            i = a;
            while(i)
                {
                j = i->Txt;
                while(*j)
                    {
                    *p++ = *j++;
                    }
                i = va_arg(ap2,const strng *);
                }
            va_end(ap2);
            *p = '\0';
            delete this->Txt;
            this->Txt = buf;
            }
        ptrdiff_t pos(int kar) const
            {
            const char * p = strchr(this->Txt,kar);
            if(p)
                return p - this->Txt;
            else
                return -1;
            }
        strng * substr(ptrdiff_t pos,ptrdiff_t len = -1) const
            {
            if(len < 0)
                len = strlen(this->Txt+pos);
            char * buf = new char[len+1];
            char * p = buf;
            char * s = this->Txt+pos;
            while(len-- > 0 && *s)
                {
                *p++ = *s++;
                }
            *p = '\0';
            strng * ret = new strng();
            ret->Txt = buf;
            return ret;
            }
        size_t length() const
            {
            return strlen(this->Txt);
            }
        const strng * salad(const strng * L,const strng * R) const
            {// Find out what is sandwiched between L and R
            const char * l = L->Txt;
            const char * r = R->Txt;
            r += strlen(r);
            const char * s = this->Txt;
            const char * e = this->Txt+strlen(this->Txt);
            while(*l && *l == *s)
                {
                ++l;
                ++s;
                }
            assert(!*l);
            while(r > R->Txt && *--e == *--r && *e) 
                {
                }
           // assert(s < e);
            strng * ret = new strng();
            ret->Txt = new char[e - s + 1];
            char * d = ret->Txt;
            while(s < e)
                {
                *d++ = *s++;
                }
            *d = '\0';
            return ret;
            }
        void trim()
            {
            const char * l = this->Txt;
            const char * r = l;
            r += strlen(r);
            while(*l == ' ' || *l == '\t' || *l == '\n' || *l == '\r')
                {
                ++l;
                }
            while(*--r == ' ' || *r == '\t' || *r == '\n' || *r == '\r')
                {
                }
            char * nTxt = new char[r - l + 2];
            char * d = nTxt;
            while(l <= r)
                {
                *d++ = *l++;
                }
            *d = '\0';
            delete this->Txt;
            this->Txt = nTxt;
            }
    };


strng * makeNode(strng ** patreps,int & nr,const strng * pat,const strng * rep,const strng * L,const strng * R,strng ** pLL,strng ** pRR);

#endif