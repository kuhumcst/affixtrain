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

#include "strng.h"
#include "settingsaffixtrain.h"


strng::strng(const char * buf)
    {
    Txt = ::dup(buf);
    ++StrngCount;
    }

strng::strng(const char * buf,size_t length)
    {
    Txt = dupn(buf,length);
    ++StrngCount;
    }


void cut(char * s,char ** a)
    {
    *a = s;
    while(*s && *s != ANY)
        ++s;
    if(*s)
        {
        *s++ = '\0';
        cut(s,++a);
        }
    else
        *++a = 0;
    }

void uncut(char ** a)
    {
    if(*++a)
        {
        (*a)[-1]= ANY;
        uncut(a);
        }
    }

bool compatibleCondition(char ** A,char ** B)
    {
    char * a = *A;
    char * b = *B;
    // Check compatibility at start of pattern
    if(a[0] == START && b[0] == START)
        {
        int i;
        for(i = 1;a[i] && b[i];++i)
            if(a[i] != b[i])
                return false;
        }
    while(A[1])
        ++A;
    while(B[1])
        ++B;
    if(a != *A || b != *B)
        {
        a = *A;
        b = *B;

        // Check compatibility at end of pattern
        size_t alast = strlen(a) - 1;
        size_t blast = strlen(b) - 1;
        if(a[alast] == END && b[blast] == END)
            {
            int i,j;
            for(i = (int)alast - 1,j = (int)blast - 1;i >= 0 && j >= 0;--i,--j)
                if(a[i] != b[j])
                    return false;
            }
        }
    return true;
    }

bool strongerCondition(char ** A,char ** B)
    {
    // Return true if all B are contained in an A.
    // The loop exits when there are no more A's.
    // If there still is a B, then this B must be unmatched, and so the
    // function must return false.
    // If there is no B left, it may still be the case that the last
    // B is not contained in an A, so we also have to check for the
    // result from the last call to strstr.
    char * sub = 0;
    char * a = *A;
    char * b = *B;
    while(b && a)
        {
        sub = 0;
        while(a)
            {
            sub = strstr(a,b);
            if(sub)
                {
                a = sub + strlen(b);
                if(!*a)
                    a = *++A;
                break;
                }
            a = *++A;
            }
        b = *++B;
        }
    return !*B && sub;
    }


edif strng::dif(strng * s) // returns Txt - s
    {
    /* returns 1 if Txt is a stronger condition (further from the root) 
       than s_Txt.
       It returns -1 if Txt is weaker than s_Txt
       It retuns 0 if the conditions are incommensurable.
    */
    char * Txt = this->Txt;
    char * s_Txt = s->Txt;
    char * A[100];
    char * B[100];
    cut(Txt,A);
    cut(s_Txt,B);
    if(compatibleCondition(A,B))
        {
        bool AB = strongerCondition(A,B);
        bool BA = strongerCondition(B,A);

        edif res = (AB && !BA) ? dif_bigger : (BA && !AB) ? dif_smaller : (AB && BA) ? dif_equal : dif_incommensurable;
        uncut(A);
        uncut(B);
        assert(res == dif_incommensurable || res == dif_smaller || res == dif_bigger || res == dif_equal);
        return res;
        }
    else
        {
        uncut(A);
        uncut(B);
        return dif_incompatible;
        }
    }


void strng::checkIntegrity()
    {
    assert(Txt);
    }

bool strng::hasWildCard()
    {
    return Txt[0] != START || Txt[strlen(Txt)-1] != END || strchr(Txt,ANY);
    }


#if 0
edif dif(char * Txt, char * s_Txt)
    {
    /* returns 1 if Txt is a stronger condition (further from the root) 
       than s_Txt.
       It returns -1 if Txt is weaker than s_Txt
       It retuns 0 if the conditions are incommensurable.
    */
    char * A[100];
    char * B[100];
    cut(Txt,A);
    cut(s_Txt,B);
    if(compatibleCondition(A,B))
        {
        bool AB = strongerCondition(A,B);
        bool BA = strongerCondition(B,A);

        edif res = (AB && !BA) ? dif_bigger : (BA && !AB) ? dif_smaller : (AB && BA) ? dif_equal : dif_incommensurable;
        uncut(A);
        uncut(B);
        assert(res == dif_incommensurable || res == dif_smaller || res == dif_bigger || res == dif_equal);
        return res;
        }
    else
        {
        uncut(A);
        uncut(B);
        return dif_incompatible;
        }
    }

#endif


