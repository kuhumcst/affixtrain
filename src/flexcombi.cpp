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

#include "flexcombi.h"
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>


static void Strrev(char * s)
    {
    char * e = s + strlen(s);
    while(s < --e)
        {
        char t = *s;
        *s++ = *e;
        *e = t;
        }
    }


static FILE * fm = NULL;

static void printpat(char ** fields,int findex,char * start,char * end)
    {
    sprintf(start+strlen(start),"%.*s",(int)(fields[1] - fields[0] - 1),fields[0]);
    fprintf(fm,"%s",start);
    for(int m = 1;2*m+3 < findex;++m)
        {
        int M = 2*m+3;
        fprintf(fm,"*%.*s",(int)(fields[M] - fields[M-1] - 1),fields[M-1]);//,vars[m].e - vars[m].s,vars[m].s);
        }
    if(findex > 2)
        {
        Strrev(end);
        size_t L = strlen(end);
        sprintf(end+L,"%.*s",(int)(fields[3] - fields[2] - 1),fields[2]);
        Strrev(end+L);
        Strrev(end);
        fprintf(fm,"*%s\t-->\t",end);
        /*Strrev(end);
        end[L] = 0;
        Strrev(end);*/
        }
    else
        fprintf(fm,"%s\t-->\t",end);

    fprintf(fm,"%.*s",(int)(fields[2] - fields[1] - 1),fields[1]);
    for(int m = 1;2*m+3 < findex;++m)
        {
        int M = 2*m+3;
        fprintf(fm,"*%.*s",(int)(fields[M+1] - fields[M] - 1),fields[M]);//,vars[m].e - vars[m].s,vars[m].s);
        }
    if(findex > 2)
        {
        fprintf(fm,"*%.*s\n",(int)(fields[4] - fields[3] - 1),fields[3]);
        }
    else
        fprintf(fm,"\n");
    }

static char * printrules(char * rules,int indent,char * max,char * start,char * end);

static char * printlist(char * p,char * e,int indent,char * start,char * end)
    {
    size_t stlen = strlen(start);
    size_t endlen = strlen(end);
    while(p < e)
        {
        p = printrules(p,indent+1,e,start,end);
        start[stlen] = '\0';
        Strrev(end);
        end[endlen] = 0;
        Strrev(end);
        }
    return e;
    }

static char * printrules(char * rules,int indent,char * max,char * start,char * end)
    {
    char * fields[44];
    ptrdiff_t index = *(int *)rules;
    if(index == 0)
        index = max - rules;
    char * p = rules + sizeof(int);
    unsigned int byt = *p;
    size_t slen = strlen(start);
    size_t elen = strlen(end);
    if(byt < 4)
        {
        // fork
        fprintf(fm,"%*s",2+indent*2,"(\n");
        static char nix[1] = "";
        char * ret = nix;
        p += sizeof(int); // skip word containing 1, 2 or 3
        start[slen] = 0;
        end[elen] = 0;
        switch(byt)
            {
            case 1:
                ret = printlist(p,max,indent,start,end);
                fprintf(fm,"%*s|\n%*s(parent)\n",2+indent*2,"",2+indent*2,"");
                break;
            case 2:
                fprintf(fm,"%*s(parent)\n%*s|\n",2+indent*2,"",2+indent*2,"");
                ret = printlist(p,max,indent,start,end);
                break;
            case 3:
                {
                char * next = p + *(int *)p;
                p += sizeof(int);
                printlist(p,next,indent,start,end);
                start[slen] = 0;
                end[elen] = 0;
                fprintf(fm,"%*s",2+indent*2,"|\n");
                ret = printlist(next,max,indent,start,end);
                break;
                }
            }
        fprintf(fm,"%*s",2+indent*2,")\n");
        return ret;
        }
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
    fprintf(fm,"%*s",indent*2,"");
    printpat(fields,findex,start,end);

    /*start[slen] = 0;
    end[elen] = 0;*/
    

    ptrdiff_t nxt = p - rules;
    nxt += sizeof(int) - 1;
    nxt /= sizeof(int);
    nxt *= sizeof(int);
    p = rules+nxt;
    return printlist(p,rules+index,indent,start,end);
    }

static bool readRules(FILE * flexrulefile,char *& buf,long int & end)
    {
    fseek(flexrulefile,0,SEEK_END);
    end = ftell(flexrulefile);
    buf = new char[end];
    rewind(flexrulefile);
    if(end != (long int)fread(buf,1,end,flexrulefile))
        {
        fclose(flexrulefile);
	    return false;
        }
    fclose(flexrulefile);
    return true;
    }


bool readRules(const char * filename,char *& buf,long int & end)
    {
    FILE * f = fopen(filename,"rb");
    if(f)
        {
        return readRules(f,buf,end);
        }
    return false;
    }

static void parse(char * buf,char * len,int & sibling,char *& p,char *& e,char *& child/*,char *& field_1,char *& field_3*/,char *& mmax)
    {
    assert(len > buf);
    sibling = *(int *)buf;
    p = buf + sizeof(int);
    e = p;
    //int t = 0;
    while(*e != '\n')
        {
        /*if(*e == '\t')
            {
            ++t;
            if(t == 2)
                field_1 = ++e;
            else if(t == 4)
                field_3 = ++e;
            else
                ++e;
            }
        else*/
            ++e;
        }
    ptrdiff_t ichild = e - buf;
    ichild += sizeof(int);
    ichild /= sizeof(int);
    ichild *= sizeof(int);
    child = buf + ichild;
    mmax = sibling ? buf + sibling : len;
    }

static int merge(char * arr, 
           char * buf1,
           char * len1,
           char * buf2,
           char * len2,
           //char * p,
           int indent
           )
    {
    int sibling1;
    int sibling2;
    char * p1;
    char * p2;
    char * e1;
    char * e2;
    char * child1;
    char * child2;
    /*    
    char * field_11;
    char * field_31;
    char * field_12;
    char * field_32;
    */
    char * mmax1;
    char * mmax2;

    parse(buf1,len1,sibling1,p1,e1,child1,/*field_11,field_31,*/mmax1);
    parse(buf2,len2,sibling2,p2,e2,child2,/*field_12,field_32,*/mmax2);

    char * k;
    int ind = 0;
//    printf("\n[%.*s] =?= [%.*s]\n",e1 - p1,p1,e2 - p2,p2);
    if(  e1 - p1 == e2 - p2 
      && !strncmp(p1,p2,e1-p1)
      )
        { // equal nodes
        // do children
        for(k = buf1;k < child1;++k)
            arr[ind++] = *k;
        //printf("\n%*s%.*s",indent,"",e1-p1,p1);
        char * arrind = arr+ind;
        *(int *)(arrind) = 0xBABE;// test
        if(child1 < mmax1)
            {
            if(child2 < mmax2)
                { // each of the two nodes has at least one child
                ind += merge(arr+ind,child1,mmax1,child2,mmax2,/*p1,*/indent+1);
                }
            else
                { // the second node has no children
                  // copy the children of the first node
                ind += sizeof(int);
                arr[ind] = 1;
                arr[ind+1] = 12;// test
                arr[ind+2] = 13;// test
                arr[ind+3] = 15;// test
                ind += sizeof(int);
                for(k = child1;k < mmax1;++k)
                    arr[ind++] = *k;
                *(int *)(arrind) = ind;
                }
            }
        else
            {
            if(child2 < mmax2)
                { // the first node has no children
                  // copy the children of the second node
                ind += sizeof(int);
                arr[ind] = 2;
                arr[ind+1] = 15;// test
                ind += sizeof(int);
                for(k = child2;k < mmax2;++k)
                    arr[ind++] = *k;
                *(int *)(arrind) = ind;
                }
            else
                ;
            }
        *(int*)(arr) = ind;
        // do siblings
        arrind = arr+ind;
        *(int *)(arrind) = 0xDEAD;// test
        if(sibling1)
            {
            if(sibling2)
                {
                ind += merge(arrind,buf1+sibling1,len1,buf2+sibling2,len2/*,p*//*,field1,field3*/,indent);
                }
            else
                {
                ind += sizeof(int);
                arr[ind] = 1;
                arr[ind+1] = 11;// test
                ind += sizeof(int);
                for(k = buf1+sibling1;k < len1;++k)
                    arr[ind++] = *k;
                *(int *)(arrind) = ind;
                }
            }
        else
            {
            if(sibling2)
                {
                ind += sizeof(int);
                arr[ind] = 2;
                arr[ind+1] = 14; // test
                ind += sizeof(int);
                for(k = buf2+sibling2;k < len2;++k)
                    arr[ind++] = *k;
                *(int *)(arrind) = ind;
                }
            else
                ;
            }
        }
    else // nodes differ
        {
        ind += sizeof(int);
        arr[ind] = 3; //second word: 3 means that there are a lhs and a rhs
        ind += sizeof(int);
        ptrdiff_t offset = len1 - buf1 + sizeof(int);
        *(int *)(arr+ind) = (int)offset; // third word: number of bytes to rhs
        ind += sizeof(int);
        for(k = buf1/*+sibling1*/;k < len1;++k)
            {
            arr[ind++] = *k;
            }
        for(k = buf2/*+sibling2*/;k < len2;++k)
            arr[ind++] = *k;
        *(int *)arr = ind;
        }
    return ind;
    }

int flexcombi(const char * bestflexrules, const char * nextbestflexrules, const char * combinedflexrules)
    {
    char * buf1;
    long int end1;
    char * buf2;
    long int end2;
    char * buf12;
    long int end12;
    if(!readRules(bestflexrules,buf1,end1))
        {
        printf("Error: Cannot open %s for reading\n",bestflexrules);
        return -1;
        }
    char start[1000] = {0};
    char end[1000] = {0};
    char one   [1000];sprintf(one,   "%s.1.txt",bestflexrules);
    char two   [1000];sprintf(two,   "%s.2.txt",nextbestflexrules);
    char onetwo[1000];sprintf(onetwo,"%s.12.txt",combinedflexrules);
    fm = fopen(one,"wb");
    printrules(buf1,0,buf1+end1,start,end);
    fclose(fm);
    if(!readRules(nextbestflexrules,buf2,end2))
        {
        printf("Error: Cannot open %s for reading\n",nextbestflexrules);
        return -1;
        }
    fm = fopen(two,"wb");
    printrules(buf2,0,buf2+end2,start,end);
    fclose(fm);
    //printf("end1 %d + end2 %d = %d\n",end1,end2,end1+end2);
    char * arr = new char[2*(end1 + end2)];
    FILE * f = fopen(combinedflexrules,"wb");
    if(!f)
        {
        printf("Error: Cannot open %s for writing\n",combinedflexrules);
        return -1;
        }
    fm = fopen("merged","wb");
    int length = merge(arr,buf1,buf1+end1,buf2,buf2+end2/*,"\t\t\t\n"*/,0);
    fclose(fm);
    //printf("length %d\n",length);
    *(int*)arr = 0;
    for(int i = 0;i < length;++i)
        fputc(arr[i],f);
    fclose(f);
    if(!readRules(combinedflexrules,buf12,end12))
        return -1;
    //printf("end12 %d\n",end12);
    start[0] = 0;
    end[0] = 0;
    /**/
    printf("\n>>\n");
    fm = fopen(onetwo,"wb");
    printrules(buf12,0,buf12+end12,start,end);
    fclose(fm);
    printf("\n<<\n");
    /**/
    return 0;
    }
