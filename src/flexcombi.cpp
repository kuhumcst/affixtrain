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

/*
Ambiguous pattern
-----------------

buf      number of bytes to next node to try, if this one fails       
            .
            .
            .
buf+4    0, 1, 2 or 3: ambiguous node. TODO Check tha value 0 cannot occur.
                        TODO change role of 0 to indicate new format.
                        first byte 0, second byte: number of alternatives

buf+8    if 3: number of bytes to second node, if 1 or 2: start of first resp. second pattern (second resp. first pattern is parent)
            .           TODO new style: if number of bytes is 4, no node follows -> use parent node
            .
            .
buf+12   if 3: start of first pattern (buf)



Unambiguous pattern
-------------------

buf      number of bytes to next pattern to try, if this one fails       
            .
            .
            .
buf+4    prefix pattern starts here
            .
            .
\t       
         prefix replacement starts here
            .
            .
            .
\t
         suffix pattern starts here
            .
\t
         suffix replacement starts here
            .
            .
\t
         first infix pattern starts here
\t
         first infix replacement starts here
...
\n       list of patterns and replacements ends here
0        padding (opt)
0        padding (opt)
0        padding (opt)
buf+*(int*)buf    start of subtree
*/

#include "strng.h"
#include "settingsaffixtrain.h"
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


static void printpat(char ** fields,int findex,char * start,char * end,FILE * fm
#if BRACMATOUTPUT
					 ,char * pattern,char * replacement
#endif
					 )
    {
#if BRACMATOUTPUT
	char * ppat = pattern;
	char * prep = replacement;
	int inc;
	strcpy(pattern,start);
	ppat += strlen(start);
#endif
    sprintf(start+strlen(start),"%.*s",(int)(fields[1] - fields[0] - 1),fields[0]);
    fprintf(fm,"%s",start);
    for(int M = 5;M < findex;M += 2)
        {
        fprintf(fm,"*%.*s",(int)(fields[M] - fields[M-1] - 1),fields[M-1]);
#if BRACMATOUTPUT
		*ppat++ = ANY;
		inc = (int)(fields[M] - fields[M-1] - 1);
		strncpy(ppat,fields[M-1],inc);
		ppat += inc;
#endif
        }
    if(findex > 2)
        {
        Strrev(end);
        size_t L = strlen(end);
        sprintf(end+L,"%.*s",(int)(fields[3] - fields[2] - 1),fields[2]);
        Strrev(end+L);
        Strrev(end);
        fprintf(fm,"*%s\t-->\t",end);
#if BRACMATOUTPUT
		*ppat++ = ANY;
		strcpy(ppat,end);
#endif
        }
    else
		{
        fprintf(fm,"%s\t-->\t",end);
#if BRACMATOUTPUT
		strcpy(ppat,end);
#endif
		}

#if BRACMATOUTPUT
	inc = (int)(fields[2] - fields[1] - 1);
	strncpy(prep,fields[1],inc);
	prep += inc;
#endif
    fprintf(fm,"%.*s",(int)(fields[2] - fields[1] - 1),fields[1]);
    for(int M = 5;M < findex;M += 2)
        {
        fprintf(fm,"*%.*s",(int)(fields[M+1] - fields[M] - 1),fields[M]);
#if BRACMATOUTPUT
		*prep++ = ANY;
		inc = (int)(fields[M+1] - fields[M] - 1);
		strncpy(prep,fields[M],inc);
		prep += inc;
#endif
        }
    if(findex > 2)
        {
        fprintf(fm,"*%.*s\n",(int)(fields[4] - fields[3] - 1),fields[3]);
#if BRACMATOUTPUT
		inc = (int)(fields[4] - fields[3] - 1);
		*prep++ = ANY;
		strncpy(prep,fields[3],inc);
		prep += inc;
		*prep = '\0';
#endif
        }
    else
		{
        fprintf(fm,"\n");
#if BRACMATOUTPUT
		*prep = '\0';
#endif
		}
    }

static char * printrules
				( char * rules
				, int indent
				, char * max
				, char * start
				, char * end
				, FILE * fm
#if BRACMATOUTPUT
				, FILE * fmbra
				, strng * L
			    , strng * R
			    , int & nr
#endif
				);

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

int prettyPrint(const char * flexrulesIn,const char * filenameOut)
    {
    char * buf;
    long int endPos;
    if(!readRules(flexrulesIn,buf,endPos))
        {
        printf("Error (prettyPrint): Cannot open %s for reading\n",flexrulesIn);
        return false;
        }
    char start[1000] = {0};
    char end[1000] = {0};
    FILE * fm = fopen(filenameOut,"wb");
    if(!fm)
        {
        printf("Error (prettyPrint): Cannot open %s for writing\n",filenameOut);
        return false;
        }
#if BRACMATOUTPUT
	char brafile[1000];
	sprintf(brafile,"%s.bra",filenameOut);
    FILE * fmbra = fopen(brafile,"wb");

	
    int nr = 0;
    strng L("");
    strng R("");
	
	printrules(buf,0,buf+endPos,start,end,fm
		,fmbra,&L,&R,nr
		);
    if(fmbra)
		fclose(fmbra);
#else
	printrules(buf,0,buf+endPos,start,end,fm);
#endif
    fclose(fm);
    return true;
    }


static char * printlist
				( char * p
				, char * e
				, int indent
				, char * start
				, char * end
				, FILE * fm
#if BRACMATOUTPUT
				, FILE * fmbra
				, strng * L
			    , strng * R
			    , int & nr
				, bool comma
#endif
				)
    {
	if(p < e)
		{
#if BRACMATOUTPUT
		if(comma)
			fprintf(fmbra,"%*s",2+indent*2,",");
#endif
		size_t stlen = strlen(start);
		size_t endlen = strlen(end);
		while(p < e)
			{
			p = printrules(p,indent+1,e,start,end,fm
#if BRACMATOUTPUT
				,fmbra,L,R,nr
#endif
				);
			start[stlen] = '\0';
			Strrev(end);
			end[endlen] = 0;
			Strrev(end);
			}
#if BRACMATOUTPUT
		fprintf(fmbra,"%*s\n",2+indent*2,"");
#endif
		}
    return e;
    }

static char * printrules
				( char * rules
				, int indent
				, char * max
				, char * start
				, char * end
				, FILE * fm
#if BRACMATOUTPUT
				, FILE * fmbra
				, strng * L
			    , strng * R
			    , int & nr
#endif
				)
    {
    char * fields[44];
    ptrdiff_t index = *(int *)rules;
    if(index == 0)
        index = max - rules;
    char * p = rules + sizeof(int);
    unsigned int byt = *p;
    size_t slen = strlen(start);
    size_t elen = strlen(end);
#if BRACMATOUTPUT
    fprintf(fmbra,"\n%*s",indent*2,"(\n");
#endif
    if(byt < 4)
        {
        // fork
#if BRACMATOUTPUT
        fprintf(fmbra,"%*s",2+indent*2,"(\n");
        fprintf(fmbra,"%*s",2+indent*2,"\n");
#endif
        fprintf(fm,"%*s",2+indent*2,"(\n");
        static char nix[1] = "";
        char * ret = nix;
        p += sizeof(int); // skip word containing 1, 2 or 3
        start[slen] = 0;
        end[elen] = 0;
        switch(byt)
            {
            case 1:
                ret = printlist(p,max,indent,start,end,fm
#if BRACMATOUTPUT
					,fmbra,L,R,nr,false
#endif
					);
                fprintf(fm,"%*s|\n%*s(parent)\n",2+indent*2,"",2+indent*2,"");
#if BRACMATOUTPUT
                fprintf(fmbra,"%*s\n%*s)\n",2+indent*2,".",2+indent*2,"(");
#endif
                break;
            case 2:
                fprintf(fm,"%*s(parent)\n%*s|\n",2+indent*2,"",2+indent*2,"");
#if BRACMATOUTPUT
				fprintf(fmbra,"%*s)\n%*s\n",2+indent*2,"(",2+indent*2,".");
#endif
                ret = printlist(p,max,indent,start,end,fm
#if BRACMATOUTPUT
					,fmbra,L,R,nr,false
#endif
					);
                break;
            case 3:
                {
                char * next = p + *(int *)p;
                p += sizeof(int);
                printlist(p,next,indent,start,end,fm
#if BRACMATOUTPUT
					,fmbra,L,R,nr,false
#endif
					);
                start[slen] = 0;
                end[elen] = 0;
                fprintf(fm,"%*s",2+indent*2,"|\n");
#if BRACMATOUTPUT
                fprintf(fmbra,"%*s",2+indent*2,".\n");
#endif
                ret = printlist(next,max,indent,start,end,fm
#if BRACMATOUTPUT
					,fmbra,L,R,nr,false
#endif
					);
                break;
                }
            }
        fprintf(fm,"%*s",2+indent*2,")\n");
#if BRACMATOUTPUT
        fprintf(fmbra,"%*s",2+indent*2,")\n");
	    fprintf(fmbra,"\n%*s",indent*2,")\n");
#endif
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
#if BRACMATOUTPUT
    fprintf(fmbra,"%*s",indent*2,"");
	char pattern[1000];
	char replacement[1000];
    printpat(fields,findex,start,end,fm,pattern,replacement);
	strng spattern(pattern);
	strng sreplacement(replacement);
	strng * nLL = 0;
	strng * nRR = 0;
	strng * patreps[100];
	unsigned int i;
	for(i = 0;i < sizeof(patreps)/sizeof(patreps[0]);++i)
		patreps[i] = 0;
	strng * snode = makeNode(patreps,nr,&spattern,&sreplacement,L,R,&nLL,&nRR);
	fprintf(fmbra,"(%s)",snode->itsTxt());
	strng nL(L);
	strng nR(nRR);
	nL.cat(nLL,(const strng *)0);
	nR.cat(R,(const strng *)0);
	delete nLL;
	delete nRR;
#else
    printpat(fields,findex,start,end,fm);
#endif
    ptrdiff_t nxt = p - rules;
    nxt += sizeof(int) - 1;
    nxt /= sizeof(int);
    nxt *= sizeof(int);
    p = rules+nxt;
#if BRACMATOUTPUT
    p = printlist(p,rules+index,indent,start,end,fm,fmbra,&nL,&nR,nr,true);
	delete snode;
    fprintf(fmbra,"\n%*s",indent*2,")\n");
#else
    p = printlist(p,rules+index,indent,start,end,fm);
#endif
	return p;
    }


static void parse(char * buf,char * len,int & sibling,char *& p,char *& e,char *& child,char *& mmax)
    {
    assert(len > buf);
    sibling = *(int *)buf;
    p = buf + sizeof(int);
    e = p;
    while(*e != '\n')
        {
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
    char * mmax1;
    char * mmax2;

    parse(buf1,len1,sibling1,p1,e1,child1,mmax1);
    parse(buf2,len2,sibling2,p2,e2,child2,mmax2);

    char * k;
    int ind = 0;
    if(  e1 - p1 == e2 - p2 
      && !strncmp(p1,p2,e1-p1)
      )
        { // equal nodes
        // do children
        for(k = buf1;k < child1;++k)
            arr[ind++] = *k;
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
                ind += merge(arrind,buf1+sibling1,len1,buf2+sibling2,len2,indent);
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
        for(k = buf1;k < len1;++k)
            {
            arr[ind++] = *k;
            }
        for(k = buf2;k < len2;++k)
            arr[ind++] = *k;
        *(int *)arr = ind;
        }
    return ind;
    }

bool flexcombi(const char * bestflexrules, const char * nextbestflexrules, const char * combinedflexrules)
    {
    char * buf1;
    long int end1;
    char * buf2;
    long int end2;
    if(!readRules(bestflexrules,buf1,end1))
        {
        printf("Error (flexcombi): Cannot open %s for reading\n",bestflexrules);
        return false;
        }
    if(!readRules(nextbestflexrules,buf2,end2))
        {
        printf("Error (flexcombi): Cannot open %s for reading\n",nextbestflexrules);
        return false;
        }
    char * arr = new char[2*(end1 + end2)];
    FILE * f = fopen(combinedflexrules,"wb");
    if(!f)
        {
        printf("Error (flexcombi): Cannot open %s for writing\n",combinedflexrules);
        return false;
        }
    int length = merge(arr,buf1,buf1+end1,buf2,buf2+end2,0);
    *(int*)arr = 0;
    for(int i = 0;i < length;++i)
        fputc(arr[i],f);
    fclose(f);
    return true;
    }
