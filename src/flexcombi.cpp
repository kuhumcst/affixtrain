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

buf -4  (Link, OPT, cannot be at start of file, in case of ambiguity) points to next link in ambiguity chain. 
        Special values:
        <=0: This is the last link. A rule follows. 
        4 or -4: There is no tree here, do parent's rewriting.
buf     Used if pattern fails. Number of bytes to 
        (A) next node (if this one fails) or 
        (B) first link in ambiguity chain of fail tree. 
        0 means: stop here       
buf+4   first bit  0: Fail branch is unambiguous, buf points to tree. (A) 
        first bit  1: Fail branch is ambiguous, buf points to chain. (B)
        second bit 0: Success branch is unambiguous, buf+8 points to tree (C)
        second bit 1: Success branch is ambiguous, buf+8 points to chain (D)
buf+8   start of rewrite rule


 (padding to word boundary)
buf+4xN Used if pattern succeeds. Number of bytes to
        (C) next child or
        (D) next link in ambiguity chain

There are no children at all if there is not made room for them. 
That happens if some rule's sibling starts right after this rule,
and also when the end of the file is reached. 

If (C), buf+4xN is new buf
If (D), buf+4xN + 4 is new buf, unless buf+4xN is 4 or -4: then the parent's rewrite rule applies.

Alternative, if buf+4 is not 0 1 2 or 3, then buf+4 is the start of the rewrite rule and there is no ambiguity (equivalent with buf+4 is 0). 

Rewrite rule
------------

buf+8   prefix pattern starts here
(or buf+4)            .
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
\n      list of patterns and replacements ends here
0       padding (opt)
0       padding (opt)
0       padding (opt)

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

static void printpatBracmat(char ** fields,int findex,char * start,char * end
                     ,char * pattern,char * replacement
                     )
    {
    sprintf(start+strlen(start),"%.*s",(int)(fields[1] - fields[0] - 1),fields[0]);
    char * ppat = pattern;
    char * prep = replacement;
    int inc;
    strcpy(pattern,start);
    ppat += strlen(start);
    for(int M = 5;M < findex;M += 2)
        {
        *ppat++ = ANY;
        inc = (int)(fields[M] - fields[M-1] - 1);
        strncpy(ppat,fields[M-1],inc);
        ppat += inc;
        }
    if(findex > 2)
        {
        Strrev(end);
        size_t L = strlen(end);
        sprintf(end+L,"%.*s",(int)(fields[3] - fields[2] - 1),fields[2]);
        Strrev(end+L);
        Strrev(end);
        *ppat++ = ANY;
        }
    strcpy(ppat,end);
 
    inc = (int)(fields[2] - fields[1] - 1);
    strncpy(prep,fields[1],inc);
    prep += inc;
    for(int M = 5;M < findex;M += 2)
        {
        *prep++ = ANY;
        inc = (int)(fields[M+1] - fields[M] - 1);
        strncpy(prep,fields[M],inc);
        prep += inc;
        }
    if(findex > 2)
        {
        inc = (int)(fields[4] - fields[3] - 1);
        *prep++ = ANY;
        strncpy(prep,fields[3],inc);
        prep += inc;
        }
    *prep = '\0';
    }

static void printpat(char ** fields,int findex,char * start,char * end,FILE * fm)
    {
    sprintf(start+strlen(start),"%.*s",(int)(fields[1] - fields[0] - 1),fields[0]);
    fprintf(fm,"%s",start);
    for(int M = 5;M < findex;M += 2)
        {
        fprintf(fm,"*%.*s",(int)(fields[M] - fields[M-1] - 1),fields[M-1]);
        }
    if(findex > 2)
        {
        Strrev(end);
        size_t L = strlen(end);
        sprintf(end+L,"%.*s",(int)(fields[3] - fields[2] - 1),fields[2]);
        Strrev(end+L);
        Strrev(end);
        fprintf(fm,"*%s\t-->\t",end);
        }
    else
        {
        fprintf(fm,"%s\t-->\t",end);
        }

    fprintf(fm,"%.*s",(int)(fields[2] - fields[1] - 1),fields[1]);
    for(int M = 5;M < findex;M += 2)
        {
        fprintf(fm,"*%.*s",(int)(fields[M+1] - fields[M] - 1),fields[M]);
        }
    if(findex > 2)
        {
        fprintf(fm,"*%.*s\n",(int)(fields[4] - fields[3] - 1),fields[3]);
        }
    else
        {
        fprintf(fm,"\n");
        }
    }

static char * printrulesBracmat
                ( char * rules
                , char * max
                , char * start
                , char * end
                , FILE * fmbra
                , strng * L
                , strng * R
                , int & nr
                , int indent
                );




struct fileBuffer
    {
    char * buf;
    long Length;
    fileBuffer() :buf(0), Length(0L){}
    ~fileBuffer(){ delete[] buf; }
    bool readRules(FILE * flexrulefile)
        {
        fseek(flexrulefile, 0, SEEK_END);
        Length = ftell(flexrulefile);
        buf = new char[Length];
        rewind(flexrulefile);
        if (Length != (long int)fread(buf, 1, Length, flexrulefile))
            {
            fclose(flexrulefile);
            return false;
            }
        fclose(flexrulefile);
        return true;
        }

    bool readRules(const char * filename)
        {
        FILE * f = fopen(filename, "rb");
        if (f)
            {
            return readRules(f);
            }
        return false;
        }
    };






static char * printlistBracmat
                ( char * p
                , char * e
                , int indent
                , char * start
                , char * end
                , FILE * fmbra
                , strng * L
                , strng * R
                , int & nr
                , bool comma
                )
    {
    if(p < e)
        {
        if(comma)
            fprintf(fmbra,"%*s",2+indent*2,",");
        size_t stlen = strlen(start);
        size_t endlen = strlen(end);
        while(p < e)
            {
            p = printrulesBracmat(p,e,start,end,fmbra,L,R,nr,indent+1);
            start[stlen] = '\0';
            Strrev(end);
            end[endlen] = 0;
            Strrev(end);
            }
        fprintf(fmbra,"%*s\n",2+indent*2,"");
        }
    return e;
    }



static char * printrules
(char * rules
, char * max
, char * start
, char * end
, FILE * fm
, int indent
);

static char * printChain
                ( char * p
                , char * e
                , int indent
                , char * start
                , char * end
                , FILE * fm
				, char * msg
                )
    {
	int index = *(int *)p;
	if(index > 0)
		{
		fprintf(fm,"%*s( %s\n",indent,"",msg);
		for(;;)
			{
			fprintf(fm,"%*s%s",indent,"","ALT\n");
			if(index == 4 || index == -4)
				fprintf(fm,"%*s",indent,"[parent]\n");
			else
				/*printlist*/printrules
					( p + sizeof(int)
					, (index == 0) ? e : (p + index)
					, start
					, end
					, fm
					, indent
					);

			if(index > 0)
				{
				p += index;
				index = *(int *)p;
				}
			else
				break;
			}
		fprintf(fm,"%*s) END %s\n",indent,"",msg);
		}
    return e;
    }

static char * printrules
                ( char * rules
                , char * max
                , char * start
                , char * end
                , FILE * fm
                , int indent
                )
    {
	if(max <= rules)
		return rules;
    ptrdiff_t index = *(int *)rules;
    if(index == 0)
        index = max - rules;
    char * p = rules + sizeof(int);
    unsigned int type = *(int*)p;
	if(type > 3)
		type = 0;
	else
		p += sizeof(int);
    size_t slen = strlen(start);
    size_t elen = strlen(end);
    char * fields[44];
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
    fprintf(fm,"%*s",indent,"");
    printpat(fields,findex,start,end,fm);
    ptrdiff_t nxt = p - rules;
    nxt += sizeof(int) - 1;
    nxt /= sizeof(int);
    nxt *= sizeof(int);
    p = rules+nxt;
	if(type & 2)
		{ // several chains of children ahead
		p = printChain
            ( p
            , rules + index
            , indent + 2
            , start
            , end
            , fm
			, " ambiguous children"
            );
		}
	else
		{
		p = /*printlist*/printrules
			( p
			, rules + index
			, start
			, end
			, fm
			, indent + 2
			);
		}
    start[slen] = '\0';
    Strrev(end);
    end[elen] = 0;
    Strrev(end);
	if(type & 1)
		{
		p = printChain
            ( rules + index
            , max
            , indent+1
            , start
            , end
            , fm
			, " ambiguous tails of children"
            );
		}
	else
		p = /*printlist*/printrules
			( rules + index
			, max
			, start
			, end
			, fm
			, indent
			);
	assert(p == max);
    return p;
    }

int prettyPrint(const char * flexrulesIn)
    {
    fileBuffer FileBuffer;
    if (!FileBuffer.readRules(flexrulesIn))
        {
        printf("Error (prettyPrint): Cannot open %s for reading\n", flexrulesIn);
        return false;
        }

    char filenameOut[1000];
    sprintf(filenameOut, "%s.pretty.txt", flexrulesIn);

    FILE * fm = fopen(filenameOut, "wb");
    if (!fm)
        {
        printf("Error (prettyPrint): Cannot open %s for writing\n", filenameOut);
        return false;
        }
    char start[1000] = { 0 };
    char end[1000] = { 0 };
    printrules(FileBuffer.buf, FileBuffer.buf + FileBuffer.Length, start, end, fm, 0);
    fclose(fm);
    return true;
    }


static char * printChainBracmat
(char * p
, char * e
, int indent
, char * start
, char * end
, FILE * fmbra
, strng * L
, strng * R
, int & nr
, char * msg
)
    {
    int index = *(int *)p;
    if (index > 0)
        {
        fprintf(fmbra, "%*s( %s.\n", indent, "", msg);
        for (;;)
            {
            fprintf(fmbra, "%*s(ALT.\n", indent, "");
            if (index == 4 || index == -4)
                fprintf(fmbra, "%*s", indent, " parent\n");
            else
                /*printlist*/printrulesBracmat
                (p + sizeof(int)
                , (index == 0) ? e : (p + index)
                , start
                , end
                , fmbra
                , L, R, nr
                , indent
                );
            fprintf(fmbra, "%*s ENDALT)\n", indent, "");

            if (index > 0)
                {
                p += index;
                index = *(int *)p;
                }
            else
                break;
            }
        fprintf(fmbra, "%*s END %s)\n", indent, "", msg);
        }
    return e;
    }


static char * printrulesBracmat
                ( char * rules
                , char * max
                , char * start
                , char * end
                , FILE * fmbra
                , strng * L
                , strng * R
                , int & nr
                , int indent
                )
    {
    if (max <= rules)
        return rules;
    ptrdiff_t index = *(int *)rules;
    if(index == 0)
        index = max - rules;
    char * p = rules + sizeof(int);
    unsigned int type = *(int*)p;
    if (type > 3)
        type = 0;
    else
        p += sizeof(int);
    size_t slen = strlen(start);
    size_t elen = strlen(end);
    fprintf(fmbra,"\n%*s",indent*2,"(\n");
#if 0
    if(type < 4)
        {
        // fork
        fprintf(fmbra,"%*s",2+indent*2,"(\n");
        fprintf(fmbra,"%*s",2+indent*2,"\n");
        static char nix[1] = "";
        char * ret = nix;
        p += sizeof(int); // skip word containing 1, 2 or 3
        start[slen] = 0;
        end[elen] = 0;
        switch(type)
            {
            case 1:
                ret = printlistBracmat(p,max,indent,start,end,fmbra,L,R,nr,false);
                fprintf(fmbra,"%*s\n%*s)\n",2+indent*2,".",2+indent*2,"(");
                break;
            case 2:
                fprintf(fmbra,"%*s)\n%*s\n",2+indent*2,"(",2+indent*2,".");
                ret = printlistBracmat(p,max,indent,start,end,fmbra,L,R,nr,false);
                break;
            case 3:
                {
                char * next = p + *(int *)p;
                p += sizeof(int);
                printlistBracmat(p,next,indent,start,end,fmbra,L,R,nr,false);
                start[slen] = 0;
                end[elen] = 0;
                fprintf(fmbra,"%*s",2+indent*2,".\n");
                ret = printlistBracmat(next,max,indent,start,end,fmbra,L,R,nr,false);
                break;
                }
            }
        fprintf(fmbra,"%*s",2+indent*2,")\n");
        fprintf(fmbra,"\n%*s",indent*2,")\n");
        return ret;
        }
#endif
    char * fields[44];
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
    fprintf(fmbra,"%*s",indent*2,"");
    char pattern[1000];
    char replacement[1000];
    printpatBracmat(fields,findex,start,end,pattern,replacement);
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
    ptrdiff_t nxt = p - rules;
    nxt += sizeof(int) - 1;
    nxt /= sizeof(int);
    nxt *= sizeof(int);
    p = rules+nxt;



    if (type & 2)
        { // several chains of children ahead
        
        p = printChainBracmat
            (p
            , rules + index
            , indent + 2
            , start
            , end
            , fmbra
            , &nL, &nR, nr
            , " ambiguous children"
            );
        }
    else
        {
        p = /*printlist*/printrulesBracmat
            (p
            , rules + index
            , start
            , end
            , fmbra
            , &nL, &nR, nr
            , indent + 2
            );
        }
    start[slen] = '\0';
    Strrev(end);
    end[elen] = 0;
    Strrev(end);
    if (type & 1)
        {
        p = printChainBracmat
            (rules + index
            , max
            , indent + 1
            , start
            , end
            , fmbra
            , &nL, &nR, nr
            , " ambiguous tails of children"
            );
        }
    else
        p = /*printlist*/printrulesBracmat
        (rules + index
        , max
        , start
        , end
        , fmbra
        , &nL, &nR, nr
        , indent
        );

//    p = printlistBracmat(p,rules+index,indent,start,end,fmbra,&nL,&nR,nr,true);
    delete snode;
    fprintf(fmbra,"\n%*s",indent*2,")\n");
    return p;
    }


int prettyPrintBracmat(const char * flexrulesIn)
    {
    fileBuffer FileBuffer;
    if (!FileBuffer.readRules(flexrulesIn))
        {
        printf("Error (prettyPrintBracmat): Cannot open %s for reading\n", flexrulesIn);
        return false;
        }

    char brafile[1000];
    sprintf(brafile, "%s.pretty.bra", flexrulesIn);
    
    FILE * fmbra = fopen(brafile, "wb");
    if (!fmbra)
        {
        printf("Error (prettyPrint): Cannot open %s for writing\n", brafile);
        return false;
        }

    char start[1000] = { 0 };
    char end[1000] = { 0 };
    
    int nr = 0;
    strng L("");
    strng R("");

    printrulesBracmat(FileBuffer.buf, FileBuffer.buf + FileBuffer.Length, start, end, fmbra, &L, &R, nr, 0);

    fclose(fmbra);
    return true;
    }

class rule
    {
    public:
    char * Start;
    char * End;
    rule(char * S, char * E) :Start(S), End(E){}
    void print(int ind)
        {
        char * s = Start;
        int dif = '>' + '|';
        int sep = '>';
        printf("%*s", ind, "");
        while (*s != '\n')
            {
            if (*s == '\t')
                {
                printf("%c",sep);
                sep = dif - sep;
                }
            else
                printf("%c", *s);
            ++s;
            }
        printf("\n");
        }
    int copy(char * arr, int ind)
        {
        strncpy(arr, Start, End - Start);
        printf("arr: %.*s %d\n", End - Start, arr, End - Start);
        return End - Start;
        }
    bool eq(rule * R)
        {
        if (End - Start != R->End - R->Start)
            return false;
        return !strncmp(Start, R->Start, End - Start);
        }
    };

class chain;
class treenode;

class oneOrMore
    {
    public:
    treenode * one;
    chain * more;
    oneOrMore(treenode * un, chain * plus) :one(un), more(plus){}
    treenode * oneOrMore::get(treenode * example);
    void merge(oneOrMore * Y);
    ~oneOrMore();
    };

class treenode
    {
    public:
        oneOrMore Sibling;
        rule Rule;
        oneOrMore Child;
        void print(int ind);
        int copy(char * arr, int ind);
        bool eq(treenode * Y)
            {
            return Rule.eq(&Y->Rule);
            }
        unsigned int type()
            {// 0: unambiguous 1: ambiguous child 2: ambiguous sibling 3: ambiguous child and ambiguous sibling
            int tp = 0;
            if (Sibling.more)
                tp += 1;
            if (Child.more)
                tp += 2;
            return tp;
            }
        treenode(char * rl, char * end, treenode * Sib, chain * ASib, treenode * Chld, chain * AChld)
            :Rule(rl,end),Sibling(Sib,ASib),Child(Chld,AChld){}
        ~treenode();
        void merge(treenode * Y)
            {
            Child.merge(&Y->Child);
            Sibling.merge(&Y->Sibling);
            }
        bool merge(oneOrMore * Y);
    };

treenode * treenodeFactory(char * buf, char * end);

class chain
    {
    public:
        chain * Alt;
        treenode * TreeNode;
        chain(char * buf, char * end)
            {
            int Next = *(int*)buf;
            if (Next > 0)
                Alt = new chain(buf + Next, end);
            else
                { 
                Alt = NULL;
                if (Next < 0)
                    Next = -Next;
                }
            buf += sizeof(int);
            TreeNode = treenodeFactory(buf, (Next == 0) ? end : buf+Next);
            }
        void append(chain * C)
            {
            if (Alt)
                Alt->append(C);
            else
                Alt = C;
            }
        chain(treenode * node, chain * alt) :Alt(alt), TreeNode(node){}
        ~chain(){ delete Alt; delete TreeNode; }
        void merge(oneOrMore * Y)
            {
            if (!TreeNode || !TreeNode->merge(Y))
                {
                if (Alt)
                    {
                    Alt->merge(Y);
                    }
                }
            }
        treenode * get(treenode * example)
            {
            if (Alt)
                {
                if (Alt->TreeNode->eq(example))
                    {
                    chain * N = Alt;
                    treenode * Result = Alt->TreeNode;
                    Alt = Alt->Alt;
                    N->Alt = 0;
                    N->TreeNode = 0;
                    delete N;
                    return Result;
                    }
                else
                    return Alt->get(example);
                }
            else
                return false;
            }
        void print(int ind);
        int copy(char * arr, int ind);
    };

treenode * oneOrMore::get(treenode * example)
    {
    if (more)
        {
        if (more->TreeNode->eq(example))
            {
            chain * N = more;
            treenode * Result = N->TreeNode;
            more = more->Alt;
            N->Alt = 0;
            N->TreeNode = 0;
            delete N;
            return Result;
            }
        else
            return more->get(example);
        }
    else if (one)
        {
        if (one->eq(example))
            {
            treenode * Result = one;
            one = 0;
            return Result;
            }
        else
            return 0;
        }
    else 
        return 0;
    }

void treenode::print(int ind)
    {
    Rule.print(ind);
    if (Sibling.one)
        {
        printf("%*s", ind, "");
        printf("Sibling:\n");
        Sibling.one->print(ind);
        printf("%*s", ind, "");
        printf("END Sibling");
        }
    if (Sibling.more)
        {
        printf("%*s", ind, "");
        printf("Sibling chain:\n");
        Sibling.more->print(ind + 2);
        printf("%*s", ind, "");
        printf("END Sibling chain\n");
        }
    if (Child.one)
        {
        printf("%*s", ind, "");
        printf("Child:\n");
        Child.one->print(ind+2);
        printf("%*s", ind, "");
        printf("END Child\n");
        }
    if (Child.more)
        {
        printf("%*s", ind+2, "");
        printf("Child chain:\n");
        Child.more->print(ind + 4);
        printf("%*s", ind+2, "");
        printf("END Child chain\n");
        }
    }

int treenode::copy(char * arr,int ind)
    {
    int * FailBranch = (int*)arr;
    int type = 0;
    arr += sizeof(int);
    if (Sibling.more)
        type += 1;
    if (Child.more)
        type += 2;
    if (type > 0)
        {
        *(int*)arr = type;
        arr += sizeof(int);
        }
    arr += Rule.copy(arr,ind);
    if (Child.one)
        {
        //	*(int *)arr = 0;
        //	arr += sizeof(int);
        arr += Child.one->copy(arr, ind + 2);
        }
    else if (Child.more)
        {
        arr += Child.more->copy(arr, ind + 4);
        }
    else
        {
        // Nothing!
        }

    if (Sibling.one)
        {
        *FailBranch = arr - (char*)FailBranch;
        arr += Sibling.one->copy(arr, ind);
        }
    else if (Sibling.more)
        {
        *FailBranch = arr - (char*)FailBranch;
        arr += Sibling.more->copy(arr, ind + 2);
        }
    else
        {
        *FailBranch = 0;
        }
    return arr - (char*)FailBranch;
    }

void chain::print(int ind)
    {
    if (TreeNode)
        {
        TreeNode->print(ind + 4);
        }
    else
        {
        printf("%*s(parent)\n", ind + 4, "");
        }
    printf("%*s", ind + 4, "");
    printf("---\n");
    if (Alt)
        Alt->print(ind);
    }

int chain::copy(char * arr, int ind)
    {
    int * Other = (int *)arr;
    arr += sizeof(int);
    if (TreeNode)
        {
        arr += TreeNode->copy(arr,ind);
        *Other = Alt ? arr - (char *)Other : 0;
        }
    else
        {
        *Other = Alt ? sizeof(int) : -(int)sizeof(int);
        }
    if (Alt)
        arr += Alt->copy(arr, ind);
    return arr - (char*)Other;
    }

treenode::~treenode()
    {
    }

oneOrMore::~oneOrMore()
    {
    delete one;
    delete more;
    }


treenode * treenodeFactory(char * buf,char * end)
    {
    if (end < buf+2*sizeof(int))
        return NULL;
    else
        {
        treenode * Sibling = NULL;
        chain * ASibling = NULL;
        treenode * Child = NULL;
        chain * AChild = NULL;

        char * Fail = buf;
        unsigned int OnFail = *(unsigned int*)buf;
        buf += sizeof(unsigned int);
        unsigned int type = *(unsigned int*)buf;
        if (type > 3)
            type = 0;
        else
            buf += sizeof(unsigned int);
		char * ChildEnd = end;
        if (OnFail)
            {
			ChildEnd = Fail + OnFail;
            if (type & 1)
                ASibling = new chain(Fail + OnFail, end);
            else
                Sibling = treenodeFactory(Fail + OnFail, end);
            }
        char * rule = buf;
        char * e = rule;
        while (*e != '\n')
            {
            ++e;
            }
        ptrdiff_t ichild = e - buf;
        ichild += sizeof(int);
        ichild /= sizeof(int);
        ichild *= sizeof(int);
        buf += ichild;
        if (type & 2)
            AChild = new chain(buf, ChildEnd);
        else
            Child = treenodeFactory(buf, ChildEnd);

        switch (type)
            {
                case 0:
                    return new treenode(rule, rule + ichild, Sibling, 0, Child, 0);
                case 1:
                    return new treenode(rule, rule + ichild, 0, ASibling, Child, 0);
                case 2:
                    return new treenode(rule, rule + ichild, Sibling, 0, 0, AChild);
                case 3:
                default:
                    return new treenode(rule, rule + ichild, 0, ASibling, 0, AChild);
            }
        }
    }

bool treenode::merge(oneOrMore * Y)
    {
    if (Y->one)
        {
        if (eq(Y->one))
            {
            merge(Y->one);
            Y->one = 0;
            return true;
            }
        }
    else if (Y->more)
        {
        treenode * Twin = Y->more->get(this);
        merge(Twin);
        return true;
        }
    return false;
    }


void oneOrMore::merge(oneOrMore * Y)
    {
    if (one)
        {
        one->merge(Y);
        if (Y->one)
            {
            more = new chain(one,new chain(Y->one,0));
            one = 0;
            }
        else if (Y->more)
            {
            more = new chain(one, Y->more);
            one = 0;
            }
        }
    else if (more)
        {
        more->merge(Y);
        if (Y->one)
            more->append(new chain(Y->one, 0));
        else if (Y->more)
            {
            more->append(Y->more);
            }
        }
    else
        {
        if (Y->one)
            {
            more = new chain(0, new chain(Y->one,0));
            }
        else if (Y->more)
            {
            more = new chain(0, Y->more);
            }
        }
    }

bool flexcombi(const char * bestflexrules, const char * nextbestflexrules, const char * combinedflexrules)
    {
    fileBuffer FileBuffer;
    fileBuffer NextFileBuffer;
    if (!FileBuffer.readRules(bestflexrules))
        {
        printf("Error (flexcombi): Cannot open %s for reading\n",bestflexrules);
        return false;
        }
    if (!NextFileBuffer.readRules(nextbestflexrules))
        {
        printf("Error (flexcombi): Cannot open %s for reading\n",nextbestflexrules);
        return false;
        }
    char * arr = new char[2 * (FileBuffer.Length + NextFileBuffer.Length)];
    FILE * f = fopen(combinedflexrules,"wb");
    if(!f)
        {
        printf("Error (flexcombi): Cannot open %s for writing\n",combinedflexrules);
        return false;
        }
    treenode * TreeNode = treenodeFactory(FileBuffer.buf, FileBuffer.buf+FileBuffer.Length);
    if (TreeNode)
        {
        printf("\nTree:\n");
        TreeNode->print(0);
        treenode * NextTreeNode = treenodeFactory(NextFileBuffer.buf, NextFileBuffer.buf + NextFileBuffer.Length);
        printf("\nNext Tree:\n");
        NextTreeNode->print(0);
        TreeNode->merge(NextTreeNode);
        printf("\nCombined Tree:\n");
        TreeNode->print(0);
        int length = TreeNode->copy(arr, 0);
        *(int*)arr = 0;
        for(int i = 0;i < length;++i)
			fputc(arr[i],f);
        }
    fclose(f);
	prettyPrint(combinedflexrules);
    prettyPrintBracmat(combinedflexrules);
    return true;
    }
