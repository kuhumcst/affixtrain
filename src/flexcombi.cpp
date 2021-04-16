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

Remark
------
Merging trees in the implemented way does not always preserve the order of the
results. In other words, lemmatising with trees T1, T2 and T3 (in that order),
giving results L1, L2 and L3, is not necessarily equivalent with lemmatising
with the merged tree, which could give result e.g., L3, L1 and L2.
The difference is more pronounced the more results are produced . 
Even a pair of ambiguous results can be in the 'wrong' order.
*/

#include "flexcombi.h"
#include "strng.h"
#include "affixtrain.h"
#include "settingsaffixtrain.h"
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>


const char * prettytxt = NULL;
const char * prettybra = NULL;

static char start[1000] = { 0 };
static char end[1000] = { 0 };
static FILE * fm = 0;
static FILE * fmbra = 0;
static int indent = 0;

typedef unsigned char typetype;

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

static void printpatBracmat(char ** fields,int findex/*,char * start,char * end*/
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
        strncpy(ppat, fields[M - 1], (size_t)inc);
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
    strncpy(prep, fields[1], (size_t)inc);
    prep += inc;
    for(int M = 5;M < findex;M += 2)
        {
        *prep++ = ANY;
        inc = (int)(fields[M+1] - fields[M] - 1);
        strncpy(prep, fields[M], (size_t)inc);
        prep += inc;
        }
    if(findex > 2)
        {
        inc = (int)(fields[4] - fields[3] - 1);
        *prep++ = ANY;
        strncpy(prep, fields[3], (size_t)inc);
        prep += inc;
        }
    *prep = '\0';
    }

static void printpat(char ** fields,int findex/*,char * start,char * end,FILE * fm*/)
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

static void printrulesBracmat
                ( char * rules
                , char * max
//                , char * start
//                , char * end
//                , FILE * fmbra
                , strng * L
                , strng * R
                , int & nr
//                , int indent
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
        long length = ftell(flexrulefile);
        if(length > 4)
            {
            rewind(flexrulefile);

            char V[5] = "";
            V[4] = 0;
            if (4 == fread(V, 1, 4, flexrulefile))
                {
                if ((V[0] & 1) && (V[3] & 1))
                    {
                    if (strcmp(V, "\rV3\r")) // Version string at beginning of file. If read as int, must have lowest bit set! Therefore \r (00001101).
                        return false;
                    length -= 4;
                    }
                else
                    {
                    rewind(flexrulefile);
                    }
                }

            Length = length;
            while (Length & 3) ++Length;
            buf = new char[(size_t)Length];
            if (length == (long int)fread(buf, 1, (size_t)length, flexrulefile))
                {
                return true;
                }
            }
        return false;
        }

    bool readRules(const char * filename)
        {
        FILE * f = fopen(filename, "rb");
        ++openfiles;
        if (f)
            {
            bool read = readRules(f);
            --openfiles;
            fclose(f);
            return read;
            }
        return false;
        }
    };





static char * printrules
(char * rules
, char * max
//, char * start
//, char * end
//, FILE * fm
//, int indent
);

static char * printChain
                ( char * p
                , char * e
//                , int indent
//                , char * start
//                , char * end
//                , FILE * fm
                , const char * msg
                )
    {
    int index = *(int *)p;
    assert(!(index & 3));
    if(index > 0)
        {
        fprintf(fm,"%*s( %s\n",indent,"",msg);
        for(;;)
            {
            fprintf(fm,"%*s%s",indent,"","ALT\n");
            if(index == 4 || index == -4)
                fprintf(fm,"%*s",indent,"[parent]\n");
            else
                printrules
                    ( p + sizeof(int)
                    , (index == 0) ? e : (p + index)
//                    , start
//                    , end
//                    , fm
//                    , indent
                    );

            if(index > 0)
                {
                p += index;
                index = *(int *)p;
                assert(!(index & 3));
                }
            else
                break;
            }
        fprintf(fm,"%*s) END %s\n",indent,"",msg);
        }
    return e;
    }

static char * printrules
(char * rules
 , char * max
 //, char * start
 //, char * end
 //, FILE * fm
 //, int indent
)
    {
    for (;;)
        {
        if (max <= rules)
            {
            assert(rules == max);
            return rules;
            }
        else
            {
            ptrdiff_t index = *(int *)rules;
            assert(!(index & 3));
            if (index == 0)
                index = max - rules;
            assert(!(index & 3));
            char * p = rules + sizeof(int);
            typetype type = *(typetype*)p;
            if (type > 3)
                type = 0;
            else
                p += sizeof(typetype);
            size_t slen = strlen(start);
            size_t elen = strlen(end);
            char * fields[44];
            fields[0] = p;
            int findex = 1;
            while (*p != '\n')
                {
                if (*p == '\t')
                    fields[findex++] = ++p;
                else
                    ++p;
                }
            fields[findex] = ++p; // p is now within 3 bytes from the next record.
            fprintf(fm, "%*s", indent, "");
            printpat(fields, findex/*, start, end, fm*/);
            ptrdiff_t nxt = p - rules;
            nxt += sizeof(int) - 1;
            nxt /= sizeof(int);
            nxt *= sizeof(int);
            p = rules + nxt;
            indent += 2;
            if (type & 2)
                { // several chains of children ahead
                p = printChain
                (p
                 , rules + index
//              , indent + 2
//              , start
//              , end
//              , fm
                 , " ambiguous children"
                );
                }
            else
                {
                p = printrules
                (p
                 , rules + index
//              , start
//              , end
//              , fm
//              , indent + 2
                );
                }
            indent -= 2;
            start[slen] = '\0';
            Strrev(end);
            end[elen] = 0;
            Strrev(end);
            if (type & 1)
                {
                ++indent;
                p = printChain
                (rules + index
                 , max
//              , indent + 1
//              , start
//              , end
//              , fm
                 , " ambiguous tails of children"
                );
                --indent;
                assert(p == max);
                return p;
                }
            else
                {
                /*
                p = printrules
                (rules + index
                , max
                , start
                , end
                , fm
                , indent
                );
                assert(p == max);
                return p;
                */
                rules += index;
                }
            }
        }
    }

int prettyPrint(const char * flexrulesIn)
    {
    fileBuffer FileBuffer;
    if (!FileBuffer.readRules(flexrulesIn))
        {
        fprintf(stderr, "Error (prettyPrint): Cannot open %s (flexrulesIn) for reading\n", flexrulesIn);
        return false;
        }

    char filenameOut[1000];
    if(prettytxt == NULL)
        prettytxt = ".pretty.txt";
    sprintf(filenameOut, "%s%s", flexrulesIn,prettytxt);

    /*FILE * */ fm = fopen(filenameOut, "wb");
    ++openfiles;
    if (!fm)
        {
        fprintf(stderr, "Error (prettyPrint): Cannot open %s (filenameOut) for writing\n", filenameOut);
        return false;
        }
//    char start[1000] = { 0 };
//    char end[1000] = { 0 };
    indent = 0;
    printrules(FileBuffer.buf, FileBuffer.buf + FileBuffer.Length/*, start, end, fm, 0*/);
    --openfiles;
    fclose(fm);
    return true;
    }


static void printChainBracmat
        ( char * p
        , char * e
//        , int indent
//        , char * start
//        , char * end
//        , FILE * fmbra
        , strng * L
        , strng * R
        , int & nr
        )
    {
    int index = *(int *)p;
    assert(!(index & 3));
    if (index > 0)
        {
        fprintf(fmbra, "%*s(\n", indent, "");
        for (;;)
            {
            fprintf(fmbra, "%*s(\n", indent, "");
            if (index == 4 || index == -4)
                fprintf(fmbra, "%*s", indent, " parent\n");
            else
                printrulesBracmat
                    ( p + sizeof(int)
                    , (index == 0) ? e : (p + index)
//                    , start
//                    , end
//                    , fmbra
                    , L, R, nr
//                    , indent
                    );
            fprintf(fmbra, "%*s )\n", indent, "");

            if (index > 0)
                {
                fprintf(fmbra, " . ");
                p += index;
                index = *(int *)p;
                assert(!(index & 3));
                }
            else
                break;
            }
        fprintf(fmbra, "%*s)\n", indent, "");
        }
    }


static void printrulesBracmat
                ( char * rules
                , char * max
//                , char * start
//                , char * end
//                , FILE * fmbra
                , strng * L
                , strng * R
                , int & nr
//                , int indent
                )
    {
    int parens = 0;
    while(max > rules)
        {
        ptrdiff_t index = *(int *)rules;
        assert(!(index & 3));
        if (index == 0)
            index = max - rules;
        assert(!(index & 3));
        char * p = rules + sizeof(int);
        typetype type = *(typetype*)p;
        if (type > 3)
            type = 0;
        else
            p += sizeof(typetype);
        size_t slen = strlen(start);
        size_t elen = strlen(end);
        fprintf(fmbra,"\n%*s",indent*2,"(\n");
        ++parens;
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
        printpatBracmat(fields,findex,/*start,end,*/pattern,replacement);
        strng spattern(pattern);
        strng sreplacement(replacement);
        strng * nLL = 0;
        strng * nRR = 0;
        strng * patreps[100];
        unsigned int i;
        for(i = 0;i < sizeof(patreps)/sizeof(patreps[0]);++i)
            patreps[i] = 0;
        strng * snode = makeNode(patreps,nr,&spattern,&sreplacement,L,R,&nLL,&nRR);
        fprintf(fmbra,"((%s) ",snode->itsTxt());
        delete snode;
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

        indent += 2;
        if (type & 2)
            { // several chains of children ahead
            if(p < rules + index)
                {
                fprintf(fmbra,",");
                printChainBracmat
                    ( p
                    , rules + index
//                    , indent + 2
//                    , start
//                    , end
//                    , fmbra
                    , &nL, &nR, nr
                    //, " ambiguous children"
                    );
                }
            }
        else
            {
            if(p < rules + index)
                {
                fprintf(fmbra,",");
                printrulesBracmat
                    ( p
                    , rules + index
//                    , start
//                    , end
//                    , fmbra
                    , &nL, &nR, nr
//                    , indent + 2
                    );
                }
            }
        indent -= 2;
        fprintf(fmbra," ) ");
        start[slen] = '\0';
        Strrev(end);
        end[elen] = 0;
        Strrev(end);
        rules += index;
        if (type & 1)
            {
            ++indent;
            printChainBracmat
                ( rules
                , max
//                , indent + 1
//                , start
//                , end
//                , fmbra
                , L, R, nr
                //, " ambiguous tails of children"
                );
            --indent;
            break;
            }
        }

    while(--parens >= 0)
        fprintf(fmbra,"\n%*s",indent*2,")\n");
    }


int prettyPrintBracmat(const char * flexrulesIn)
    {
    fileBuffer FileBuffer;
    if (!FileBuffer.readRules(flexrulesIn))
        {
        fprintf(stderr, "Error (prettyPrintBracmat): Cannot open %s (flexrulesIn) for reading\n", flexrulesIn);
        return false;
        }

    char brafile[1000];
    if(prettybra == NULL)
        prettybra = ".pretty.bra";
    sprintf(brafile, "%s%s", flexrulesIn,prettybra);
    
    /*FILE * */fmbra = fopen(brafile, "wb");
    ++openfiles;
    if (!fmbra)
        {
        fprintf(stderr, "Error (prettyPrintBracmat): Cannot open %s (brafile) for writing\n", brafile);
        return false;
        }

//    char start[1000] = { 0 };
//    char end[1000] = { 0 };
    
    int nr = 0;
    strng L("");
    strng R("");

    indent = 0;
    printrulesBracmat(FileBuffer.buf, FileBuffer.buf + FileBuffer.Length/*, start, end, fmbra*/, &L, &R, nr/*, 0*/);

    --openfiles;
    fclose(fmbra);
    return true;
    }

class rule
    {
    public:
    char * Start;
    char * End;
    rule(char * S, char * E) :Start(S), End(E){}
    void print(FILE * fp,int ind)
        {
        char * s = Start;
        int dif = '>' + '|';
        int sep = '>';
        fprintf(fp, "%*s", ind, "");
        while (*s != '\n')
            {
            if (*s == '\t')
                {
                fprintf(fp, "%c",sep);
                sep = dif - sep;
                }
            else
                fprintf(fp, "%c", *s);
            ++s;
            }
        fprintf(fp, "\n");
        }
    ptrdiff_t copy(char * arr/*, int ind*/)
        {
        ptrdiff_t n = End - Start;
        int p[1];
        ptrdiff_t diff = (ptrdiff_t)(((char*)p - arr - n) % sizeof(int));
        strncpy(arr, Start, (size_t)n);
        return n + diff;
        }
    bool eq(rule * R)
        {
        if (End - Start != R->End - R->Start)
            return false;
        return !strncmp(Start, R->Start, (size_t)(End - Start));
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
    treenode * get(treenode * example);
    void setSib(treenode * Sib) { one = Sib; }
    void merge(oneOrMore * Y);
    ~oneOrMore();
    };

class treenode
    {
    public:
        rule Rule;
        oneOrMore Sibling;
        oneOrMore Child;
        void print(FILE * fp,int ind);
        ptrdiff_t copy(char * arr, int ind);
        bool eq(treenode * Y)
            {
            return Rule.eq(&Y->Rule);
            }
        typetype type()
            {// 0: unambiguous 1: ambiguous child 2: ambiguous sibling 3: ambiguous child and ambiguous sibling
            int tp = 0;
            if (Sibling.more)
                tp += 1;
            if (Child.more)
                tp += 2;
            return (typetype)tp;
            }
        treenode(char * rl, char * anend, chain * ASib, treenode * Chld, chain * AChld)
            :Rule(rl,anend),Sibling(0,ASib),Child(Chld,AChld){}
        void setSib(treenode * Sib) { Sibling.setSib(Sib); }
        ~treenode();
        void merge(treenode * Y)
            {
            Child.merge(&Y->Child);
            Sibling.merge(&Y->Sibling);
            }
        bool merge(oneOrMore * Y);
    };

void treenodeFactory(char * buf, char * end, treenode ** NodePointerNeedingValue);

class chain
    {
    public:
        chain * Alt;
        treenode * TreeNode;
        chain(char * buf, char * anend)
            {
            int Next = *(int*)buf;
            if (Next > 0)
                Alt = new chain(buf + Next, anend);
            else
                { 
                Alt = NULL;
                if (Next < 0)
                    Next = -Next;
                }
            buf += sizeof(int);
            /*TreeNode =*/ treenodeFactory(buf, (Next == 0) ? anend : buf+Next,&TreeNode);
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
                return 0;
            }
        void print(FILE * fp,int ind);
        ptrdiff_t copy(char * arr, int ind);
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

void treenode::print(FILE * fp, int ind)
    {
    Rule.print(fp,ind);
    if (Sibling.one)
        {
        fprintf(fp, "%*s", ind, "");
        fprintf(fp, "Sibling:\n");
        Sibling.one->print(fp,ind);
        fprintf(fp, "%*s", ind, "");
        fprintf(fp, "END Sibling");
        }
    if (Sibling.more)
        {
        fprintf(fp, "%*s", ind, "");
        fprintf(fp, "Sibling chain:\n");
        Sibling.more->print(fp,ind + 2);
        fprintf(fp, "%*s", ind, "");
        fprintf(fp, "END Sibling chain\n");
        }
    if (Child.one)
        {
        fprintf(fp, "%*s", ind, "");
        fprintf(fp, "Child:\n");
        Child.one->print(fp,ind+2);
        fprintf(fp, "%*s", ind, "");
        fprintf(fp, "END Child\n");
        }
    if (Child.more)
        {
        fprintf(fp, "%*s", ind+2, "");
        fprintf(fp, "Child chain:\n");
        Child.more->print(fp,ind + 4);
        fprintf(fp, "%*s", ind+2, "");
        fprintf(fp, "END Child chain\n");
        }
    }

ptrdiff_t treenode::copy(char * arr,int ind)
    {
    int * FailBranch = (int*)arr;
    typetype type = 0;
    arr += sizeof(int);
    if (Sibling.more)
        type += 1;
    if (Child.more)
        type += 2;
    if (type > 0)
        {
        *(typetype*)arr = type;
        arr += sizeof(typetype);
        }
    arr += Rule.copy(arr/*, ind*/);
    if (Child.one)
        {
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
        *FailBranch = (int)(arr - (char*)FailBranch);
        arr += Sibling.one->copy(arr, ind);
        }
    else if (Sibling.more)
        {
        *FailBranch = (int)(arr - (char*)FailBranch);
        arr += Sibling.more->copy(arr, ind + 2);
        }
    else
        {
        *FailBranch = 0;
        }
    return arr - (char*)FailBranch;
    }

void chain::print(FILE * fp,int ind)
    {
    if (TreeNode)
        {
        TreeNode->print(fp,ind + 4);
        }
    else
        {
        fprintf(fp, "%*s(parent)\n", ind + 4, "");
        }
    fprintf(fp, "%*s", ind + 4, "");
    fprintf(fp, "---\n");
    if (Alt)
        Alt->print(fp,ind);
    }

ptrdiff_t chain::copy(char * arr, int ind)
    {
    int * Other = (int *)arr;
    arr += sizeof(int);
    if (TreeNode)
        {
        arr += TreeNode->copy(arr,ind);
        *Other = Alt ? (int)(arr - (char *)Other) : 0;
        }
    else
        {
        *Other = (int)(Alt ? sizeof(int) : -(int)sizeof(int));
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


void treenodeFactory(char * buf,char * anend, treenode ** NodePointerNeedingValue)
    {
    assert(NodePointerNeedingValue != 0);
    for (;;)
        {
        assert(((anend - buf) & 3) == 0);
        if (anend < buf + 2 * sizeof(int))
            {
            *NodePointerNeedingValue = NULL;
            return;
            }
        else
            {
            //treenode * Sibling = NULL;
            chain * ASibling = NULL;
            treenode * Child = NULL;
            chain * AChild = NULL;

            char * Fail = buf;
            unsigned int OnFail = *(unsigned int*)buf;
            buf += sizeof(unsigned int);
            typetype type = *(typetype*)buf;
            if (type > 3)
                type = 0;
            else
                buf += sizeof(typetype);

            char * rule = buf;
            char * e = rule;
            while (*e != '\n')
                {
                ++e;
                }

            ptrdiff_t ichild = e - Fail;
            ichild += sizeof(int);
            ichild /= sizeof(int);
            ichild *= sizeof(int);

            buf = Fail + ichild;

            if (OnFail == 0)
                {
                if (type & 2)
                    {
                    AChild = new chain(buf, anend);
                    *NodePointerNeedingValue = new treenode(rule, e + 1, 0, 0, AChild);
                    return;
                    }
                else
                    {
                    treenodeFactory(buf, anend, &Child);
                    *NodePointerNeedingValue = new treenode(rule, e + 1, 0, Child, 0);
                    return;
                    }
                }
            else
                {
                if (type & 1)
                    {
                    if (type & 2)
                        {
                        AChild = new chain(buf, Fail + OnFail);
                        ASibling = new chain(Fail + OnFail, anend);
                        *NodePointerNeedingValue = new treenode(rule, e + 1, ASibling, 0, AChild);
                        return;
                        }
                    else
                        {
                        treenodeFactory(buf, Fail + OnFail, &Child);
                        ASibling = new chain(Fail + OnFail, anend);
                        *NodePointerNeedingValue = new treenode(rule, e + 1, ASibling, Child, 0);
                        return;
                        }
                    }
                else
                    {
                    if (type & 2)
                        {
                        Child = 0;
                        AChild = new chain(buf, Fail + OnFail);
                        }
                    else
                        {
                        treenodeFactory(buf, Fail + OnFail, &Child);
                        AChild = 0;
                        }
                    treenode * ret = new treenode(rule, e + 1, 0, Child, AChild);
                    *NodePointerNeedingValue = ret;
                    NodePointerNeedingValue = &ret->Sibling.one;
                    buf = Fail + OnFail;
                    }
                }
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

/* NOTE: if cutoff (-c option) is > 0, then the value of combinedflexrules is
not the name of the final result. Instead, for each cutoff level <n> from 0 to
the value of the -c option a directory <n> is created containing the flexrules
pruned to the cutoff level <n>. Each of these flexrule files has the same file
name. See the trainRules function in affixtrain.cpp.
*/
bool flexcombi(const char * bestflexrules, const char * nextbestflexrules, const char * combinedflexrules)
    {
    fileBuffer FileBuffer;
    fileBuffer NextFileBuffer;
    if (!FileBuffer.readRules(bestflexrules))
        {
        fprintf(stderr, "Error (flexcombi): Cannot open %s for reading\n",bestflexrules);
        return false;
        }

    char * arr;
    if (strcmp(bestflexrules, nextbestflexrules))
        {
        if (!NextFileBuffer.readRules(nextbestflexrules))
            {
            fprintf(stderr, "Error (flexcombi): Cannot open %s for reading\n", nextbestflexrules);
                return false;
            }
        arr = new char[(size_t)(2 * (FileBuffer.Length + NextFileBuffer.Length))];
        }
    else
        arr = new char[(size_t)(2 * (FileBuffer.Length))];

    FILE * f = fopen(combinedflexrules, "wb");
    ++openfiles;
    if(!f)
        {
        fprintf(stderr, "Error (flexcombi): Cannot open %s for writing\n",combinedflexrules);
        return false;
        }
    fprintf(f,"\r%.2s\r","V31");
    treenode * TreeNode = NULL;
    treenodeFactory(FileBuffer.buf, FileBuffer.buf + FileBuffer.Length, &TreeNode);
    if (TreeNode)
        {
        if (strcmp(bestflexrules, nextbestflexrules))
            {
            treenode* NextTreeNode = NULL;
            treenodeFactory(NextFileBuffer.buf, NextFileBuffer.buf + NextFileBuffer.Length, &NextTreeNode);
            TreeNode->merge(NextTreeNode);
            }
        ptrdiff_t length = TreeNode->copy(arr, 0);
        *(int*)arr = 0;
        for(ptrdiff_t i = 0;i < length;++i)
            fputc(arr[i],f);
        }
    --openfiles;
    fclose(f);
    return true;
    }
