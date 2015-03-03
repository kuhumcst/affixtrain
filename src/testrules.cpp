#include "testrules.h"
#include "applyaffrules.h"
#include "settingsaffixtrain.h"
#include "affixtrain.h"
#include "optionaff.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef WIN32
#define HOME ""
#ifdef ATHOME
#define BASEDIR ""
#define RESOURCEDIR ""
#else
#define BASEDIR ""
#define RESOURCEDIR ""
#endif
#define BINDIR ""
#define SLASH "\\"
#define QUOTE "\"" 
#else
#define HOME ""
#define BASEDIR HOME ""
#define RESOURCEDIR HOME ""
#define BINDIR HOME ""
#define SLASH "/" 
#define QUOTE "\'" 
#endif
#define SEPARATOR "--"



#define ISO 0
#define OnCE 0
#define ENSURELEMMAPRESENTASFULLFORM 0 // If set to 1 creates doublets! Set to 0 and run again on output from first run to remove doublets.
#define FORCESAMECASING 0
#if FORCESAMECASING
#include "utf8func.h"
#endif
#define CLUMPSENSITIVE 0

#define N_FOLDCROSSVALIDATION 0//1
#if N_FOLDCROSSVALIDATION
#define N_ 10
#endif

static bool REMOVE;
static bool TAGGED;
static bool COUNTLEMMAS;
static bool WRITEWEIRD;
static int NOTWEIRD; // if less than NOTWEIRD characters in the full form
              // and the lemma are the same, the word is written
              // to the file weird_xx.txt, where xx is the language code.

static int FRACTION_LOW; // min 0
static int FRACTION_HIGH; // max 10000
static int CUTOFF_LOW;
static int CUTOFF_HIGH;
static int MAXCOUNT;
static const int STEPSIZE = 100;

static bool TRAINTEST;// = false;// Set to true if the lemmatisation test has to be done on the training data

static int globmaxcount = 0;

static const char * XTRf(optionStruct * Options)
    {
    static char xtr[10];
    sprintf(xtr,"C%d%s%s%s"
        ,Options->expectedCutoff()
        ,Options->redo() ? "R" : ""
        ,Options->suffixOnly() ? "S" : ""
        ,Options->doweights() ? "W" : ""
        );
    return xtr;
    }

static const char * LGf(optionStruct * Options)
    {
    return Options->extra();
    }

static const char * lemmalistef(optionStruct * Options)
    {
    return Options->wordList();
    }

static const char * reducedlemmalistef(optionStruct * Options)
    {
    static char rwl[1000];
    sprintf(rwl,"%s.reduced",Options->wordList());
    return rwl;
    }

static const char * newStyleRulesf(optionStruct * Options)
    {
    static char nsr[1000];
    sprintf(nsr,"%s.%s.newstyle",Options->wordList(),Options->extra());
    return nsr;
    }

static char tr_pathname[256];

class lineab
    {
    private:
        int N;
        double Sxy,Sx,Sy,Sx2;
    public:
        lineab()
            {
            N = 0;
            Sxy = Sx = Sy = Sx2 = 0;
            }

        void add(double rules,double trainsize)
            {
            double drules = log(rules);
            double dtrainsize = log(trainsize);
            if(dtrainsize == 0.0)
                fprintf(stderr,"add(rules %f,trainsize %f); dtrainsize %f\n",rules,trainsize,dtrainsize);
            Sx += dtrainsize;
            Sx2 += dtrainsize*dtrainsize;
            Sy += drules;
            Sxy += dtrainsize*drules;
            ++N;
            }

        double b()
            {
            if((double)N*Sx2 > Sx*Sx)
                return ((double)N*Sxy-Sy*Sx)/((double)N*Sx2-Sx*Sx);
            else
                {
                fprintf(stderr,"N %d Sx2 %f Sx %f N*Sx2 - Sx*Sx %f\n",N,Sx2,Sx,(double)N*Sx2-Sx*Sx);
                return 0;
                }
            }

        double a()
            {
            if(N > 0)
                {
                double B = b();
                return (Sy-Sx*B)/(double)N;
                }
            else
                return 0;
            }
    };

#define CUTOFFS 10

static char Lemmas[256];
static char Weird[256];

static void init(bool TrainTest,const char *& XT,const char *& TT)
    {
    REMOVE = false;
    TAGGED = false;
    COUNTLEMMAS = true;
    WRITEWEIRD = true;
    NOTWEIRD = 2; // if less than NOTWEIRD characters at the beginning of
    // the full form and the lemma are the same, the word is written
    // to the file weird_xx.txt, where xx is the language code.
#if N_FOLDCROSSVALIDATION
    FRACTION_LOW = 10000 - 10000/N_;
    FRACTION_HIGH = FRACTION_LOW;
#else
    switch(STEPSIZE)
        {
        case 100:
            FRACTION_LOW  = 154;//9856;//9000;//9856;//154;//154;//9700;//9000;//9900; // min 0
            //FRACTION_LOW  =       9856;//9000;//9856;//154;//154;//9700;//9000;//9900; // min 0
            break;
        case 25:
            FRACTION_LOW  = 171;
            break;
        default:
            FRACTION_LOW  = 1000;
            ;
        }

    FRACTION_HIGH = 10000;//9900; // max 10000
#endif
    CUTOFF_LOW = 0;
    CUTOFF_HIGH = 5;
    MAXCOUNT = 100;


    // Set to 0 if the lemmatisation test has to be done with words that are not in the training data
    TRAINTEST = TrainTest;
    
    if(TRAINTEST)
        {
        XT = "_train";
        FRACTION_LOW = 10000;
        FRACTION_HIGH = 10000;
        }
    else
        {
        XT = "";
        }

    if(TAGGED)
        {
        TT = "_tag";
        }
    else
        TT = "";
    }

class Count
    {
    public:
        Count * next;
        long refcnt;
        long cnt;
        long correct;
        long wrong;
        void print(FILE *fp,double divisor)
            {
            if(refcnt == 0)
                fprintf(fp,"%ld %8ld %14.6f %14.6f %2.2f %14.6f %14.6f\n",refcnt,cnt,(double)correct/divisor,(double)wrong/divisor,100.0*(double)correct/((double)correct+(double)wrong),(double)correct/divisor,(double)wrong/divisor);
            else            
                fprintf(fp,"%ld %ld %14.6f %14.6f %2.2f %14.6f %14.6f\n",refcnt,cnt,(double)correct/divisor,(double)wrong/divisor,100.0*(double)correct/((double)correct+(double)wrong),(double)correct/(divisor*cnt*refcnt),(double)wrong/(divisor*cnt*refcnt));
            if(next)
                next->print(fp,divisor);
            }
        Count(long refcnt):next(NULL),refcnt(refcnt),cnt(1),correct(0),wrong(0)
            {
            }
        ~Count()
            {
            delete next;
            }
        void Correct(long Cnt)
            {
            if(Cnt == refcnt)
                {
                correct++;
                }
            else if(next)
                {
                next->Correct(Cnt);
                }
            }
        void Wrong(long Cnt)
            {
            if(Cnt == refcnt)
                {
                wrong++;
                }
            else if(next)
                {
                next->Wrong(Cnt);
                }
            }
    };

static Count * counts = NULL;

void Correct(long Cnt)
    {
    if(counts)
        {
        counts->Correct(Cnt);
        }
    }

void Wrong(long Cnt)
    {
    if(counts)
        {
        counts->Wrong(Cnt);
        }
    }

static void trim(char * s)
    {
    char * d = s + strlen(s);
    while(d > s && (*--d == ' ' || *d == '\t' || *d == '\r' || *d == '\n'))
        *d = '\0';
    }

struct line
    {
    bool ambiguous:1;
    char * s;
    line():ambiguous(false),s(NULL){}
    ~line(){delete s;}
    };

static int cmpamb(const void * a,const void * b)
    {
    line * A = *(line **)a;
    line * B = *(line **)b;
    char * as = A->s;
    char * bs = B->s;
    while(*as && *as == *bs && *as != '\t')
        {
        ++as;
        ++bs;
        }
    if(*as == '\t' && *bs == '\t')
        {
        A->ambiguous = true;
        B->ambiguous = true;
        }
    return strcmp(A->s,B->s);
    }

static int mystrcmp(const void * a,const void * b)
    {
    return strcmp(((line *)a)->s, ((line * )b)->s);
    }

#if CLUMPSENSITIVE
struct clump
    {
#if N_FOLDCROSSVALIDATION
    int start;
    int ceil:31;
    unsigned int done:1;
#else
    int start;
    int ceil;
#endif
    };
#endif

static int fileRead(line * lines,
#if CLUMPSENSITIVE
                    clump clumps[],
#endif
                    FILE * fpi,int columnfull,int columnbase,int columnPOS,int sep)
    {
    int linecnt = 0;
#if CLUMPSENSITIVE
    int clumpcnt = 0;
#endif
    int kar;
    int oldkar = 0;
    char buf[256];
    char * pbuf = buf;
    while((kar = fgetc(fpi)) != EOF)
        {
        if(kar == '\n')
            {
#if CLUMPSENSITIVE
            if(oldkar == kar)
                { // clump found
                if(clumpcnt > 0)
                    clumps[clumpcnt].start = clumps[clumpcnt-1].ceil;
                clumps[clumpcnt].ceil = linecnt;
#if N_FOLDCROSSVALIDATION
                clumps[clumpcnt].done = 0;
#endif
                ++clumpcnt;
                }
            else
#endif
                {
                *pbuf++ = '\0';
                pbuf = buf;
                trim(buf);
                if(*buf)
                    {
                    char * col[3] = {NULL,NULL,NULL};
                    col[0] = buf;
                    col[1] = strchr(buf,sep/*'\t'*/);
                    if(col[1])
                        {
                        *col[1] = '\0';
                        ++col[1];
                        col[2] = strchr(col[1],sep/*'\t'*/);
                        if(col[2])
                            {
                            *col[2] = '\0';
                            ++col[2];
                            }
                        }
                    char * f = columnfull ? col[columnfull - 1] : NULL;
                    char * b = columnbase ? col[columnbase - 1] : NULL;
                    const char * cf = f;
                    const char * cb = b;
                    if(f && b)
                        {
                        if(!strcmp(b,"="))
                            {
                            b = f;
                            cb = b;
                            }
                        char * t = columnPOS ? col[columnPOS - 1] : NULL;
#if FORCESAMECASING
                        size_t len = 0;
                        if(isUpperUTF8(b))
                            {
                            if(isAllUpper(b,0))
                                {
                                if(!isAllUpper(f,0))
                                    {
                                    cf = changeCase(f,false,len);
                                    }
                                }
                            else
                                {
                                if(!isUpperUTF8(f))
                                    {
                                    cf = changeCase(f,false,len);
                                    static char * bu = NULL;
                                    if(bu)
                                        delete [] bu;
                                    bu = new char[len+len];
                                    strcpy(bu,cf);
                                    bool UTF8 = true;
                                    const char *pcf = bu;
                                    getUTF8char(pcf,UTF8);
                                    len = 0;
                                    cf = changeCase(pcf,true,len);
                                    strcpy(bu+(pcf-bu),cf);
                                    cf = bu;
                                    }
                                }
                            }
                        else
                            {
                            if(isUpperUTF8(f))
                                {
                                cf = changeCase(f,true,len);
                                }
                            }
#endif
                        if(*cf && *cb && !strchr(cf,' ') && !strchr(cb,' '))
                            {
                            if(TAGGED)
                                {
                                if(*t)
                                    {
                                    lines[linecnt].s = new char[strlen(cf)+strlen(cb)+strlen(t)+3];
                                    sprintf(lines[linecnt].s,"%s\t%s\t%s",cf,cb,t);
                                    lines[linecnt].ambiguous = false;
                                    ++linecnt;
                                    }
                                }
                            else
                                {
                                lines[linecnt].s = new char[strlen(cf+1)+strlen(cb)+2];
                                sprintf(lines[linecnt].s,"%s\t%s",cf+1,cb);
                                lines[linecnt].ambiguous = true;
                                lines[linecnt].s = new char[strlen(cf)+strlen(cb)+2];
                                sprintf(lines[linecnt].s,"%s\t%s",cf,cb);
                                lines[linecnt].ambiguous = false;
                                ++linecnt;
                                }
                            }
                        }
                    else
                        {
                        printf("Format error in line %d\n",linecnt);
                        return -1;
                        }
                    }
                }
            }
        else
            *pbuf++ = (char)kar;
        oldkar = kar;
        }
    line ** plines = new line * [linecnt];
    int i;
    for(i = 0;i < linecnt;++i)
        {
        plines[i] = lines+i;
        }
    qsort(plines,linecnt,sizeof(line *),cmpamb);
    delete [] plines;
    int ambi = 0;
    for(i = 0;i < linecnt;++i)
        {
        if(lines[i].ambiguous)
            {
            //printf("%s\n",lines[i].s);
            ++ambi;
            }
        }
    //printf("%d/%d = %14.6f%%\n",ambi,linecnt,(100.0*ambi)/linecnt);
    return linecnt;
    }

static int removeDuplicateLines(int linecnt,line * lines)
    {
    int i;
    int j;
    for(i = 0,j = i + 1;j < linecnt;++j)
        {
        if(strcmp(lines[i].s,lines[j].s))
            {
            ++i;
            if(i != j)
                {
                assert(lines[i].s == NULL);
                lines[i] = lines[j];
                lines[j].s = 0;
                }
            }
        else
            {
            //printf("lines[%d] %s == lines[%d] %s\n",i,lines[i],j,lines[j]);
            delete lines[j].s;
            lines[j].s = 0;
            }
        }
    return i+1;
    }

static char * trailnumber(char * s,int & number)
/* Removes a number at the end of the line (dict_de_without_doubles_random).
   If there are more than one number, separated by tab (STOok.txt), 
   only the last is removed. */
    {
    char * pi = s + strlen(s);
    number = 0;
    int pos = 1;
    while(--pi > s && '0' <= *pi && *pi <= '9')
        {
        number += (*pi - '0')*pos;
        pos *= 10;
        }
    //assert(*pi != '\t'); //How can there be such? If the word-lemma pair consists of numbers! (res.sl/trainingpairs.txt
    if(*pi == '\t')
        {
        *pi = '\0';
        return pi;
        }
    return NULL;
    }

static int removeDuplicateLinesTrailingNumber(int linecnt,line * lines)
    {
    int i;
    int j;
    int numberi = 0;
    char * pi = trailnumber(lines[0].s,numberi);
    // pi != NULL if a tab has been removed
    for(i = 0,j = i + 1;j < linecnt;++j)
        {
        int numberj = 0;
        char * pj = trailnumber(lines[j].s,numberj);
        if(strcmp(lines[i].s,lines[j].s))
            {
            if(pi)
                {
                *pi = '\t';
                pi = NULL;
                }
            ++i;
            if(i != j)
                {
                assert(lines[i].s == NULL);
                lines[i] = lines[j];
                lines[j].s = 0;
                /*delete [] lines[i];
                lines[i] = lines[j];*/
                }
            numberi = numberj;
            pi = pj;
            }
        else //if(pi && pj)
            {
            if(pj)
                *pj = '\t';
            delete [] lines[j].s;
            lines[j].s = 0;
            /*lines[i] = lines[j];
            pi = pj;*/
            }
        }
    if(pi)
        *pi = '\t';
    return i+1;
    }

static void randomix(line * lines,int linecnt)
    {
    int i;
    line tmp;
    for(i = 0;i < linecnt;++i)
        {
        int j = (int)((double)rand() * (double)linecnt / (double)RAND_MAX);
        if(j != i && j > 0 && j < linecnt)
            {
            tmp = lines[i];
            lines[i] = lines[j];
            lines[j] = tmp;
            }
        }
    tmp.s = NULL;
    }

#if CLUMPSENSITIVE
static void randomix(clump * clumps,int clumpcnt)
    {
    int i;
    for(i = 0;i < clumpcnt;++i)
        {
        int j = (int)((double)rand() * (double)clumpcnt / (double)RAND_MAX);
        if(j != i && j > 0 && j < clumpcnt)
            {
            int start = clumps[i].start;
            int ceil = clumps[i].ceil;
            clumps[i] = clumps[j];
            clumps[j].start = start;
            clumps[j].ceil = ceil;
            }
        }
    }
#endif

static void countLemmas()
    {
    FILE * fplemma = fopen(Lemmas,"rb");
    ++openfiles;
    if(fplemma)
        {
        int kar;
        int linecnt = 0;
        char buf[256];
        char * pbuf = buf;
        //int i;
        //int j;
        while((kar = fgetc(fplemma)) != EOF)
            {
            if(kar == '\n')
                ++linecnt;
            }
        rewind(fplemma);
        line * lines2;
        lines2 = new line[linecnt];
        linecnt = 0;
        while((kar = fgetc(fplemma)) != EOF)
            {
            if(kar == '\n')
                {
                *pbuf++ = '\0';
                pbuf = buf;
                trim(buf);
                if(*buf)
                    {
                    lines2[linecnt].s = new char[strlen(buf)+1];
                    strcpy(lines2[linecnt].s,buf);
                    ++linecnt;
                    }
                }
            else
                *pbuf++ = (char)kar;
            }
        --openfiles;
        fclose(fplemma);
        qsort(lines2,linecnt,sizeof(line),mystrcmp);
        //printf("linecnt lemma %d\n",linecnt);
        int newlinecnt = removeDuplicateLines(linecnt,lines2);
        if(newlinecnt != linecnt)
            {
            linecnt = newlinecnt;
            //printf(" After removing duplicate lines linecnt = %d",linecnt);
            }
        fplemma = fopen(Lemmas,"wb");
        ++openfiles;
        if(fplemma)
            {
            int i;
            for(i = 0;i < linecnt;++i)
                {
                fprintf(fplemma,"%s\n",lines2[i].s);
                }
            --openfiles;
            fclose(fplemma);
            fplemma = NULL;
            }
        delete [] lines2;
        //            getchar();
        }
    else
        fprintf(stderr,"Cannot open \"%s\" for reading\n",Lemmas);
    }


static int readlines(int columnfull,int columnbase,int columnPOS/*,bool TAGGED*/,line *& lines,int sep
#if CLUMPSENSITIVE
                     ,clump *& clumps,int & clumpcnt
#endif
                     ,optionStruct * Options
                     )
    {
    FILE * fpi = NULL;
    int kar;
    int linecnt = 0;
    fpi = fopen(lemmalistef(Options),"r");
    ++openfiles;
    if(!fpi)
        {
        fprintf(stderr,"Cannot open \"%s\" for reading\n",lemmalistef(Options));
        return 0;
        }

    //srand( (unsigned)time( NULL ) );
    srand(1); // 20090513, makes results independent of the order in which the languages are processed
    int oldkar = 0;
    while((kar = fgetc(fpi)) != EOF)
        {
        if(kar == '\n')
            {
            if(oldkar == '\n')
                {
#if CLUMPSENSITIVE
                ++clumpcnt;
#endif
                }
            else
                ++linecnt;
            }
        oldkar = kar;
        }
    assert(!lines);
    rewind(fpi);
    /*
#if CLUMPSENSITIVE
    printf("lines %d clumps %d\n",linecnt,clumpcnt);
#else
    printf("lines %d\n",linecnt);
#endif
    */
    lines = new line[linecnt];
#if CLUMPSENSITIVE
    if(clumpcnt)
        {
        clumps = new clump[clumpcnt];
        clumps[0].start = 0;
        }
    else
        clumps = NULL;
#endif
    //printf("fileRead. Open files:%d\n",::openfiles);

    linecnt = fileRead(lines,
#if CLUMPSENSITIVE
        clumps,
#endif
        fpi,columnfull,columnbase,columnPOS,sep);
    if(linecnt < 0)
        {
        printf("Error in file \"%s\"\n",lemmalistef(Options));
        exit(-1);
        }
    //printf("fileRead DONE\n");
    --openfiles;
    fclose(fpi);
#if CLUMPSENSITIVE
    if(clumpcnt < 2)
#endif
        {
        //printf("clumpcnt < 2\n");
        /*
        FILE * ff = fopen("liste.wri","wb");
        ++openfiles;
        if(ff)
            {
            for(int h = 0;h < linecnt;++h)
                {
                fprintf(ff,"%s\n",lines[h].s ? lines[h].s : "NULLLL");
                }
            --openfiles;
            fclose(ff);
            }
        else
            printf("Cannot open %s for writing\n","liste.wri");
        */
        qsort(lines,linecnt,sizeof(line),mystrcmp);
        //printf("linecnt %d\n",linecnt);
        int newlinecnt = removeDuplicateLines(linecnt,lines);
        if(newlinecnt != linecnt)
            {
            linecnt = newlinecnt;
            //printf(" After removing duplicate lines linecnt = %d\n",linecnt);
            }

        newlinecnt = removeDuplicateLinesTrailingNumber(linecnt,lines);
        assert(newlinecnt == linecnt);
        if(newlinecnt != linecnt)
            {
            linecnt = newlinecnt;
            //printf(" After removing duplicate lines followed by number linecnt = %d\n",linecnt);
            }
        FILE * fplemma = NULL;
        if(COUNTLEMMAS)
            {
            fplemma = fopen(Lemmas,"wb");
            ++openfiles;
            if(!fplemma)
                COUNTLEMMAS = false;
            }

        FILE * fpweird = NULL;
        if(WRITEWEIRD)
            {
            fpweird = fopen(Weird,"wb");
            ++openfiles;
            if(!fpweird)
                WRITEWEIRD = false;
            }
            
        FILE * fpo = fopen(reducedlemmalistef(Options),"wb");
        ++openfiles;
        if(fpo)
            {
            int i;
            for(i = 0;i < linecnt;++i)
                {
                fprintf(fpo,"%s\n",lines[i].s);
#if ENSURELEMMAPRESENTASFULLFORM
                char * ptab = strchr(lines[i].s,'\t');
                assert(ptab);
                fprintf(fpo,"%s%s\n",ptab+1,ptab);
#endif
                if(WRITEWEIRD)
                    {
                    if(columnfull > 0 && columnbase > 0)
                        {
                        char * col[3];
                        col[0] = lines[i].s;
                        col[1] = strchr(lines[i].s,'\t');
                        assert(col[1]);
                        //if(col[1])
                            {
                            ++col[1];
                            if(columnfull > 2 || columnbase > 2)
                                {
                                col[2] = strchr(col[1],'\t');
                                if(col[2])
                                    {
                                    ++col[2];
                                    }
                                }
                            char * b = col[columnbase - 1];
                            char * f = col[columnfull - 1];
                            assert(f < b);
                            int j = 0;

                            while(*b && *b == *f)
                                {
                                ++j;
                                ++b;
                                ++f;
                                }
                            if(j < NOTWEIRD)
                                fprintf(fpweird,"%s\n",lines[i].s);
                            
                            }
                        //else 
                          //  fprintf(fpweird,"%s\n",lines[i].s);
                        }
                    }
                if(COUNTLEMMAS)
                    {
                    char * tab1 = strchr(lines[i].s,'\t');
                    if(tab1)
                        {
                        ++tab1;
                        char * tab2 = strchr(tab1,'\t');
                        if(tab2)
                            {
                            *tab2 = '\0';
                            fprintf(fplemma,"%s\n",tab1);
                            *tab2 = '\t';
                            }
                        else
                            fprintf(fplemma,"%s\n",tab1);
                        }
                    }
                
                }
            --openfiles;
            fclose(fpo);
            fpo = NULL;
            if(WRITEWEIRD)
                {
                --openfiles;
                fclose(fpweird);
                fpweird = NULL;
                }
            if(COUNTLEMMAS)
                {
                assert(fplemma);
                //if(fplemma)
                --openfiles;
                fclose(fplemma);
                fplemma = NULL;
                countLemmas();
                }   
            }
    //    getchar();
        if(fplemma)
            {
            --openfiles;
            fclose(fplemma);
            fplemma = NULL;
            }
        if(fpweird)
            {
            --openfiles;
            fclose(fpweird);
            fpweird = NULL;
            }

        randomix(lines,linecnt);
        assert(fplemma == NULL);
        //if(fplemma)
          //  fclose(fplemma);
        assert(fpweird == NULL);
        //if(fpweird)
          //  fclose(fpweird);
        }
#if CLUMPSENSITIVE
    else
        {
        printf("randomix(%d)\n",clumpcnt);
        randomix(clumps,clumpcnt);
        }
#endif
    return linecnt;
    }


static bool IsSpace(int k)
    {
    return k == ' ' || k == '\t' || k == '\r' || k == '\n' || !k;
    }

static void writeTestAndControl(line * Line,FILE * ftest,FILE * fcontrol,int columntest, int columncontrol, int columntag)
    {
    bool slashprinted = false;
    bool printstop = false;
    int kar;
    int m = 0;
    int col = 1;
    char * s = Line->s;
    if(Line->ambiguous)
        fputc('|',fcontrol);
    while((kar = s[m++]) != '\0')
        {
        if(kar == '\t')
            {
            ++col;
            printstop = false;
            }
        else if(kar == ',')
            {
            printstop = true;
            }
        else if(!printstop)
            {
            if(col == columntest)
                fputc(kar,ftest);
            else if(col == columncontrol)
                fputc(kar,fcontrol);
            else if(TAGGED && col == columntag)
                {
                if(!slashprinted)
                    {
                    slashprinted = true;
                    fputc('/',ftest);
                    fputc('/',fcontrol);
                    }
                fputc(kar,ftest);
                fputc(kar,fcontrol);
                }
            }
        }
    fputc('\n',ftest);
    fputc('\n',fcontrol);
    }


static void printTrain(line * Line,FILE * fptrain,FILE * fptraintest,FILE * fptraincontrol,int columntest,int columncontrol,int columntag)
    {
    if(!IsSpace(Line->s[0]))
        {
        fprintf(fptrain,"%s\n",Line->s);
        if(TRAINTEST)
            {
            writeTestAndControl(Line,fptraintest,fptraincontrol,columntest,columncontrol,columntag);
            }
        }
    else
        {
        printf("%s\n",Line->s);
        }
    }


#if CLUMPSENSITIVE
static int splitLemmaliste
                   (int clumpcnt
                   ,char * training
                   ,char * test
                   ,char * control
                   ,int columntest
                   ,int columncontrol
                   ,int columntag
                   ,int fraction
                   ,char * traintest,char * traincontrol
                   ,clump * clumps
                   ,line * lines
                   ) // 0 <= fraction <= 10000
    {
    FILE * fptrain = NULL;
    FILE * fptest = NULL;
    FILE * fpcontrol = NULL;
    FILE * fptraintest = NULL;
    FILE * fptraincontrol = NULL;
    fptrain = fopen(training,"wb");
    ++openfiles;
//    int lines = 1000000;
    int trainlines = 0;
    int trainclumps = 0;
    int testclumps = 0;
    trainclumps = (int)((double)clumpcnt * (double)fraction / 10000.0);

    testclumps = clumpcnt - trainclumps;
    if(TRAINTEST)
        {
        fptraintest = fopen(traintest,"wb");
        ++openfiles;
        fptraincontrol = fopen(traincontrol,"wb");
        ++openfiles;
        }

    if(clumpcnt && fptrain)
        {
//        int kar;
//        int col = 1;
        //        bool slashprinted = false;
        //        bool printstop = false;
        fptest = fopen(test,"wb");
        ++openfiles;
        fpcontrol = fopen(control,"wb");
        ++openfiles;
        if(fptest && fpcontrol && (!TRAINTEST || (fptraintest && fptraincontrol)))
            {
            while(testclumps > 0 || trainclumps > 0)
                {
                int k = 0;
                while(k < clumpcnt && (testclumps > 0 || trainclumps > 0))
                    {
                    int rd = rand();
                    if(  testclumps == 0 
                      || (  trainclumps > 0
#if N_FOLDCROSSVALIDATION
                         && !clumps[k].done 
#endif
                         && (double)rd * 10000.0 < (double)RAND_MAX * (double)fraction
                         )
                      ) // write to training file
                        {
                        assert(trainclumps > 0);
                        for(int L = clumps[k].start;L < clumps[k].ceil;++L)
                            printTrain(lines+L,fptrain,fptraintest,fptraincontrol,columntest,columncontrol,columntag);
                        trainlines += clumps[k].ceil - clumps[k].start;
#if N_FOLDCROSSVALIDATION
                        clumps[k].done = 1;
#endif
                        --trainclumps;
                        }
                    else
                        {
                        for(int L = clumps[k].start;L < clumps[k].ceil;++L)
                            writeTestAndControl(lines+L,fptest,fpcontrol,columntest,columncontrol,columntag);
                        --testclumps;
                        }
                    ++k;
                    }
                }
            --openfiles;
            fclose(fptest);     fptest = NULL;
            --openfiles;
            fclose(fpcontrol);  fpcontrol = NULL;
            }
        }
    if(fptest){--openfiles;fclose(fptest);}
    if(fpcontrol){--openfiles;fclose(fpcontrol);}
    if(fptraintest){--openfiles;fclose(fptraintest);}
    if(fptraincontrol){--openfiles;fclose(fptraincontrol);}
    if(fptrain){--openfiles;fclose(fptrain);
    return trainlines;
    }
#endif

static int splitLemmaliste
                   (int linecnt
                   ,char * training
                   ,char * test
                   ,char * control
                   ,int columntest
                   ,int columncontrol
                   ,int columntag
                   ,int fraction
                   ,char * traintest,char * traincontrol
                   ,line * lines
                   ) // 0 <= fraction <= 10000
    {
    FILE * fptrain = NULL;
    FILE * fptest = NULL;
    FILE * fpcontrol = NULL;
    FILE * fptraintest = NULL;
    FILE * fptraincontrol = NULL;
    fptrain = fopen(training,"wb");
    ++openfiles;
//    int lines = 1000000;
    int trainlines = 0;
    int ret;
    int testlines = 0;
    if(fraction < 10000)
        trainlines = (int)((double)linecnt * (double)fraction / 10000.0);
    else
        trainlines = linecnt;

    testlines = linecnt - trainlines;
    ret = trainlines;
    if(TRAINTEST)
        {
        fptraintest = fopen(traintest,"wb");
        ++openfiles;
        fptraincontrol = fopen(traincontrol,"wb");
        ++openfiles;
        }

    if(linecnt && fptrain)
        {
        if(fraction < 10000)
            {
            fptest = fopen(test,"wb");
            ++openfiles;
            fpcontrol = fopen(control,"wb");
            ++openfiles;
            if(fptest && fpcontrol && (!TRAINTEST || (fptraintest && fptraincontrol)))
                {
                while(testlines > 0 || trainlines > 0)
                    {
                    int k = 0;
                    while(k < linecnt && (testlines > 0 || trainlines > 0))
                        {
                        assert(!IsSpace(lines[k].s[0]));
                        int rd = rand();
                        if(  testlines == 0 
                          || (  trainlines > 0
                             && (double)rd * 10000.0 < (double)RAND_MAX * (double)fraction && trainlines > 0
                             )
                          ) // write to training file
                            {
                            assert(trainlines > 0);
                            printTrain(lines+k,fptrain,fptraintest,fptraincontrol,columntest,columncontrol,columntag);
                            trainlines--;
                            }
                        else
                            {
                            assert(testlines > 0);
                            if(!IsSpace(lines[k].s[0]))
                                {
                                writeTestAndControl(lines+k,fptest,fpcontrol,columntest,columncontrol,columntag);
                                }
                            testlines--;
                            }
                        ++k;
                        }
                    }
                --openfiles;
                fclose(fptest);     fptest = NULL;
                --openfiles;
                fclose(fpcontrol);  fpcontrol = NULL;
                }
            }
        else
            {
            int k = 0;
            while(k < linecnt)
                {
                printTrain(lines+k,fptrain,fptraintest,fptraincontrol,columntest,columncontrol,columntag);
                ++k;
                }
            }
        }
    if(fptest)          {--openfiles;fclose(fptest);}
    if(fpcontrol)       {--openfiles;fclose(fpcontrol);}
    if(fptraintest)     {--openfiles;fclose(fptraintest);}
    if(fptraincontrol)  {--openfiles;fclose(fptraincontrol);}
    if(fptrain)         {--openfiles;fclose(fptrain);}
    return ret;
    }

typedef int ambty[3]; // counts number of times that first lemma is correct, second lemma is correct and third or higher lemma is correct

struct decision
    {
    int true_amb;       //tp
    int false_amb;      //fp
    int false_not_amb;  //fn
    int true_not_amb;   //tn
    decision():true_amb(0),false_amb(0),false_not_amb(0),true_not_amb(0){}
    double tp()
        {
        return 0.5*true_amb;
        }
    double fp()
        {
        return false_amb;
        }
    double fn()
        {
        return 0.5*false_not_amb;
        }
    double tn()
        {
        return true_not_amb;
        }
    double relevant()
        {
        return tp() + fn();
        }
    double retrieved()
        {
        return tp() + fp();
        }
    double precision()
        {
        return retrieved() > 0.0 ? tp()/retrieved() : 0.0;
        }
    double recall()
        {
        return tp()/relevant();
        }
    };

static void compare(const char * output,const char * control,int & same,int & different,ambty & ambiguous,const char * controlResult,const char * test,int & ambiguousRules,decision & Decision)
    {
    FILE * f1 = fopen(output,"r");
    ++openfiles;
    FILE * f2 = fopen(control,"r");
    ++openfiles;
    FILE * f3 = fopen(controlResult,"w");
    ++openfiles;
    FILE * f4 = fopen(test,"r");
    ++openfiles;
    if(f1 && f2 && f3 && f4)
        {
        char b1[256];
        char b2[256];
        char b4[256];
        while(fgets(b1,sizeof(b1),f1) && fgets(b2,sizeof(b2),f2) && fgets(b4,sizeof(b4),f4))
            {
            char * pb2 = b2;
            bool bambiguous;
            if(*pb2 == '|')
                {
                ++pb2;
                bambiguous = true;
                }
            else
                bambiguous = false;
            long Ref = 0;
            char * hash;
            char * s;
            bool found;
            char * komma = strchr(pb2,',');
            char * tab;
            if(komma && komma[1])
                *komma = '\0';
            else
                trim(pb2);

            trim(b1);
            trim(b4);

            int c = -1;
            bool ambiguousRule = false;
            for(s = b1,found = false;!found && s;)
                {
                ++c;
                tab = NULL;
                hash = strchr(s,'#');
                char * bar;
                if(hash)
                    {
                    *hash = '\0';
                    Ref = strtol(hash+1,NULL,10);
                    if(!strcmp(s,pb2))
                        {
                        *hash = '#';
                        found = true;
                        }
                    else
                        {
                        *hash = '#';
                        bar = strchr(hash,'|');
                        if(bar)
                            s = bar + 1;
                        else
                            s = NULL;
                        }
                    }
                else
                    {
                    bar = strchr(s,'|');
                    if(bar)
                        {
                        ambiguousRule = true;
                        *bar = '\0';
                        }
                    tab = strchr(s,'\t');
                    if(tab)
                        *tab = '\0';
                    if(!strcmp(s,pb2))
                        {
                        if(bar)
                            *bar = '|';
                        if(tab)
                            *tab = '\t';
                        found = true;
                        }
                    else
                        {
                        if(bar)
                            {
                            *bar = '|';
                            s = bar + 1;
                            }
                        else
                            s = NULL;
                        if(tab)
                            *tab = '\t';
                        }
                    }
                }

            if(found)
                {
                if(strchr(b1,'|'))
                    {
                    assert(ambiguousRule);
                    if(c > 2)
                        c = 2;
                    ++ambiguous[c];
                    Correct(Ref);
                    if(f3)
                        fprintf(f3,"|\t%ld\t%s\t%s\t%s\n",Ref,b4,b1,b2);
                    }
                else
                    {
                    assert(!ambiguousRule);
                    ++same;
                    Correct(Ref);
                    if(f3)
                        fprintf(f3,"+\t%ld\t%s\t%s\t%s\n",Ref,b4,b1,b2);
                    }
                }
            else
                {
                ++different;
                Wrong(Ref);
                if(f3)
                    fprintf(f3,"-\t%ld\t%s\t%s\t%s\n",Ref,b4,b1,b2);
                }

            if(ambiguousRule)
                {
                ++ambiguousRules;
                if(bambiguous)
                    Decision.true_amb++;
                else
                    Decision.false_amb++;
                }
            else
                {
                if(bambiguous)
                    Decision.false_not_amb++;
                else
                    Decision.true_not_amb++;
                }

            if(komma)
                *komma = ',';
            }
        }
    else
        {
        if(!f1)
            {
            printf("\n!f1 (output %s) PRESS ENTER\n",output);
            getchar();
            }
        if(!f2)
            {
            printf("\n!f2 (control %s) PRESS ENTER\n",control);
            getchar();
            }
        if(!f3)
            {
            printf("\n!f3 (controlResult %s) PRESS ENTER\n",controlResult);
            getchar();
            }
        if(!f4)
            {
            printf("\n!f4 (test %s) PRESS ENTER\n",test);
            getchar();
            }
        }
    if(f1)
        {
        --openfiles;
        fclose(f1);
        }
    if(f2)
        {
        --openfiles;
        fclose(f2);
        }
    if(f3)
        {
        --openfiles;
        fclose(f3);
        }
    if(f4)
        {
        --openfiles;
        fclose(f4);
        }
    if(REMOVE)
        {
        remove(output);
        remove(controlResult);
        }
    }

static void doTheTest
    (/*const char * executable
    ,*/const char * traintest
    ,const char * test
    ,const char * output
    ,const char * traincontrol
    ,const char * controlResult
    ,const char * control
    ,const char * flexrules
    ,int & same
    ,int & different
    ,ambty & ambiguous
    ,int & ambiguousRules
    ,decision & Decision)
    {
    /*
    char command[1000];

    sprintf(command,"%s -L -l%c -q- -U- -H2 -m0 -t%c -c"
         QUOTE
         "$B\\n"
         QUOTE
         " -B"
         QUOTE
         "$w%s"
         QUOTE
         " -f %s -i %s -o %s"
        ,executable
        ,'-'
        ,TAGGED ? ' ' : '-'
        ,TAGGED ? "/$t" : ""
        ,flexrules
        ,TRAINTEST ? traintest : test
        ,output);

    printf("%s",command);
    printf("\nOpen files:%d\n",::openfiles);
    // C:\utils\cstlemma -L -q- -U- -H2 -m0 -t- -c'$B\n' -B'$w' -f D:\projects\train\flex\el_parms0_154_1.txt0 -i D:\projects\train\test\el_parms0_154_1.txt -o D:\projects\train\sout\el_0_154_1.txt
    system(command);
    */
    lemmatiseFile(TRAINTEST ? traintest : test,flexrules,output);
    if(REMOVE)
        {
        remove(flexrules);
        }
    same = 0;
    different = 0;
    ambiguous[0] = ambiguous[1] = ambiguous[2] = 0;
    ambiguousRules = 0;
    if(TRAINTEST)
        {
        compare(output,traincontrol,same,different,ambiguous,controlResult,traintest,ambiguousRules,Decision);
        }
    else
        {
        compare(output,control,same,different,ambiguous,controlResult,test,ambiguousRules,Decision);
        }
    }

typedef enum {esam,eamb0,eamb1,eamb2,edif,eresceil} eres;
class stddev // standard deviations
    {
    private:
        typedef int tri[eresceil];
        int N;
        int maxN;
        int *x;
        tri *X;
    public:
        void datum(int x)
            {
            assert(N < maxN);
            this->x[N++] = x;
            }
        void datum(int s,ambty a,int d)
            {
            assert(N < maxN);
            this->X[N][esam] = s;
            this->X[N][eamb0] = a[0];
            this->X[N][eamb1] = a[1];
            this->X[N][eamb2] = a[2];
            this->X[N][edif] = d;
            N++;
            }
        stddev()
            {
            assert(globmaxcount > 0);
            maxN = globmaxcount;
            N = 0;
            x = new int[maxN];
            X = new tri[maxN];
            }
        ~stddev()
            {
            delete [] x;
            delete [] X;
            }
        double calculate()
            {
            if(N < 2)
                return 0.0;
            double mean = 0.0;
            int i;
            for(i = 0;i < N;++i)
                mean += x[i];
            mean /= N;
            double sumdeviation2 = 0;
            for(i = 0;i < N;++i)
                {
                double deviation = x[i] - mean;
                sumdeviation2 += deviation*deviation;
                }
            double variance = sumdeviation2/(N-1); // N-1: sample standard deviation. N: standard deviation
            return sqrt(variance);
            }
        double calculate(eres r)
            {
            if(N < 2)
                return 0.0;
            double meanresult = 0.0;
            int i;
            double * result = new double[N];
            for(i = 0;i < N;++i)
                {
                int cases = X[i][esam] + X[i][eamb0] + X[i][eamb1] + X[i][eamb2] + X[i][edif];
                result[i] = (double)X[i][r] / (double)cases;
                meanresult += result[i];
                }
            meanresult /= N;
            double sumdeviation2 = 0;
            for(i = 0;i < N;++i)
                {
                double deviation = (result[i] - meanresult);
                sumdeviation2 += deviation*deviation;
                }
            delete [] result;
            double variance = sumdeviation2/(N-1); // N-1: sample standard deviation. N: standard deviation
            return sqrt(variance);
            }
    };

static void printResults(FILE * fptally,int flexcount,int same,ambty ambiguous,int different,int ambiguousRules)
    {
    int iall = same + ambiguous[0] + ambiguous[1] + ambiguous[2] +  different;
    if(iall > 0)
        {
        double all = 0.01*iall;
        assert(all > 0);
        fprintf(fptally,"#rules:%d #same:%d #ambiguous:%d %d %d #different:%d %14.6f %14.6f %14.6f %14.6f %14.6f (amb %d all %d amb%% %14.6f)\n"
            ,flexcount,same,ambiguous[0],ambiguous[1],ambiguous[2],different,(double)same/all,(double)ambiguous[0]/all,(double)ambiguous[1]/all,(double)ambiguous[2]/all,(double)different/all,ambiguousRules,iall,(double)ambiguousRules/all);
        }
    else
        fprintf(fptally,"same + ambiguous[0] + ambiguous[1] + ambiguous[2] +  different = %d\n"
            ,same + ambiguous[0] + ambiguous[1] + ambiguous[2] +  different);
    fflush(fptally);
    }

typedef struct countStruct
    {
    int tflexcount;
    int tsame;
    int tdifferent;
    ambty tambiguous;
    int same;
    ambty ambiguous;
    int different;
    int ambiguousRules;
    int tambiguousRules;
    decision Decision;
    countStruct()
        {
        tflexcount = tsame = tdifferent = same = different = ambiguousRules = tambiguousRules = 0;
        for(int i = 0;i < 3;++i)
            {
            tambiguous[i] = 0;
            ambiguous[i] = 0;
            }
        }
    } countStruct;

class counting
    {
    private:
        stddev StandardDev;
    public:
        countStruct n;
        counting()
            {            
            }
        void testing
            ( int cutoff
            , const char * Affixrules
            , const char * output
            //, const char * executable
            , const char * traintest
            , const char * test
            , const char * traincontrol
            , const char * controlResultN
            , const char * control
            , FILE * fptally
            , optionStruct * Options
            )
            {
            long nflexcount = 0;
            char textualRules[200];
            sprintf(textualRules,"%snumberOfRules_%s.pass%d_%d.txt",Options->tempDir(),LGf(Options),Options->redo() ? 2 : 1,cutoff);
            FILE * ft = fopen(textualRules,"rb");
            ++openfiles;
            if(ft)
                {
                char Line[100];
                if(fgets(Line,sizeof(Line),ft) != 0)
                    {
                    nflexcount = strtol(Line,NULL,10);
                    n.tflexcount += nflexcount;
                    }
                --openfiles;
                fclose(ft);
                }
            else
                {
                fprintf(stderr,"Cannot open \"%s\" for reading\n",textualRules);
                }

            doTheTest
                (/*executable
                ,*/traintest
                ,test
                ,output
                ,traincontrol
                ,controlResultN
                ,control
                ,Affixrules
                ,n.same
                ,n.different
                ,n.ambiguous
                ,n.ambiguousRules
                ,n.Decision
                );
            this->StandardDev.datum(n.same,n.ambiguous,n.different);
            /*printf("n.same %d + n.different %d + n.ambiguous %d %d %d = nsum %d\n"
                ,n.same,n.different,n.ambiguous[0],n.ambiguous[1],n.ambiguous[2],n.same+n.different+n.ambiguous[0]+n.ambiguous[1]+n.ambiguous[2]);*/
            n.tsame += n.same;
            n.tdifferent += n.different;
            for(int j = 0;j < 3;++j)
                {
                n.tambiguous[j] += n.ambiguous[j];
                }
            n.tambiguousRules += n.ambiguousRules;
            fprintf(fptally,"cutoff:%d ",cutoff);
            printResults(fptally,nflexcount,n.same,n.ambiguous,n.different,n.ambiguousRules);
            fflush(fptally);
            if(counts)
                {
                delete counts;
                counts = NULL;
                }
            }
        void printing(int maxcount,FILE * fptab,int ttrainlines)
            {
            double ntot = n.tsame + n.tambiguous[0] + n.tambiguous[1] + n.tambiguous[2] + n.tdifferent;

            fprintf(fptab,  "%10d %14.6f %14.6f %14.6f "
                ,maxcount
                ,(double)ttrainlines/(double)maxcount
                ,(double)n.tflexcount/(double)maxcount
                ,ttrainlines == 0 ? 100.0 : 100.0*(double)n.tflexcount/(double)ttrainlines
                );
            fprintf(fptab,  "%14.6f %14.6f %14.6f %14.6f %14.6f "
                ,(double)n.tsame/(double)maxcount
                ,(double)n.tambiguous[0]/(double)maxcount
                ,(double)n.tambiguous[1]/(double)maxcount
                ,(double)n.tambiguous[2]/(double)maxcount
                ,(double)n.tdifferent/(double)maxcount
                );
            fprintf(fptab,  "%14.6f %14.6f %14.6f %14.6f %14.6f "
                ,ntot > 0 ? 100.0*this->StandardDev.calculate(esam) : 0
                ,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb0) : 0
                ,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb1) : 0
                ,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb2) : 0
                ,ntot > 0 ? 100.0*this->StandardDev.calculate(edif) : 0
                );
            fprintf(fptab,  "%14.6f %14.6f %14.6f %14.6f %14.6f %14.6f "
                ,ntot > 0 ? 100.0*(double)n.tsame/ntot : 0
                ,ntot > 0 ? 100.0*(double)n.tambiguous[0]/ntot : 0
                ,ntot > 0 ? 100.0*(double)n.tambiguous[1]/ntot : 0
                ,ntot > 0 ? 100.0*(double)n.tambiguous[2]/ntot : 0
                ,ntot > 0 ? 100.0*(double)n.tdifferent/ntot : 0
                ,ntot > 0 ? 100.0*(double)n.tambiguousRules/ntot : 0
                );
            fprintf(fptab,  "%14.6f %14.6f %14.6f %14.6f "
                ,ntot > 0 ? 100.0*(double)this->n.Decision.false_amb/ntot : 0
                ,ntot > 0 ? 100.0*(double)this->n.Decision.false_not_amb/ntot : 0
                ,ntot > 0 ? 100.0*(double)this->n.Decision.true_amb/ntot : 0
                ,ntot > 0 ? 100.0*(double)this->n.Decision.true_not_amb/ntot : 0
                );
            fprintf(fptab,  "%14.6f %14.6f\n"
                ,this->n.Decision.precision()
                ,this->n.Decision.recall()
                );
            fflush(fptab);
            }
        void validationcolumn(char ** cell,int cutoff,int maxcount,int ttrainlines)
            {
            double ntot = n.tsame + n.tambiguous[0] + n.tambiguous[1] + n.tambiguous[2] + n.tdifferent;
            char * f1 = "%14d ";
            char * f2 = "%14.6f ";
            for(int i = 0;i < 20;++i)
                cell[i] = new char[20];
            sprintf(cell[0],f1,cutoff);
            sprintf(cell[1],f2,(double)n.tflexcount/(double)maxcount);
            sprintf(cell[2],f2,ttrainlines == 0 ? 100.0 : 100.0*(double)n.tflexcount/(double)ttrainlines);
            sprintf(cell[3],f2,ntot > 0 ? 100.0*(double)n.tsame/ntot : 0.0);
            sprintf(cell[4],f2,ntot > 0 ? 100.0*(double)n.tambiguous[0]/ntot : 0.0);
            sprintf(cell[5],f2,ntot > 0 ? 100.0*(double)n.tambiguous[1]/ntot : 0.0);
            sprintf(cell[6],f2,ntot > 0 ? 100.0*(double)n.tambiguous[2]/ntot : 0.0);
            sprintf(cell[7],f2,ntot > 0 ? 100.0*(double)n.tdifferent/ntot : 0.0);
            sprintf(cell[8],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(esam) : 0.0);
            sprintf(cell[9],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb0) : 0.0);
            sprintf(cell[10],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb1) : 0.0);
            sprintf(cell[11],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb2) : 0.0);
            sprintf(cell[12],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(edif) : 0.0);
            sprintf(cell[13],f2,ntot > 0 ? 100.0*(double)n.tambiguousRules/ntot : 0.0);
            sprintf(cell[14],f2,ntot > 0 ? 100.0*(double)this->n.Decision.false_amb/ntot : 0.0);
            sprintf(cell[15],f2,ntot > 0 ? 100.0*(double)this->n.Decision.false_not_amb/ntot : 0.0);
            sprintf(cell[16],f2,ntot > 0 ? 100.0*(double)this->n.Decision.true_amb/ntot : 0.0);
            sprintf(cell[17],f2,ntot > 0 ? 100.0*(double)this->n.Decision.true_not_amb/ntot : 0.0);
            sprintf(cell[18],f2,this->n.Decision.precision());
            sprintf(cell[19],f2,this->n.Decision.recall());
            }
        void printingNice(int cutoff,double fraction,int maxcount,FILE * fptab,int ttrainlines,optionStruct * Options)
            {
            double ntot = n.tsame + n.tambiguous[0] + n.tambiguous[1] + n.tambiguous[2] + n.tdifferent;
            fprintf(fptab,"\nNew algorithm, least wrongly lemmatised (MIN(diff)).\n");
            fprintf(fptab,"Suffix only     %s\n",Options->suffixOnly() ? "yes" : "no");
            fprintf(fptab,"Redo training   %s\n",Options->redo() ? "yes" : "no");
            fprintf(fptab,"cutoff          %d\n",cutoff);
            fprintf(fptab,"fraction        %14.6f\n",fraction);
            fprintf(fptab,"iterations      %14d\n",maxcount);
            fprintf(fptab,"trainlines     %14.6f\n",(double)ttrainlines/(double)maxcount);
            fprintf(fptab,"rules          %14.6f\n",(double)n.tflexcount/(double)maxcount);
            fprintf(fptab,"rules%%         %14.6f\n",ttrainlines == 0 ? 100.0 : 100.0*(double)n.tflexcount/(double)ttrainlines);
            fprintf(fptab,"same%%stdev     %14.6f\n",ntot > 0 ? 100.0*this->StandardDev.calculate(esam) : 0.0);
            fprintf(fptab,"ambi1%%stdev    %14.6f\n",ntot > 0 ? 100.0*this->StandardDev.calculate(eamb0) : 0.0);
            fprintf(fptab,"ambi2%%stdev    %14.6f\n",ntot > 0 ? 100.0*this->StandardDev.calculate(eamb1) : 0.0);
            fprintf(fptab,"ambi3%%stdev    %14.6f\n",ntot > 0 ? 100.0*this->StandardDev.calculate(eamb2) : 0.0);
            fprintf(fptab,"diff%%stdev     %14.6f\n",ntot > 0 ? 100.0*this->StandardDev.calculate(edif) : 0.0);
            fprintf(fptab,"same%%          %14.6f\n",ntot > 0 ? 100.0*(double)n.tsame/ntot : 0.0);
            fprintf(fptab,"ambi1%%         %14.6f\n",ntot > 0 ? 100.0*(double)n.tambiguous[0]/ntot : 0.0);
            fprintf(fptab,"ambi2%%         %14.6f\n",ntot > 0 ? 100.0*(double)n.tambiguous[1]/ntot : 0.0);
            fprintf(fptab,"ambi3%%         %14.6f\n",ntot > 0 ? 100.0*(double)n.tambiguous[2]/ntot : 0.0);
            fprintf(fptab,"diff%%          %14.6f\n",ntot > 0 ? 100.0*(double)n.tdifferent/ntot : 0.0);
            fprintf(fptab,"amb.rules%%     %14.6f\n",ntot > 0 ? 100.0*(double)n.tambiguousRules/ntot : 0.0);
            fprintf(fptab,"false_amb%%     %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Decision.false_amb/ntot : 0.0);
            fprintf(fptab,"false_not_amb%% %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Decision.false_not_amb/ntot : 0.0);
            fprintf(fptab,"true_amb%%      %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Decision.true_amb/ntot : 0.0);
            fprintf(fptab,"true_not_amb%%  %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Decision.true_not_amb/ntot : 0.0);
            fprintf(fptab,"precision      %14.6f\n",this->n.Decision.precision());
            fprintf(fptab,"recall         %14.6f\n",this->n.Decision.recall());
            fflush(fptab);
            }
    };


void trainAndTest
        (int linecnt
#if CLUMPSENSITIVE
        ,int clumpcnt
#endif
        ,const char * XT
        ,const char * TT
        ,line * lines
#if CLUMPSENSITIVE
        ,clump * clumps
#endif
         ,optionStruct * Options
        )
    {
    lineab AffixLine[CUTOFFS];
    lineab SuffixLine[CUTOFFS];
    char formatprefix[256]          ;sprintf(formatprefix,          BASEDIR "%s%%s" SEPARATOR "%s_%s%s%s%s%s_%%s.txt",Options->tempDir(),LGf(Options),XTRf(Options),(Options->suffixOnly() ? "_suffix" : "_affix"),XT,TT,(Options->redo() ? "redone" : "singleshot"));

    char formatTraining[256]        ;sprintf(formatTraining,        formatprefix,"training"     ,"%d_%d");        // the training words
    char formatTest[256]            ;sprintf(formatTest,            formatprefix,"test"         ,"%d_%d");        // the test words (<> training words)
    char formatTrainTest[256]       ;sprintf(formatTrainTest,       formatprefix,"test"         ,"_train_%d_%d"); // the test words (= training words)
    char formatControl[256]         ;sprintf(formatControl,         formatprefix,"control"      ,"%d_%d");        // the correct answers (<> training words)
    char formatTrainControl[256]    ;sprintf(formatTrainControl,    formatprefix,"control"      ,"_train_%d_%d"); // the correct answers (= training words)

    char formatFlexrules[256]       ;sprintf(formatFlexrules,       formatprefix,"flex"         ,"%d_%d");
    char formatSOutput[256]         ;sprintf(formatSOutput,         formatprefix,"sout"         ,"%d_%d_%d");
    char formatOutput[256]          ;sprintf(formatOutput,          formatprefix,"out"          ,"%d_%d_%d");
    char formatControlResult[256]   ;sprintf(formatControlResult,   formatprefix,"resultSuffix" ,"%d_%d_cutoff_%d");

    char formatControlResultN[256]  ;sprintf(formatControlResultN,  formatprefix,"resultAffix"  ,"%d_%d_cutoff_%d");

    char formatLemmas[256]          ;sprintf(formatLemmas,          formatprefix,"lemmas"       ,"");
    char formatWeird[256]           ;sprintf(formatWeird,           formatprefix,"weird"        ,"");

    char formatTally[256]           ;sprintf(formatTally,           formatprefix,"tally"        ,"");
    char formatTab[256]             ;sprintf(formatTab,             formatprefix,"tab"          ,"");
    char tally[256];
    char tab[256];
    sprintf(tally,"%s",formatTally);
    sprintf(tab,"%s",formatTab);
    FILE * fptally;    
    FILE * fptab;    
    fptally = fopen(tally,"wb");
    ++openfiles;
    if(fptally)
        {
        fprintf(fptally,"Lemmatiser training+test\n");
        --openfiles;
        fclose(fptally);
        fptally = 0;
        }
    fptab = fopen(tab,"wb");
    ++openfiles;
    if(fptab)
        {
        fprintf(fptab,"Lemmatiser training+test\n");
        --openfiles;
        fclose(fptab);
        fptab = 0;
        }

    int fraction = 1;// 0 <= fraction <= 10000
    for(fraction = FRACTION_LOW;fraction <= FRACTION_HIGH;)
        {
        if(!fptally)
            {
            fptally = fopen(tally,"ab");
            ++openfiles;
            }
        if(fptally)
            {
            fprintf(fptally,"\nfraction:%d/10000\n",fraction);
            --openfiles;
            fclose(fptally);
            fptally = 0;
            }
        if(!fptab)
            {
            fptab = fopen(tab,"ab");
            ++openfiles;
            }
        if(fptab)
            {
            fprintf(fptab,"cutoff  fraction  iterations    trainlines     affixrules         affix%%      a-same        a-ambiguous                                 a-different    a-same-stddev%% a-amb-stddev%%                               a-diff-stddev%%   a-same%%       a-ambiguous%%                                a-different%%   a-amb.rules%%    a_false_amb a_false_not_amb    a_true_amb  a_true_not_amb    a_precision     a_recall\n");
            --openfiles;
            fclose(fptab);
            fptab = 0;
            }
        int ttrainlines = 0;
        int maxcount;
#if N_FOLDCROSSVALIDATION
        maxcount = N_;
#else
        maxcount = (int)(5000000.0*sqrt(1.0/((double)fraction*(10000.0-fraction)*(double)linecnt)));
        if(maxcount < 1)
            maxcount = 1;
        else if(maxcount > 30)
            maxcount = 30;
#endif
        char paramfile[100];
        sprintf(paramfile,"%s%s%s%s%s.txt",LGf(Options),(Options->suffixOnly() ? "_suffix" : "_affix"),XTRf(Options),TT,(Options->redo() ? "redone" : "singleshot"));
        globmaxcount = maxcount;
        counting c[CUTOFFS]; // indexed by cutoff
        int count = 1;
        for( count = 1
            ; count <= maxcount
            ; ++count
            )
            {
            char test[256];
            char Tcontrol[256];
            char Ttraincontrol[256];
            char traintest[256];
            char controlResult[256];
            char controlResultN[256];
            char Soutput[256];
            char output[256];

            sprintf(Lemmas,"%s",formatLemmas);
            sprintf(Weird,"%s",formatWeird);
            if(!fptally)
                {
                fptally = fopen(tally,"ab");
                ++openfiles;
                }
            if(fptally)
                {
                fprintf(fptally,"\ncount:%d/10000\n",count);
                --openfiles;
                fclose(fptally);
                fptally = 0;
                }
            if(!counts)
                counts = new Count(0); // for those words that didn't match any rule.

            char Ttraining[256];
            sprintf(Ttraining,formatTraining,fraction,count);
            sprintf(test,formatTest,fraction,count);
            sprintf(Tcontrol,formatControl,fraction,count);
            if(TRAINTEST )
                {
                sprintf(traintest,formatTrainTest,fraction,count);
                sprintf(Ttraincontrol,formatTrainControl,fraction,count);
                }
            static char suffixrulesArr[256];
            sprintf(suffixrulesArr,formatFlexrules,fraction,count);

            int trainlines = 0;
#if CLUMPSENSITIVE
            if(clumpcnt > 1)
                trainlines = splitLemmaliste(clumpcnt,training,test,control,1,2,3,fraction,traintest,traincontrol,clumps,lines);
            else
#endif
                trainlines = splitLemmaliste(linecnt,Ttraining,test,Tcontrol,1,2,3,fraction,traintest,Ttraincontrol,lines);
            if(trainlines > 0)
                {
                ttrainlines += trainlines;

                optionStruct testOptions(*Options);
                testOptions.seti(Ttraining);
                testOptions.setc(CUTOFF_HIGH);
                testOptions.seto(newStyleRulesf(Options));
                testOptions.sete(LGf(Options));
                testOptions.setn("123");
                testOptions.setf(XTRf(Options));
                testOptions.setP(paramfile);

                testOptions.completeArgs();
                trainRules("",&testOptions);
                if(REMOVE)
                    {
                    remove(Ttraining);
                    //printf("removed training files %s\n",Ttraining);
                    }
                }

            if(ttrainlines > 0)
                {
                for(int cutoff = CUTOFF_LOW;cutoff <= CUTOFF_HIGH;++cutoff)
                    {
                    //printf("cutoff %d <= %d <= %d\n",CUTOFF_LOW,cutoff,CUTOFF_HIGH);
                    sprintf(controlResult,formatControlResult,fraction,count,cutoff);
                    sprintf(controlResultN,formatControlResultN,fraction,count,cutoff);
                    char Affixrules[250];
                    const char * lastslash = strrchr(newStyleRulesf(Options),*SLASH);
                    const char * filename;
                    if(lastslash)
                        {
                        filename = lastslash + 1;
                        sprintf(Affixrules,"%.*s%d%s%s",(int)(filename - newStyleRulesf(Options)),newStyleRulesf(Options),cutoff,SLASH,filename);
                        }
                    else
                        sprintf(Affixrules,"%d%c%s",cutoff,DIRSEP,newStyleRulesf(Options));

                    sprintf(Soutput,formatSOutput,cutoff,fraction,count);
                    sprintf(output,formatOutput,cutoff,fraction,count);
                    //printf("Affixrules:%s\n",Affixrules);
                    //const char executable[] = "cstlemma";
                    FILE * fptally = fopen(tally,"ab");
                    ++openfiles;
                    if(fptally)
                        {
                        c[cutoff].testing(cutoff,Affixrules,output/*,executable*/,traintest,
                            test,Ttraincontrol,controlResultN,Tcontrol,fptally,Options);
                        --openfiles;
                        fclose(fptally);
                        }
                    else
                        {
                        fprintf(stderr,"Cannot open \"%s\" for appending. No testing is done.\n",tally);
                        }
#if OnCE
                    break;
#endif
                    }
                }
            if(REMOVE)
                {
                remove(test);
                if(TRAINTEST)
                    {
                    remove(traintest);
                    remove(Ttraincontrol);
                    }
                remove(Tcontrol);
                }
#if OnCE
            break; // STOP after after first results
#endif
            }
        if(ttrainlines > 0)
            {
            fptab = fopen(tab,"ab");
            ++openfiles;
            if(fptab)
                {
                for(int cutoff = CUTOFF_LOW;cutoff <= CUTOFF_HIGH;++cutoff)
                    {
                    fprintf(fptab,"%d %14.6f ",cutoff,(double)fraction/10000.0);
                    c[cutoff].printing(maxcount,fptab,ttrainlines);
                    }
                for(int cutoff = CUTOFF_LOW;cutoff <= CUTOFF_HIGH;++cutoff)
                    {
                    AffixLine [cutoff].add((double)c[cutoff].n.tflexcount/(double)maxcount,(double)ttrainlines/(double)maxcount);
                    fprintf (fptab,"cutoff %d Affix  a %14.6lf b %14.6lf: N(rules)=%14.6lf*N(trainpairs)^%lf\n"
                        ,cutoff,AffixLine [cutoff].a(),AffixLine [cutoff].b()
                        ,exp(AffixLine[cutoff].a())
                        ,AffixLine[cutoff].b()
                        );
                    }
                fflush(fptab);
                int lowestntdifferent = -1;
                int bestcutoff = 0;
                for(int cutoff = CUTOFF_LOW;cutoff <= CUTOFF_HIGH;++cutoff)
                    {
                    if(lowestntdifferent == -1 || lowestntdifferent > c[cutoff].n.tdifferent)
                        {
                        bestcutoff = cutoff;
                        lowestntdifferent = c[cutoff].n.tdifferent;
                        }
                    }
                char * texts[] =    {
                                    "cutoff         ",
                                    "rules          ",
                                    "rules%         ",
                                    "same%          ",
                                    "ambi1%         ",
                                    "ambi2%         ",
                                    "ambi3%         ",
                                    "diff%          ",
                                    "same%stdev     ",
                                    "ambi1%stdev    ",
                                    "ambi2%stdev    ",
                                    "ambi3%stdev    ",
                                    "diff%stdev     ",
                                    "amb.rules%     ",
                                    "false_amb%     ",
                                    "false_not_amb% ",
                                    "true_amb%      ",
                                    "true_not_amb%  ",
                                    "precision      ",
                                    "recall         "
                                    };
                int nocol = 7;
                int norow = sizeof(texts)/sizeof(texts[0]);
                int nocell = norow*nocol;
                char * cell[sizeof(texts)/sizeof(texts[0]) * 7];
                for(int j = 0;j < norow;++j)
                    cell[j] = texts[j];
                for(int k = 0;k < nocol - 1;++k)
                    {
                    c[k].validationcolumn(cell+(k+1)*norow,k,maxcount,ttrainlines);
                    }
                int len = 0;
                for(int m = 0;m < nocell;++m)
                    len += strlen(cell[m]);
                char * report = new char[len + 3*norow + 1];
                report[0] = 0;
                for(int row = 0;row < norow;++row)
                    {
                    strcat(report,"; ");
                    for(int col = 0;col < 7;++col)
                        {
                        int ind = row + col*norow;
                        strcat(report,cell[ind]);
                        }
                    strcat(report,"\n");
                    }
                Options->printArgFile(report);
                delete [] report;
                for(int row2 = 0;row2 < norow;++row2)
                    {
                    for(int col = 1;col < 7;++col)
                        {
                        delete [] cell[row2+col*norow];
                        }
                    }
                c[bestcutoff].printingNice(bestcutoff,fraction,maxcount,fptab,ttrainlines,Options);
                FILE * fp = fopen(paramfile,"r");
                ++openfiles;
                if(fp)
                    {
                    fputc('\n',fptab);
                    int K;
                    while((K = fgetc(fp)) != EOF)
                        fputc(K,fptab);
                    --openfiles;
                    fclose(fp);
                    }
                --openfiles;
                fclose(fptab);
                fptab = 0;
                }
            }
        switch(STEPSIZE)
            {
            case 100:
                if(fraction < 1)
                    ++fraction;
                else
                    fraction *= 2; 
                break;
            case 50:
                if(fraction < 2)
                    ++fraction;
                else 
                    {
                    fraction *= 3; 
                    fraction /= 2;
                    }
                break;
            default:
            case 25:
                if(fraction < 4)
                    ++fraction;
                else 
                    {
                    fraction *= 5; 
                    fraction /= 4;
                    }
                break;
            case 10:
                if(fraction < 10)
                    ++fraction;
                else 
                    {
                    fraction *= 11; 
                    fraction /= 10;
                    }
                break;
            case 5:
                if(fraction < 20)
                    ++fraction;
                else 
                    {
                    fraction *= 21; 
                    fraction /= 20;
                    }
                break;
            }
#if OnCE
        break; // STOP after after first results
#endif
        }
    }


static int readFile
        (line *& lines
#if CLUMPSENSITIVE
        ,clump *& clumps
        ,int & clumpcnt
#endif
        ,int sep
        ,optionStruct * Options
        )
    {
    return readlines(1,2,3,lines,sep
#if CLUMPSENSITIVE
        ,clumps,clumpcnt
#endif
        ,Options
        );
    }


static int doit(optionStruct * Options)
    {
    const char * XT;
    const char * TT;
    init(Options->trainTest(),XT,TT);
    int sep = '\t';

    if(XTRf(Options) == NULL)
        {
        printf("Rule comparison function (variable XTRf(Options)) not specified for language %s\n",LGf(Options));
        exit(1);
        }
    
    line * lines = NULL;
#if CLUMPSENSITIVE
    clump * clumps = NULL;
    int clumpcnt = 0;
#endif
    int linecnt = readFile
        (lines
#if CLUMPSENSITIVE
        ,clumps
#endif
#if CLUMPSENSITIVE
        ,clumpcnt
#endif
        ,sep
        ,Options
        );
    if(linecnt)
        {
        trainAndTest
            (linecnt
#if CLUMPSENSITIVE
            ,clumpcnt
#endif
            ,XT
            ,TT
            ,lines
#if CLUMPSENSITIVE
            ,clumps
#endif
            ,Options
            );






        delete [] lines;
        lines = NULL;
        }
    return 0;
    }

int testrules(optionStruct * Options)
    {
    if(lemmalistef(Options) && LGf(Options) && XTRf(Options))
        {
        doit(Options);
        return 0;
        }

    return 0;
    }
