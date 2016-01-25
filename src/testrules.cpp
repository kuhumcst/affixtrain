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
#define ENSURELEMMAPRESENTASFULLFORM 0 // If set to 1 creates doublets! Set to 0 and run again on output from first run to remove doublets.
#define FORCESAMECASING 0
#if FORCESAMECASING
#include "utf8func.h"
#endif

static bool TAGGED;

static int FRACTION_LOW; // min 0
static int FRACTION_HIGH; // max 10000
static int CUTOFF_LOW;
static int CUTOFF_HIGH;
static int MAXCOUNT;
static const int STEPSIZE = 100;

static int globmaxcount = 0;

static const char * XTRf(optionStruct * Options)
    {
    static char xtr[10];
    switch(Options->getWeightFunction())
        {
        case econstant:
            sprintf(xtr,"%s%sXC"
                ,Options->redo() ? "R" : ""
                ,Options->suffixOnly() ? "S" : ""
                );
            break;
        case esupport:
            sprintf(xtr,"C%d%s%sXW"
                ,Options->expectedCutoff()
                ,Options->redo() ? "R" : ""
                ,Options->suffixOnly() ? "S" : ""
                );
            break;
        case eentropy:
            sprintf(xtr, "%s%sXE"
                    , Options->redo() ? "R" : ""
                    , Options->suffixOnly() ? "S" : ""
                    );
            break;
        case edepth:
            sprintf(xtr,"%s%sXD"
                ,Options->redo() ? "R" : ""
                ,Options->suffixOnly() ? "S" : ""
                );
            break;
        case esize:
            sprintf(xtr,"%s%sXS"
                ,Options->redo() ? "R" : ""
                ,Options->suffixOnly() ? "S" : ""
                );
            break;
        default:
            *xtr = 0;
        }
    return xtr;
    }

static const char * LGf(optionStruct * Options)
    {
    return Options->extra();
    }

static const char * lemmalistef(optionStruct * Options)
    {
    return Options->wordLemmaList();
    }

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
        bool computable()
            {
            return N >= 2;
            }
    };

#define CUTOFFS 10

static char Lemmas[256];
static char Weird[256];

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
    ~line(){delete [] s;}
    };

static int mystrcmp(const void * a,const void * b)
    {
    return strcmp(((line *)a)->s, ((line * )b)->s);
    }

struct clump
    {
    line * start;
    int linecnt;
    };

static void countLinesAndClumps(FILE * fpi,int & linecnt,int & clumpcnt)
    {
    linecnt = 0;
    clumpcnt = 0;
    int kar;
    int oldkar = '\n';
    bool nonEmptyLineSeen = false;
    bool separatorLineSeen = false;
    while((kar = fgetc(fpi)) != EOF)
        {
        if(kar == '\r')
            continue;
        if(kar == '\n')
            {
            if(oldkar == '\n')
                {
                if(!separatorLineSeen && nonEmptyLineSeen)
                    { // first empty line after clump found
                    nonEmptyLineSeen = false;
                    separatorLineSeen = true;
                    }
                }
            else
                {
                if(!nonEmptyLineSeen)
                    {
                    ++clumpcnt;
                    nonEmptyLineSeen = true;
                    separatorLineSeen = false;
                    }
                ++linecnt;
                }
            }
        oldkar = kar;
        }
    }

static int fileRead(line * lines,
                    clump clumps[],
                    FILE * fpi,int columnfull,int columnbase,int columnPOS,int sep)
    {
    char buf[8192];
    char * pbuf = buf;
    int linecnt = 0;
    int clumpcnt = 0;
    int kar;
    int oldkar = 0;
    bool nonEmptyLineSeen = false;
    bool separatorLineSeen = false;
    while((kar = fgetc(fpi)) != EOF)
        {
        if(kar == '\r')
            continue; // DOS file. Ignore CR. (Under Windows, CR is ignored per default.)
        if(kar == '\n')
            {
            if(oldkar == kar)
                {
                if(!separatorLineSeen)
                    { // clump found
                    if(nonEmptyLineSeen)
                        {
                        assert(clumpcnt > 0);
                        nonEmptyLineSeen = false;
                        separatorLineSeen = true;
                        clumps[clumpcnt-1].linecnt = lines+linecnt - clumps[clumpcnt-1].start;
                        }
                    }
                }
            else
                {
                *pbuf++ = '\0';
                pbuf = buf;
                trim(buf);
                if(*buf)
                    {
                    if(!nonEmptyLineSeen)
                        {
                        clumps[clumpcnt].start = lines + linecnt;
                        ++clumpcnt;
                        nonEmptyLineSeen = true;
                        separatorLineSeen = false;
                        }
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
                else
                    {
                    printf("Line %d contains only whitespace.\n",linecnt+clumpcnt);
                    return -1;
                    }
                }
            }
        else
            {
            if((size_t)(pbuf - buf + 2) > sizeof(buf))
                {
                *pbuf = 0;
                printf("Buffer buf too small for line starting with [%s] in function fileRead in file testrules.c\n",buf);
                exit(-1);
                }
            *pbuf++ = (char)kar;
            }
        oldkar = kar;
        }
    clumps[clumpcnt-1].linecnt = lines+linecnt - clumps[clumpcnt-1].start;

    return linecnt;
    }

static int removeDuplicateLines(clump * Clump)
    {
    qsort(Clump->start,Clump->linecnt,sizeof(line),mystrcmp);
    int i;
    int j;
    for(i = 0,j = i + 1;j < Clump->linecnt;++j)
        {
        if(strcmp(Clump->start[i].s,Clump->start[j].s))
            {
            ++i;
            if(i != j)
                {
                assert(Clump->start[i].s == NULL);
                Clump->start[i] = Clump->start[j];
                Clump->start[j].s = 0;
                }
            }
        else
            {
            delete Clump->start[j].s;
            Clump->start[j].s = 0;
            }
        }
    Clump->linecnt = i+1;
    for(int j = 0;j < Clump->linecnt - 1;++j)
        {
        line * A = Clump->start+j;
        line * B = Clump->start+j+1;
        char * as = A->s;
        char * bs = B->s;
        while(*as && *as == *bs && *as != '\t')
            {
            ++as;
            ++bs;
            }
        if(*as == '\t' && *bs == '\t')
            {
            ++as;
            ++bs;
            while(*as && (*as != '\t') && (*as != '\r') && (*as != '\n') && (*as == *bs))
                {
                ++as;
                ++bs;
                }
            if(*as != *bs)
                {
                A->ambiguous = true;
                B->ambiguous = true;
                }
            }
        }
    return Clump->linecnt;
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

static void randomix(clump * clumps,int clumpcnt)
    {
    int i;
    for(i = 0;i < clumpcnt;++i)
        {
        int j = (int)((double)rand() * (double)clumpcnt / (double)RAND_MAX);
        if(j != i && j > 0 && j < clumpcnt)
            {
            clump tmp = clumps[i];
            clumps[i] = clumps[j];
            clumps[j] = tmp;
            }
        }
    }

static int readlines(int columnfull,int columnbase,int columnPOS,line *& lines,int sep
                     ,clump *& clumps,int & clumpcnt
                     ,optionStruct * Options
                     )
    {
    FILE * fpi = NULL;
    int Linecnt = 0;
    //fpi = fopen(lemmalistef(Options),"r");
    fpi = fopen(lemmalistef(Options),"rb");
    ++openfiles;
    if(!fpi)
        {
        fprintf(stderr,"Cannot open \"%s\" for reading\n",lemmalistef(Options));
        return 0;
        }

    srand(1);
    if (Options->verbose())
        printf("readlines countLinesAndClumps\n");
    countLinesAndClumps(fpi,Linecnt,clumpcnt);
    if (Options->verbose())
        printf("readlines countLinesAndClumps: %d lines\n",Linecnt);
    assert(!lines);
    rewind(fpi);
    lines = new line[Linecnt];
    if(clumpcnt)
        {
        if (Options->verbose())
            printf("readlines countLinesAndClumps: %d clumps\n",clumpcnt);
        clumps = new clump[clumpcnt];
        clumps[0].start = 0;
        }
    else
        clumps = NULL;

    if (Options->verbose())
        printf("readlines fileRead\n");
    Linecnt = fileRead(lines,clumps,fpi,columnfull,columnbase,columnPOS,sep);
    if(Linecnt < 0)
        {
        printf("Error in file \"%s\"\n",lemmalistef(Options));
        exit(-1);
        }
    --openfiles;
    fclose(fpi);
    if(clumpcnt < 2)
        {
        clump Clump = {lines,Linecnt};
        if (Options->verbose())
            printf("readlines removeDuplicateLines\n");
        removeDuplicateLines(&Clump);
        if (Options->verbose())
            printf("readlines randomix\n");
        randomix(Clump.start,Clump.linecnt);
        lines = Clump.start;
        Linecnt = Clump.linecnt;
        }
    else
        {
        if (Options->verbose())
            printf("readlines removeDuplicateLines in %d clumps\n",clumpcnt);
        for(int m = 0; m < clumpcnt;++m)
            removeDuplicateLines(clumps+m);

        if (Options->verbose())
            printf("readlines randomix\n");
        randomix(clumps,clumpcnt);
        }
    if (Options->verbose())
        printf("readlines DONE\n");
    return Linecnt;
    }


static bool IsSpace(int k)
    {
    return k == ' ' || k == '\t' || k == '\r' || k == '\n' || !k;
    }

static void writeTestAndControl(line * Line,FILE * ftest,FILE * fcontrol,int columntest, int columncontrol, int columntag)
    {
    bool slashprinted = false;
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
            }
        else
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


static void printTrain(line * Line,FILE * fptrain,FILE * fptraintest,FILE * fptraincontrol,int columntest,int columncontrol,int columntag,bool TrainTest)
    {
    if(!IsSpace(Line->s[0]))
        {
        fprintf(fptrain,"%s\n",Line->s);
        if(TrainTest)
            {
            writeTestAndControl(Line,fptraintest,fptraincontrol,columntest,columncontrol,columntag);
            }
        }
    else
        {
        printf("%s\n",Line->s);
        }
    }

static int TenthCnt = -1;
static int Tenth = -1;

static void tenth(int t)
    {
    Tenth = t;
    Tenth = Tenth % 10;
    TenthCnt = -1;
    }

static int aTenth(void)
    {
    ++TenthCnt;
    if((TenthCnt+Tenth) % 10 == 0)
        return RAND_MAX;
    else
        return 0;
    }

static int splitLemmaList
   (int clumpcnt
   ,char * training
   ,char * test
   ,char * control
   ,int columntest
   ,int columncontrol
   ,int columntag
   ,int fraction
   ,char * traintest
   ,char * traincontrol
   ,clump * clumps
   ,bool TrainTest
   ,bool TenFoldXValidation
   ) // 0 <= fraction <= 10000
    {
    FILE * fptrain = NULL;
    FILE * fptest = NULL;
    FILE * fpcontrol = NULL;
    FILE * fptraintest = NULL;
    FILE * fptraincontrol = NULL;
    fptrain = fopen(training,"wb");
    ++openfiles;
    int trainlines = 0;
    int trainclumps = 0;
    int testclumps = 0;
    int (*Rand)(void);
    if(TenFoldXValidation)
        {
        fraction = 9000;
        Rand = aTenth;
        }
    else
        Rand = rand;
    trainclumps = (int)((double)clumpcnt * (double)fraction / 10000.0);

    testclumps = clumpcnt - trainclumps;
    if(TrainTest)
        {
        fptraintest = fopen(traintest,"wb");
        ++openfiles;
        fptraincontrol = fopen(traincontrol,"wb");
        ++openfiles;
        }

    if(clumpcnt && fptrain)
        {
        fptest = fopen(test,"wb");
        ++openfiles;
        fpcontrol = fopen(control,"wb");
        ++openfiles;
        if(fptest && fpcontrol && (!TrainTest || (fptraintest && fptraincontrol)))
            {
            while(testclumps > 0 || trainclumps > 0)
                {
                int k = 0;
                while(k < clumpcnt && (testclumps > 0 || trainclumps > 0))
                    {
                    int rd = Rand();
                    int L = 0;
                    if(  testclumps == 0 
                      || (  trainclumps > 0
                         && (double)rd * 10000.0 < (double)RAND_MAX * (double)fraction
                         )
                      ) // write to training file
                        {
                        assert(trainclumps > 0);
                        for(;L < clumps[k].linecnt;++L)
                            printTrain(clumps[k].start+L,fptrain,fptraintest,fptraincontrol,columntest,columncontrol,columntag,TrainTest);
                        trainlines += clumps[k].linecnt;
                        --trainclumps;
                        }
                    else
                        {
                        for(;L < clumps[k].linecnt;++L)
                            writeTestAndControl(clumps[k].start+L,fptest,fpcontrol,columntest,columncontrol,columntag);
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
    if(fptrain){--openfiles;fclose(fptrain);}
    return trainlines;
    }

static int splitLemmaList
   (int linecnt
   ,char * training
   ,char * test
   ,char * control
   ,int columntest
   ,int columncontrol
   ,int columntag
   ,int fraction
   ,char * traintest
   ,char * traincontrol
   ,line * lines
   ,bool TrainTest
   ,bool TenFoldXValidation
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
    int (*Rand)(void);
    if(TenFoldXValidation)
        {
        fraction = 9000;
        Rand = aTenth;
        }
    else
        Rand = rand;

    if(fraction < 10000)
        trainlines = (int)((double)linecnt * (double)fraction / 10000.0);
    else
        trainlines = linecnt;

    testlines = linecnt - trainlines;
    ret = trainlines;
    if(TrainTest)
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
            if(fptest && fpcontrol && (!TrainTest || (fptraintest && fptraincontrol)))
                {
                while(testlines > 0 || trainlines > 0)
                    {
                    int k = 0;
                    while(k < linecnt && (testlines > 0 || trainlines > 0))
                        {
                        assert(!IsSpace(lines[k].s[0]));
                        int rd = Rand();
                        if(  testlines == 0 
                          || (  trainlines > 0
                             && (double)rd * 10000.0 < (double)RAND_MAX * (double)fraction && trainlines > 0
                             )
                          ) // write to training file
                            {
                            assert(trainlines > 0);
                            printTrain(lines+k,fptrain,fptraintest,fptraincontrol,columntest,columncontrol,columntag,TrainTest);
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
                printTrain(lines+k,fptrain,fptraintest,fptraincontrol,columntest,columncontrol,columntag,TrainTest);
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

struct evaluation
    {
    int same;
    int different;
    int ambiguousRules;
    ambty ambiguous;
    decision Decision;
    evaluation(): same(0),different(0),ambiguousRules(0)
        {
        for(size_t j=0;j < sizeof(ambiguous)/sizeof(ambiguous[0]);++j)
            ambiguous[j] = 0;
        }
    };

static evaluation compare(const char * output, const char * control, const char * controlResult, const char * test, optionStruct * options)
    {
    evaluation Evaluation;
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
            char * tab;
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
                if(hash && hash > s && '0' <= hash[1] && hash[1] <= 9)
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
                    if(bar && bar > s && bar[1])
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
                    ++Evaluation.ambiguous[c];
                    Correct(Ref);
                    if(f3)
                        fprintf(f3,"|\t%ld\t%s\t%s\t%s\n",Ref,b4,b1,b2);
                    }
                else
                    {
                    assert(!ambiguousRule);
                    ++Evaluation.same;
                    Correct(Ref);
                    if(f3)
                        fprintf(f3,"+\t%ld\t%s\t%s\t%s\n",Ref,b4,b1,b2);
                    }
                }
            else
                {
                ++Evaluation.different;
                Wrong(Ref);
                if(f3)
                    fprintf(f3,"-\t%ld\t%s\t%s\t%s\n",Ref,b4,b1,b2);
                }

            if(ambiguousRule)
                {
                ++Evaluation.ambiguousRules;
                if(bambiguous)
                    Evaluation.Decision.true_amb++;
                else
                    Evaluation.Decision.false_amb++;
                }
            else
                {
                if(bambiguous)
                    Evaluation.Decision.false_not_amb++;
                else
                    Evaluation.Decision.true_not_amb++;
                }
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
    if(options->remove())
        {
        remove(output);
        remove(controlResult);
        }
    return Evaluation;
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

static void printResults(FILE * fptally,int flexcount,evaluation Evaluation)
    {
    int iall = Evaluation.same + Evaluation.ambiguous[0] + Evaluation.ambiguous[1] + Evaluation.ambiguous[2] +  Evaluation.different;
    if(iall > 0)
        {
        double all = 0.01*iall;
        assert(all > 0);
        fprintf
            (fptally,"#rules:%d #same:%d #ambiguous:%d %d %d #different:%d %14.6f %14.6f %14.6f %14.6f %14.6f (amb %d all %d amb%% %14.6f)\n"
            ,flexcount
            ,Evaluation.same
            ,Evaluation.ambiguous[0]
            ,Evaluation.ambiguous[1]
            ,Evaluation.ambiguous[2]
            ,Evaluation.different
            ,(double)Evaluation.same/all
            ,(double)Evaluation.ambiguous[0]/all
            ,(double)Evaluation.ambiguous[1]/all
            ,(double)Evaluation.ambiguous[2]/all
            ,(double)Evaluation.different/all
            ,Evaluation.ambiguousRules
            ,iall
            ,(double)Evaluation.ambiguousRules/all
            );
        }
    else
        fprintf
            (fptally,"same + ambiguous[0] + ambiguous[1] + ambiguous[2] +  different = %d\n"
            ,Evaluation.same + Evaluation.ambiguous[0] + Evaluation.ambiguous[1] + Evaluation.ambiguous[2] +  Evaluation.different
            );
    fflush(fptally);
    }

typedef struct countStruct
    {
    int tflexcount;
    int tsame;
    int tdifferent;
    ambty tambiguous;
    evaluation Evaluation;
    int tambiguousRules;
    countStruct()
        {
        tflexcount = tsame = tdifferent = tambiguousRules = 0;
        for(size_t i = 0;i < sizeof(tambiguous)/sizeof(tambiguous[0]);++i)
            {
            tambiguous[i] = 0;
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
            , const char * output
            , const char * traintest
            , const char * test
            , const char * traincontrol
            , const char * controlResult
            , const char * control
            , FILE * fptally
            , countAndWeight * Counts
            , optionStruct * Options
            , bool TrainTest
            )
            {
            long nflexcount = 0;
            nflexcount = Counts[cutoff].getNnodes();
            n.tflexcount += nflexcount;
            if(TrainTest)
                {
                n.Evaluation = compare(output,traincontrol,controlResult,traintest,Options);
                }
            else
                {
                n.Evaluation = compare(output,control,controlResult,test,Options);
                }
            this->StandardDev.datum(n.Evaluation.same,n.Evaluation.ambiguous,n.Evaluation.different);
            n.tsame += n.Evaluation.same;
            n.tdifferent += n.Evaluation.different;
            for(size_t j = 0;j < sizeof(n.tambiguous)/sizeof(n.tambiguous[0]);++j)
                {
                n.tambiguous[j] += n.Evaluation.ambiguous[j];
                }
            n.tambiguousRules += n.Evaluation.ambiguousRules;
            fprintf(fptally,"cutoff:%d ",cutoff);
            printResults(fptally,nflexcount,n.Evaluation);
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
                ,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.false_amb/ntot : 0
                ,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.false_not_amb/ntot : 0
                ,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.true_amb/ntot : 0
                ,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.true_not_amb/ntot : 0
                );
            fprintf(fptab,  "%14.6f %14.6f\n"
                ,this->n.Evaluation.Decision.precision()
                ,this->n.Evaluation.Decision.recall()
                );
            fflush(fptab);
            }
        void validationcolumn(const char ** cell, int cutoff, int maxcount, int ttrainlines, lineab * AffixLine, int norow)
            {
            double ntot = n.tsame + n.tambiguous[0] + n.tambiguous[1] + n.tambiguous[2] + n.tdifferent;
            const char * f1 = "%14d ";
            const char * f2 = "%14.6f ";
            char ** ell = new char *[norow];
            for(int i = 0;i < norow;++i)
                cell[i] = ell[i] = new char[22];
            sprintf(ell[0],f1,cutoff);
            sprintf(ell[1],f2,(double)n.tflexcount/(double)maxcount);
            sprintf(ell[2],f2,ttrainlines == 0 ? 100.0 : 100.0*(double)n.tflexcount/(double)ttrainlines);
            sprintf(ell[3],f2,ntot > 0 ? 100.0*(double)n.tsame/ntot : 0.0);
            sprintf(ell[4],f2,ntot > 0 ? 100.0*(double)n.tambiguous[0]/ntot : 0.0);
            sprintf(ell[5],f2,ntot > 0 ? 100.0*(double)n.tambiguous[1]/ntot : 0.0);
            sprintf(ell[6],f2,ntot > 0 ? 100.0*(double)n.tambiguous[2]/ntot : 0.0);
            sprintf(ell[7],f2,ntot > 0 ? 100.0*(double)n.tdifferent/ntot : 0.0);
            sprintf(ell[8],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(esam) : 0.0);
            sprintf(ell[9],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb0) : 0.0);
            sprintf(ell[10],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb1) : 0.0);
            sprintf(ell[11],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(eamb2) : 0.0);
            sprintf(ell[12],f2,ntot > 0 ? 100.0*this->StandardDev.calculate(edif) : 0.0);
            sprintf(ell[13],f2,ntot > 0 ? 100.0*(double)n.tambiguousRules/ntot : 0.0);
            sprintf(ell[14],f2,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.false_amb/ntot : 0.0);
            sprintf(ell[15],f2,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.false_not_amb/ntot : 0.0);
            sprintf(ell[16],f2,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.true_amb/ntot : 0.0);
            sprintf(ell[17],f2,ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.true_not_amb/ntot : 0.0);
            sprintf(ell[18],f2,this->n.Evaluation.Decision.precision());
            sprintf(ell[19],f2,this->n.Evaluation.Decision.recall());
            sprintf(ell[20],"%6.3f*N^%4.3f ", exp(AffixLine->a()), AffixLine->b());//0.056414*N^0.799693
            delete[]ell;
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
            fprintf(fptab,"false_amb%%     %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.false_amb/ntot : 0.0);
            fprintf(fptab,"false_not_amb%% %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.false_not_amb/ntot : 0.0);
            fprintf(fptab,"true_amb%%      %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.true_amb/ntot : 0.0);
            fprintf(fptab,"true_not_amb%%  %14.6f\n",ntot > 0 ? 100.0*(double)this->n.Evaluation.Decision.true_not_amb/ntot : 0.0);
            fprintf(fptab,"precision      %14.6f\n",this->n.Evaluation.Decision.precision());
            fprintf(fptab,"recall         %14.6f\n",this->n.Evaluation.Decision.recall());
            fflush(fptab);
            }
    };


void createReport(counting c[CUTOFFS],int maxcount,int ttrainlines,lineab AffixLine[CUTOFFS],int noOfSteps,int lastFraction,optionStruct * Options)
    {
    const char * texts[] =    {
        "prun. thrshld. ",
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
        "\n"
        ";Evaluation of prediction of ambiguity (whether a word has more than one possible lemma)\n"
        ";---------------------------------------------------------------------------------------\n"
        "; amb.rules%     ",
        "false_amb%     ",
        "false_not_amb% ",
        "true_amb%      ",
        "true_not_amb%  ",
        "precision      ",
        "recall         ",
        "\n"
        "; Power law relating the number of rules in the decision tree to the number of examples in the training data\n"
        ";----------------------------------------------------------------------------------------------------------\n"
        "; #rules =       "
        };
    int nocol = 7;
    int norow = sizeof(texts)/sizeof(texts[0]);
    int presentRows = norow;
    const char * introduction;
    char postScriptum[500];
    if(AffixLine[0].computable())
        {
        introduction = "; Lemmatization results for data that is not part of the training data.\n";
            
        sprintf(postScriptum,
            "; The number of rules can be estimated from the number of training examples by\n"
            "; a power law. See the last line in the table above, which is based on %d\n"
            "; different samples from the total available training data mass varying in size\n"
            "; from %.2f %% to %.2f %%",noOfSteps,FRACTION_LOW/100.0,lastFraction/100.0);
        }
    else
        {
        introduction = 
            "; Lemmatization results for all data in the training set.\n"
            "; For pruning threshold 0 there may be no errors (diff%).\n"
            "; (This is not necessarily true for external trainers/lemmatizers that\n"
            "; e.g. do not return all possible solutions in case of ambiguity.)\n"
            ;
        postScriptum[0] = 0;
        --presentRows;
        }
    int nocell = norow*nocol;
    const char * cell[sizeof(texts)/sizeof(texts[0]) * 7];
    for(int j = 0;j < norow;++j)
        cell[j] = texts[j];
    for(int k = 0;k < nocol - 1;++k)
        {
        c[k].validationcolumn(cell + (k + 1)*norow, k, maxcount, ttrainlines, AffixLine + k, norow);
        }
    int len = 0;
    for(int m = 0;m < nocell;++m)
        len += strlen(cell[m]);
    char * report = new char[len + 3*norow + 1];
    report[0] = 0;
    for(int row = 0;row < presentRows;++row)
        {
        strcat(report,"; ");
        for(int col = 0;col < 7;++col)
            {
            int ind = row + col*norow;
            strcat(report,cell[ind]);
            }
        strcat(report,"\n");
        }
    Options->printEvaluation(introduction,report,postScriptum);
    delete [] report;
    for(int row2 = 0;row2 < norow;++row2)
        {
        for(int col = 1;col < 7;++col)
            {
            delete [] cell[row2+col*norow];
            }
        }
    }

int nextFraction(int fraction)
    {
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
    return fraction;
    }


void trainAndTest
        (int linecnt
        ,int clumpcnt
        ,const char * XT
        ,const char * TT
        ,line * lines
        ,clump * clumps
        ,optionStruct * Options
        ,bool TrainTest
        )
    {
    lineab AffixLine[CUTOFFS];

    char formatprefix[256]          ;sprintf(formatprefix,          BASEDIR "%s%%s" SEPARATOR "%s_%s%s%s%s%s_%%s.txt",Options->tempDir(),LGf(Options),XTRf(Options),(Options->suffixOnly() ? "_suffix" : "_affix"),XT,TT,(Options->redo() ? "redone" : "singleshot"));

    char formatTraining[256]        ;sprintf(formatTraining,        formatprefix,"training"     ,"%d_%d");        // the training words
    char formatTest[256]            ;sprintf(formatTest,            formatprefix,"test"         ,"%d_%d");        // the test words (<> training words)
    char formatTrainTest[256]       ;sprintf(formatTrainTest,       formatprefix,"test"         ,"_train_%d_%d"); // the test words (= training words)
    char formatControl[256]         ;sprintf(formatControl,         formatprefix,"control"      ,"%d_%d");        // the correct answers (<> training words)
    char formatTrainControl[256]    ;sprintf(formatTrainControl,    formatprefix,"control"      ,"_train_%d_%d"); // the correct answers (= training words)

    char formatFlexrules[256]       ;sprintf(formatFlexrules,       formatprefix,"flex"         ,"%d_%d");
    char formatSOutput[256]         ;sprintf(formatSOutput,         formatprefix,"sout"         ,"%d_%d_%d");
    char formatOutput[256]          ;sprintf(formatOutput,          formatprefix,"out"          ,"%d_%d_%d");

    char formatcontrolResult[256]  ;sprintf(formatcontrolResult,  formatprefix,"resultAffix"  ,"%d_%d_cutoff_%d");

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
    int noOfSteps;
    for(fraction = FRACTION_LOW,noOfSteps = 1;fraction <= FRACTION_HIGH;fraction = nextFraction(fraction),++noOfSteps)
        {
        if (Options->verbose())
            printf("Test: take fraction %d\n",fraction);
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
        if(Options->tenfoldCrossValidation())
            {
            maxcount = 10;
            }
        else
            {
            maxcount = (int)(5000000.0*sqrt(1.0/((double)fraction*(10000.0-fraction)*(double)linecnt)));
            if(maxcount < 1)
                maxcount = 1;
            else if(maxcount > 30)
                maxcount = 30;
            }
        char paramfile[100];
        sprintf(paramfile,"%s%s%s%s%s.txt",LGf(Options),(Options->suffixOnly() ? "_suffix" : "_affix"),XTRf(Options),TT,(Options->redo() ? "redone" : "singleshot"));
        globmaxcount = maxcount;
        counting c[CUTOFFS]; // indexed by cutoff
        countAndWeight * Counts = new countAndWeight[CUTOFF_HIGH+1];
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
            char Soutput[256];
            char output[256];
            if(Options->tenfoldCrossValidation())
                tenth(count - 1);
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
            if(TrainTest )
                {
                sprintf(traintest,formatTrainTest,fraction,count);
                sprintf(Ttraincontrol,formatTrainControl,fraction,count);
                }
            static char suffixrulesArr[256];
            sprintf(suffixrulesArr,formatFlexrules,fraction,count);

            int trainlines = 0;
            if(clumpcnt > 1)
                trainlines = splitLemmaList(clumpcnt,Ttraining,test,Tcontrol,1,2,3,fraction,traintest,Ttraincontrol,clumps,TrainTest,Options->tenfoldCrossValidation());
            else
                trainlines = splitLemmaList(linecnt,Ttraining,test,Tcontrol,1,2,3,fraction,traintest,Ttraincontrol,lines,TrainTest,Options->tenfoldCrossValidation());
            if(trainlines > 0)
                {
                ttrainlines += trainlines;

                optionStruct testOptions(*Options);
                testOptions.seti(Ttraining);
                testOptions.setc(CUTOFF_HIGH);
                testOptions.seto(Options->flexrules());
                testOptions.sete(LGf(Options));
                testOptions.setn("123");
                testOptions.setf(XTRf(Options));
                testOptions.setP(paramfile);

                testOptions.completeArgs();
                trainRules("",&testOptions,Counts);
                }
            if (Options->remove())
                {
                remove(Ttraining);
                }

            if(ttrainlines > 0)
                {
                for(int cutoff = CUTOFF_LOW;cutoff <= CUTOFF_HIGH;++cutoff)
                    {
                    sprintf(controlResult,formatcontrolResult,fraction,count,cutoff);
                    char Affixrules[250];
                    const char * lastslash = strrchr(Options->flexrules(),*SLASH);
                    const char * filename;
                    if(lastslash)
                        {
                        filename = lastslash + 1;
                        sprintf(Affixrules,"%.*s%d%s%s",(int)(filename - Options->flexrules()),Options->flexrules(),cutoff,SLASH,filename);
                        }
                    else
                        sprintf(Affixrules,"%d%c%s",cutoff,DIRSEP,Options->flexrules());
                    

                    sprintf(Soutput,formatSOutput,cutoff,fraction,count);
                    sprintf(output,formatOutput,cutoff,fraction,count);
                    FILE * fptally = fopen(tally,"ab");
                    ++openfiles;
                    if(fptally)
                        {
                        if(Options->externalLemmatizer())
                            {
                            char * command = new char[strlen(Options->externalLemmatizer())+strlen(TrainTest ? traintest : test)+strlen(Affixrules)+strlen(output)+5];
                            sprintf(command,"%s %s %s %s",Options->externalLemmatizer(),TrainTest ? traintest : test,Affixrules,output);
                            if(system(command))
                                fprintf(stderr,"Cannot execute system command \"%s\".\n",command);
                            delete [] command;
                            }
                        else
                            lemmatiseFile(TrainTest ? traintest : test,Affixrules,output);
                        if (Options->remove()) { remove(Affixrules); }
                        c[cutoff].testing(cutoff,output,traintest,
                            test,Ttraincontrol,controlResult,Tcontrol,fptally,Counts,Options,TrainTest);
                        --openfiles;
                        fclose(fptally);
                        }
                    else
                        {
                        fprintf(stderr,"Cannot open \"%s\" for appending. No testing is done.\n",tally);
                        }
                    }
                }
            if (Options->remove())
                {
                remove(test);
                if(TrainTest)
                    {
                    remove(traintest);
                    remove(Ttraincontrol);
                    }
                remove(Tcontrol);
                }
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
                    fprintf (fptab,"cutoff %d Affix  a %14.6f b %14.6f: N(rules)=%14.6f*N(trainpairs)^%f\n"
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
                if(nextFraction(fraction) > FRACTION_HIGH)
                    {
                    createReport(c,maxcount,ttrainlines,AffixLine,noOfSteps,fraction,Options);
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
       
        delete [] Counts;
        }
    }


static int readFile
        (line *& lines
        ,clump *& clumps
        ,int & clumpcnt
        ,int sep
        ,optionStruct * Options
        )
    {
    return readlines(1,2,3,lines,sep,clumps,clumpcnt,Options);
    }


static int readFileAndTrainAndTest(optionStruct * Options,bool TrainTest)
    {
    const char * XT;
    const char * TT;
    if (Options->verbose())
        printf("Test initialization\n");
    CUTOFF_LOW = 0;
    CUTOFF_HIGH = Options->cutoff();
    TAGGED = false;
    if(Options->tenfoldCrossValidation())
        {
        FRACTION_LOW = 9000;
        FRACTION_HIGH = FRACTION_LOW;
        MAXCOUNT = 10;
        TrainTest = false;
        }
    else
        {
        switch(STEPSIZE)
            {
            case 100:
                FRACTION_LOW  = 154;
                break;
            case 25:
                FRACTION_LOW  = 171;
                break;
            default:
                FRACTION_LOW  = 1000;
                ;
            }
        FRACTION_HIGH = 10000;
        MAXCOUNT = 100;
        }

    if(TrainTest)
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
    int sep = '\t';

    if(XTRf(Options) == NULL)
        {
        printf("Rule comparison function (variable XTRf(Options)) not specified for language %s\n",LGf(Options));
        exit(1);
        }
    
    line * lines = NULL;
    clump * clumps = NULL;
    int clumpcnt = 0;
    if (Options->verbose())
        printf("Test: read file\n");
    int linecnt = readFile
        (lines
        ,clumps
        ,clumpcnt
        ,sep
        ,Options
        );
    if (Options->verbose())
        printf("Test: read file DONE\n");
    if(linecnt)
        {
        if(Options->verbose())
            printf("Test: train and then test\n");
        trainAndTest
            (linecnt
            ,clumpcnt
            ,XT
            ,TT
            ,lines
            ,clumps
            ,Options
            ,TrainTest
            );

        if (Options->verbose())
            printf("Test: train and then test DONE\n");
        delete [] lines;
        lines = NULL;
        }
    return 0;
    }

int testrules(optionStruct * Options)
    {
    if(lemmalistef(Options) && LGf(Options) && XTRf(Options))
        {
        if(Options->trainTest() && !Options->tenfoldCrossValidation())
            {
            if (Options->verbose())
                printf("trainTest.\n");
            readFileAndTrainAndTest(Options,true);
            }
        if(Options->test())
            {
            if (Options->verbose())
                printf("OOVTest.\n");
            readFileAndTrainAndTest(Options,false);
            }
        return 0;
        }

    return 0;
    }
