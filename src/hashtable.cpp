#include "hashtable.h"
#include "vertex.h"

int HashCount = 0;

long nextprime(unsigned long g)
    {
    int i;
    unsigned long smalldivisor;
    static int byte[12]=
        {1,  2,  2,  4,    2,    4,    2,    4,    6,    2,  6};
    /*2-3,3-5,5-7,7-11,11-13,13-17,17-19,19-23,23-29,29-1,1-7*/
    unsigned long bigdivisor;
    if(g & 1)
        ++g; /* even -> uneven */
    smalldivisor = 2;
    i = 0;
    while((bigdivisor = g / smalldivisor) >= smalldivisor)
        {
        if(bigdivisor * smalldivisor == g)
            {
            ++g;
            smalldivisor = 2;
            i = 0;
            }
        else
            {
            smalldivisor += byte[i];
            if(++i > 10)
                i = 3;
            }
        }
    return g;
    }

void hashTable::rehash(int loadFactor/*1-100*/)
    {
    long oldsize = hash_size;
    hash_size = nextprime((100 * record_count)/loadFactor);
    vertex ** new_hash_table = new vertex * [hash_size];
    unsigned long i;
    for(i = 0;i < hash_size;++i)
        new_hash_table[i] = 0;
    if(hash_table)
        {
        for(i = oldsize;i > 0;)
            {
            vertex * r = hash_table[--i];
            while(r)
                {
                long Key = key(r->itsTxt());
                vertex ** head = new_hash_table + Key;
                vertex * Next = r->getNext();
                r->setNext(*head);
                r->setHead(head);
                *head = r;
                r = Next;
                }
            }
        }
    delete [] hash_table;
    hash_table = new_hash_table;
    }

hashTable::hashTable(long size):record_count(0),next(0)
    {
    hash_size = nextprime(size);
    hash_table = new vertex * [hash_size];
    unsigned long i;
    for(i = 0;i < hash_size;++i)
        hash_table[i] = 0;
    ++HashCount;
    }

hashTable::~hashTable()
    {
    unsigned long i;
    for(i = 0;i < hash_size;++i)
        {
        vertex * ps = hash_table[i];
        while(ps)
            {
            ps->Hash = 0;
            ps = ps->getNext();
            }
        }
    delete [] hash_table;
    // do we want to delete the strngs as well?
    delete next;
    --HashCount;
    }

long casesensitivehash(const char * cp)
    {
    long hash_temp = 0;
    while (*cp != '\0')
        {
        if(hash_temp < 0)
            hash_temp = (hash_temp << 1) +1;
        else
            hash_temp = hash_temp << 1;
        hash_temp ^= *cp;
        ++cp;
        }
    return hash_temp;
    }

long hashTable::key(const char * ckey)
    {
    return ((unsigned long)casesensitivehash(ckey)) % hash_size;
    }

vertex * hashTable::find(const char * ckey,vertex **& head)
    {
    long Key = key(ckey);
    head = hash_table + Key;
    vertex * p;
    for(p = *head;p && strcmp(p->itsTxt(),ckey);p = p->getNext())
        ;
    return p;
    }

vertex ** hashTable::convertToList()
    {
    vertex ** ret = new vertex * [record_count];
    unsigned long i;
    unsigned long n = 0;
    for(i = 0;i < hash_size;++i)
        {
        vertex * ps = hash_table[i];
        while(ps)
            {
            ps->incRefCnt();
            ret[n++] = ps;
            ps = ps->getNext();
            }
        }
    return ret;
    }

unsigned long hashTable::forall(forallfunc fnc)
    {
    unsigned long i;
    unsigned long n = 0;
    for(i = 0;i < hash_size;++i)
        {
        vertex * ps = hash_table[i];
        while(ps)
            {
            ++n;
            (ps->*fnc)();
            ps = ps->getNext();
            }
        }
    return n;
    }


vertex * hashTable::getVertex(vertex * Rule,bool & New)
    {
    vertex * ret = 0;
    New = false;
    vertex ** head;
    long lf = loadfactor();
    if(lf > 100)
        {
        rehash(60);
        }
    vertex * p = find(Rule->Pattern.itsTxt(),head);
    if(p)
        p = p->findReplacement(Rule);
    if(p)
        {
        ret = p;
        }
    else
        {
        ret = new vertex(Rule,this);
        ret->setNext(*head);
        *head = ret;
        ret->setHead(head);
        New = true;
        inc();
        }
    return ret;
    }
