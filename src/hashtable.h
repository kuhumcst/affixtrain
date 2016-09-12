#ifndef HASHTABLE_H
#define HASHTABLE_H

class vertex;
typedef  void (vertex::*forallfunc)();

class hashTable
    {
    private:
        unsigned long hash_size;
        unsigned long record_count;
        vertex ** hash_table;
        hashTable * next;
    public:
        unsigned long recordCount()const
            {
            return record_count;
            }
        vertex ** convertToList();
        hashTable * getNext()
            {
            return next;
            }
        hashTable *& getNextR()
            {
            return next;
            }
        void inc()
            {
            ++record_count;
            }
        void dec()
            {
            --record_count;
            }
        void rehash(int loadFactor/*1-100*/);
        unsigned long loadfactor()
            {
            if(record_count < 10000000UL)
                return (100 * record_count) / hash_size;
            else
                return record_count / (hash_size/100);
            }
        hashTable(long size);
        ~hashTable();
        long key(const char * ckey);
        vertex * find(const char * ckey,vertex **& head);
        unsigned long forall(forallfunc fnc);
        vertex * getVertex(vertex * Rule,bool & New);
        bool deleteVertex(vertex * Rule);
    };

#endif
