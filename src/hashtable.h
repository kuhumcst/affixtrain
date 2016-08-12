#ifndef HASHTABLE_H
#define HASHTABLE_H

class vertex;
typedef  void (vertex::*forallfunc)();

class hashTable
    {
    private:
        long hash_size;
        long record_count;
        vertex ** hash_table;
        hashTable * next;
    public:
        long recordCount()const
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
        long loadfactor()
            {
            if(record_count < 10000000L)
                return (100 * record_count) / hash_size;
            else
                return record_count / (hash_size/100);
            }
        hashTable(long size);
        ~hashTable();
        long key(const char * ckey);
        vertex * find(const char * ckey,vertex **& head);
        int forall(forallfunc fnc);
        vertex * getVertex(vertex * Rule,bool & New);
        bool deleteVertex(vertex * Rule);
    };

#endif
