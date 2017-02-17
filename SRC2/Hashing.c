#include "Hashing.h"

// HashInitialize(T)会把传入的hash table清空
void HashInitialize(HashTable * T)
{
    int i;

    for (i = 0; i < HashTableSize; i++) {
        T->Entry[i].Hash = UINT_MAX;
        T->Entry[i].Cost = MINUS_INFINITY;
    }
    T->Count = 0;
}

/*
  HashInsert(T,H,Cost)函数向哈希表T中插入H(key)和Cost(value)
  Cost表示可行解的权重
  用双散列法解决冲突
  如果原来的H(key)已经存有cost，那么这样处理：如果新的cost大于等于原来的cost，就会把
  原来的那个cost替换掉。
 */
void HashInsert(HashTable * T, unsigned Hash, GainType Cost)
{
    int i = Hash % HashTableSize;
    if (T->Count >= MaxLoadFactor * HashTableSize) {
        if (Cost > T->Entry[i].Cost)
            return;
    } else {
        int p = Hash % 97 + 1;
        while (T->Entry[i].Cost != MINUS_INFINITY)
            if ((i -= p) < 0)
                i += HashTableSize;
        T->Count++;
    }
    T->Entry[i].Hash = Hash;
    T->Entry[i].Cost = Cost;
}

/*
   HashSearch(T,H,Cost)
   如果hash表T中含有H(key)和Cost(value)(Cost代表可行解的权重),这个函数会返回1
   否则返回0
 */

int HashSearch(HashTable * T, unsigned Hash, GainType Cost)
{
    int i, p;

    i = Hash % HashTableSize;
    p = Hash % 97 + 1;
    while ((T->Entry[i].Hash != Hash || T->Entry[i].Cost != Cost)
            && T->Entry[i].Cost != MINUS_INFINITY)
        if ((i -= p) < 0)
            i += HashTableSize;
    return T->Entry[i].Hash == Hash;
}
