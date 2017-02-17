#include "Segment.h"
#include "LKH.h"

/*
 NormalizeNodeList()函数用来交换节点的Suc和Pred，然后得到一个循环的双向链表
 */

void NormalizeNodeList()
{
    Node *t1, *t2;

    t1 = FirstNode;
    do {
        t2 = SUC(t1);
        t1->Pred = PRED(t1);
        t1->Suc = t2;
    }
    while ((t1 = t2) != FirstNode);
}
