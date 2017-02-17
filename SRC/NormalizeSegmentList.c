#include "LKH.h"

/*
  NormalizeSegmentList()函数用来交换一些小段的Suc和Pred，然后形成一个循环的双向链表
  这个函数会被LKH算法调用
 */

void NormalizeSegmentList()
{
    Segment *s1, *s2;

    s1 = FirstSegment;
    do {
        if (!s1->Parent->Reversed)
            s2 = s1->Suc;
        else {
            s2 = s1->Pred;
            s1->Pred = s1->Suc;
            s1->Suc = s2;
        }
    }
    while ((s1 = s2) != FirstSegment);
}
