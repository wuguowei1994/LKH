#include "LKH.h"

/* 
 * The IsCommonEdge function is used to test if an edge, (ta,tb), 
 * is common to the set of tours to be merged.
 *
 * If the edge is a common edge the function returns 1; otherwise 0.
 */
//IsCommonEdge()用来检测(ta,tb)这条边是否被一条需要合并的路径包含
int IsCommonEdge(const Node * ta, const Node * tb)
{
    int i;

    if (MergeTourFiles < 2)
        return 0;
    for (i = 0; i < MergeTourFiles; i++)
        if (ta->MergeSuc[i] != tb && tb->MergeSuc[i] != ta)
            return 0;
    return 1;
}
