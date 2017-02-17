#include "LKH.h"

/* 
 * The Forbidden function is used to test if an edge, (ta,tb), 
 * is one of the forbidden edges (C(ta, tb) == M) in a solution of
 * an asymmetric traveling saleman problem. 
 * If the edge is forbidden, the function returns 1; otherwise 0.
 */
// Forbidden()函数用来判断atsp问题中，(ta,tb)这条边是否是一条权重为无穷的边
// 如果权重为无穷，就会返回1;否则返回0
// 补充:如果问题为tsp问题，会直接返回0
int Forbidden(const Node * ta, const Node * tb)
{
    return ProblemType == ATSP &&
        (ta->Id <= DimensionSaved) == (tb->Id <= DimensionSaved);
}
