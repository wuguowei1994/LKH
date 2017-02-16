#include "LKH.h"

/*
  Minimum1TreeCost()函数会返回最小1-tree的cost
  叶节点:没有子节点的点。
  分支节点:除根节点外，有子节点的点。
  最小1-tree的概念见论文第18页

  一个节点的V-value等于它的度数减去2
  Norm等于所有节点V-value的平方和.Norm主要用来评测这个最小1-tree是否存在满足条件的最优环路。
  如果Norm为零，那么1-tree就包含一个满足条件的最优环路，否则不是。
 */

/*
   实参说明:程序第一次和最后一次调用这个函数时传入的实参为0。其他阶段实参全为1
 */ 
GainType Minimum1TreeCost(int Sparse)
{
    Node *N, *N1 = 0;
    GainType Sum = 0;
    int Max = INT_MIN;

    MinimumSpanningTree(Sparse);
    N = FirstNode;
    //把所有节点的V初始化为-2;所有节点的Pi初始化为0
    do {
        N->V = -2;
        Sum += N->Pi;
    }
    while ((N = N->Suc) != FirstNode);
    Sum *= -2;
    while ((N = N->Suc) != FirstNode) {
        N->V++;
        N->Dad->V++;
        Sum += N->Cost;
        N->Next = 0;
    }
    FirstNode->Dad = FirstNode->Suc;
    FirstNode->Cost = FirstNode->Suc->Cost;
    do {
        if (N->V == -1) {
            /*
             Connect(N, Max, Sparse)
             参数说明:N是最小生成树上度数为1的节点，Max是查找过程中边权重的临界值(一旦找到一条小于Max的边就终止查找)
             功能说明:确定一条从N节点出发的最短边(这条边不在原最小生成树中)，N1的Next指针将指向这条最短边的终点，N1的NextCost记录
             了这条最短边的权重
            */
            Connect(N, Max, Sparse);
            if (N->NextCost > Max) {
                N1 = N;
                Max = N->NextCost;
            }
        }
    }
    while ((N = N->Suc) != FirstNode);
    N1->Next->V++;
    N1->V++;
    Sum += N1->NextCost;
    Norm = 0;
    //Norm等于所有节点V-value的平方和
    do
        Norm += N->V * N->V;
    while ((N = N->Suc) != FirstNode);
    if (N1 == FirstNode)
        // 如果N1为FirstNode，就把N1置为0
        N1->Suc->Dad = 0;
    else {
        FirstNode->Dad = 0;
        //把节点N1置于节点FirstNode的前面
        Precede(N1, FirstNode);
        FirstNode = N1;
    }
    // 这个if语句没有满足过(从始至终)
    if (Norm == 0) {
        for (N = FirstNode->Dad; N; N1 = N, N = N->Dad)
            Follow(N, N1);
        for (N = FirstNode->Suc; N != FirstNode; N = N->Suc)
            N->Dad = N->Pred;
        FirstNode->Suc->Dad = 0;
    }
    return Sum;
}
