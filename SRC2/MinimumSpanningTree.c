#include "LKH.h"
#include "Heap.h"

/*
 * The MinimumSpanningTree function determines a minimum spanning tree using 
 * Prim's algorithm.
 *
 * At return the Dad field of each node contains the father of the node, and 
 * the Cost field contains cost of the corresponding edge. The nodes are 
 * placed in a topological ordered list, i.e., for any node its father precedes 
 * the node in the list. The fields Pred and Suc of a node are pointers to the 
 * predecessor and successor node in this list.
 *
 * The function can be used to determine a minimum spanning tree in a dense 
 * graph, or in a sparse graph (a graph determined by a candidate set).
 *
 * When the graph is sparse a priority queue, implemented as a binary heap, 
 * is used  to speed up the determination of which edge to include next into 
 * the tree. The Rank field of a node is used to contain its priority (usually 
 * equal to the shortest distance (Cost) to nodes of the tree).        
 */
/*
  MinimumSpanningTree()函数使用Prim算法计算一个最小生成树
  每个点会得到一个父节点，同时用cost数组记录相应的权重，这些节点被放置在一个拓扑结构的链表中。
  换句话说，每一个节点的父节点就是它在这个拓扑结构的链表中的前一个节点。指针pred指向了这个节点的前面一个节点
  指针suc指向了这个节点的后一个节点。
  这个函数不仅可以被用来计算稠密图的最小生成树，也可以用来计算稀疏图的最小生成树。(请注意，这里的图不是以输入文件为基础，而是以候选集为基础的)
 */
void MinimumSpanningTree(int Sparse)
{
    Node *Blue;         /* 指向的是被加入到最小生成树中的最后一个节点*/
    Node *NextBlue = 0; /* 指向的是被加入到最小生成树中的最后一个节点(临时的记录一下)*/
    Node *N;
    Candidate *NBlue;  /*Candidate是一个结构体，表示候选边，它包含指针*To，表示这条边的终点，Cost(int)表示这条边的权重，Alpha(int)表示它的alpha值*/
    int d;
    Blue = N = FirstNode;
    Blue->Dad = 0;              /* The root of the tree has no father */
    /*
    Bule->CandidateSet：只在程序初始化的时候为空
    Sparse：只在程序初始化和结束时为0
    下面的这个if语句除了在程序初始化和程序结束的时候不满足，其他情况下都满足
    CandidateSet是一个指针，它指向了Bule这个节点的Candidate*/
    if (Sparse && Blue->CandidateSet) {
        //这张图是稀疏图
        //把所有的节点压入堆栈
        Blue->Loc = 0;          //blue节点不在堆栈中
        while ((N = N->Suc) != FirstNode) {
            N->Dad = Blue;
            N->Cost = N->Rank = INT_MAX;
            HeapLazyInsert(N);
        }
        //更新blue节点的所有邻居
        for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
            if (FixedOrCommon(Blue, N)) {
                N->Dad = Blue;
                N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                N->Rank = INT_MIN;
                HeapSiftUp(N);
            } else if (!Blue->FixedTo2 && !N->FixedTo2) {
                N->Dad = Blue;
                N->Cost = N->Rank = NBlue->Cost + Blue->Pi + N->Pi;
                HeapSiftUp(N);
            }
        }
        // 循环终止条件：所有节点都被放入树中
        while ((NextBlue = HeapDeleteMin())) {
            Follow(NextBlue, Blue);
            Blue = NextBlue;
            //更新blue节点的所有邻居
            for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
                if (!N->Loc)
                    continue;
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                    N->Rank = INT_MIN;
                    HeapSiftUp(N);
                } else if (!Blue->FixedTo2 && !N->FixedTo2 &&
                           (d = NBlue->Cost + Blue->Pi + N->Pi) < N->Cost)
                {
                    N->Dad = Blue;
                    N->Cost = N->Rank = d;
                    HeapSiftUp(N);
                }
            }
        }
    }
     // 这个else只会在程序初始化和结束的时候进入
    else {
        // 这张图是稠密图
        // 初始化所有节点的cost为正无穷
        while ((N = N->Suc) != FirstNode)
            N->Cost = INT_MAX;
        // 循环终止条件:所有节点都出现在最小生成树中
        while ((N = Blue->Suc) != FirstNode) {
            int Min = INT_MAX;
            //更新链表中除blue节点外的所有节点
            do {
                // FixedOrCommon(a,b):判断(a,b)是否在同一条fixed edge上，或者(a,b)这条边是否被一条即将合并的路径包含
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = D(Blue, N);
                    NextBlue = N;
                    Min = INT_MIN;
                } else {
                    if (!Blue->FixedTo2 && !N->FixedTo2 &&
                        !Forbidden(Blue, N) &&
                        (!c || c(Blue, N) < N->Cost) &&
                        (d = D(Blue, N)) < N->Cost) {
                        N->Cost = d;
                        N->Dad = Blue;
                    }
                    if (N->Cost < Min) {
                        Min = N->Cost;
                        NextBlue = N;
                    }
                }
            }
            while ((N = N->Suc) != FirstNode);
            Follow(NextBlue, Blue);
            Blue = NextBlue;
        }
    }
}
