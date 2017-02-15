#include "LKH.h"

/*
 * Let T be a minimum spanning tree on the graph, and let N1 be a node of
 * degree one in T. The Connect function determines a shortest edge emanating
 * from N1, but not in T. At return, the Next field of N1 points to the end 
 * node of the edge, and its NextCost field contains the cost of the edge. 
 * However, the search for the shortest edge is stopped if an edge shorter 
 * than a specified threshold (Max) is found.
*/

/*
  假设T是图上的一个最小生成树，N1是这棵树上度数为1的一个节点。Connect()函数就是用来确定一条
  从N1节点出发的最短边，注意这条边不在T中。函数运行结束后，N1的Next指针指向这条最短边的终点，
  NextCost会存储这条最短边的权重。如果在查找过程中发现了一条比临界值(Max)更短的边，函数将结束查找。
 */
void Connect(Node * N1, int Max, int Sparse)
{
    Node *N;
    Candidate *NN1;
    int d;
    N1->Next = 0;
    N1->NextCost = INT_MAX;
    if (!Sparse || N1->CandidateSet == 0 ||
        N1->CandidateSet[0].To == 0 || N1->CandidateSet[1].To == 0) {
        // 在稠密图中寻找需要的边
        N = FirstNode;
        do {
            if (N == N1 || N == N1->Dad || N1 == N->Dad)
                continue;
            if (FixedOrCommon(N1, N)) {
                N1->NextCost = D(N1, N);
                N1->Next = N;
                return;
            }
            if (!N1->FixedTo2 && !N->FixedTo2 &&
                !Forbidden(N1, N) &&
                (!c || c(N1, N) < N1->NextCost) &&
                (d = D(N1, N)) < N1->NextCost) {
                N1->NextCost = d;
                if (d <= Max)
                    return;
                N1->Next = N;
            }
        }
        while ((N = N->Suc) != FirstNode);
    } else {
        //在稀疏图中寻找需要的边
        for (NN1 = N1->CandidateSet; (N = NN1->To); NN1++) {
            if (N == N1->Dad || N1 == N->Dad)
                continue;
            if (FixedOrCommon(N1, N)) {
                N1->NextCost = NN1->Cost + N1->Pi + N->Pi;
                N1->Next = N;
                return;
            }
            if (!N1->FixedTo2 && !N->FixedTo2 &&
                (d = NN1->Cost + N1->Pi + N->Pi) < N1->NextCost) {
                N1->NextCost = d;
                if (d <= Max)
                    return;
                N1->Next = N;
            }
        }
    }
}
