#include "LKH.h"

/*
 * The GenerateCandidates function associates to each node a set of incident
 * candidate edges. The candidate edges of each node are sorted in increasing
 * order of their Alpha-values.
 *
 * The parameter MaxCandidates specifies the maximum number of candidate edges
 * allowed for each node, and MaxAlpha puts an upper limit on their
 * Alpha-values.
 *
 * A non-zero value of Symmetric specifies that the candidate set is to be
 * complemented such that every candidate edge is associated with both its
 * two end nodes (in this way MaxCandidates may be exceeded).
 *
 * The candidate edges of each node is kept in an array (CandidatSet) of
 * structures. Each structure (Candidate) holds the following information:
 *
 *      Node *To    : pointer to the other end node of the edge
 *      int Cost    : the cost (length) of the edge
 *      int Alpha   : the alpha-value of the edge
 *
 * The algorithm for computing Alpha-values in time O(n^2) and space O(n)
 * follows the description in
 *
 *      Keld Helsgaun,
 *      An Effective Implementation of the Lin-Kernighan Traveling
 *      Salesman Heuristic,
 *      Report, RUC, 1998.
 */
/*
  GenerateCandidates()函数会把每个节点和它的入度候选边集合关联起来。每个节点的候选边们会按照候选边的Alpha值排序(升序)。
  MaxCandidates:节点候选边的最大个数,50
  MaxAlpha:候选边Alpha的上限,3536200
  Symmetric:非零表示候选边的集合要被扩充，以至于让每一条候选边都被它的两个端点关联
            0表示候选边集合不需要被扩充
*/
static int Max(const int a, const int b)
{
    return a > b ? a : b;
}
// 第一次调用时GenerateCandidates(50, 3536200, 1);第二次调用时GenerateCandidates(5,15614,0)
void GenerateCandidates(int MaxCandidates, GainType MaxAlpha,
                        int Symmetric)
{
    Node *From, *To;
    Candidate *NFrom, *NN;
    int a, d, Count;
    //下面这两个if都不会执行,TraceLevel=1
    if (TraceLevel >= 2)
        printff("Generating candidates ... ");
    if (MaxAlpha < 0 || MaxAlpha > INT_MAX)
        MaxAlpha = INT_MAX;
    //初始化每个节点的候选边
    //先清空
    FreeCandidateSets();
    From = FirstNode;
    do
        // #define Mark LastV   在广度优先搜索中标记一个节点
        // 把每一个节点的Mark置为0
        From->Mark = 0;
    while ((From = From->Suc) != FirstNode);
    // MaxCandidates=50
    if (MaxCandidates > 0) {
        do {
            assert(From->CandidateSet =
                       (Candidate *) malloc((MaxCandidates + 1) *
                                            sizeof(Candidate)));
            // 先把所有节点的to指针置为0
            From->CandidateSet[0].To = 0;
        }
        while ((From = From->Suc) != FirstNode);
    }
    // 这个else不会进入
    else {
        AddTourCandidates();
        do {
            if (!From->CandidateSet)
                eprintf("MAX_CANDIDATES = 0: No candidates");
        } while ((From = From->Suc) != FirstNode);
        return;
    }
    /*
    这里是一个双层do while循环，外层从from节点开始，内层从to节点开始。
    循环结束后会把每个节点和它的入度候选边关联起来。每个节点的候选边们会按照候选边的Alpha值排序(升序)。
    */
    /* Loop for each node, From */
    do {
        NFrom = From->CandidateSet;
        // From只会在第一次进来的时候不等于FirstNode
        if (From != FirstNode) {
            // 每个节点的Beta先置为负无穷
            From->Beta = INT_MIN;
            for (To = From; To->Dad != 0; To = To->Dad) {
                /*
                如果节点to和它在最小生成树中的父节点在同一条fixed edge上，或者(to,to->Dad)这条边属于一条将要被合并的路径上的边，
                就把to->Dad->Beta=to->Beta;
                否则就把to->Dad->Beta置为to->Beta和To->Cost中的最大值
                 */
                To->Dad->Beta =
                    !FixedOrCommon(To, To->Dad) ?
                    Max(To->Beta, To->Cost) : To->Beta;
                // 这里其实相当于To->Dad->LastV = From;
                To->Dad->Mark = From;
            }
        }
        Count = 0;
        /* Loop for each node, To */
        To = FirstNode;
        do {
            if (To == From)
                continue;
            /*
             这里问号左边的条件不满足，
             所以会运行D(From, To)，这个函数是C.c文件中的D_EXPLICIT(Node * Na, Node * Nb)
             而c(From, To)是C.c文件中的C_EXPLICIT(Node * Na, Node * Nb)
             D(From, To):先比较From节点和To节点的id，返回的是id较大者的cost[i]加上这两个节点的Pi值,这里i等于id较小者的id
             */
            d = c && !FixedOrCommon(From, To) ? c(From, To) : D(From, To);


            if (From == FirstNode)
                a = To == From->Dad ? 0 : d - From->NextCost;
            else if (To == FirstNode)
                a = From == To->Dad ? 0 : d - To->NextCost;
            else {
                // To->LastV!=From
                if (To->Mark != From)
                    To->Beta =
                        !FixedOrCommon(To, To->Dad) ?
                        Max(To->Dad->Beta, To->Cost) : To->Dad->Beta;
                a = d - To->Beta;
            }
            if (FixedOrCommon(From, To))
                a = INT_MIN;
            else {
                if (From->FixedTo2 || To->FixedTo2 || Forbidden(From, To))
                    continue;
                if (InInputTour(From, To)) {
                    a = 0;
                    if (c)
                        d = D(From, To);
                } else if (c) {
                    if (a > MaxAlpha ||
                            (Count == MaxCandidates &&
                             (a > (NFrom - 1)->Alpha ||
                              (a == (NFrom - 1)->Alpha
                               && d >= (NFrom - 1)->Cost))))
                        continue;
                    if (To == From->Dad) {
                        d = From->Cost;
                        a = 0;
                    } else if (From == To->Dad) {
                        d = To->Cost;
                        a = 0;
                    } else {
                        a -= d;
                        a += (d = D(From, To));
                    }
                }
            }
            // IsPossibleCandidate()用来计算(from,to)这条边是否有可能属于最优解中的一条边。如果这条边有可能是，返回1否则返回0
            if (a <= MaxAlpha && IsPossibleCandidate(From, To)) {
                // 向From->CandidateSet中插入新的候选边
                NN = NFrom;
                while (--NN >= From->CandidateSet) {
                    if (a > NN->Alpha || (a == NN->Alpha && d >= NN->Cost))
                        break;
                    *(NN + 1) = *NN;
                }
                NN++;
                NN->To = To;
                NN->Cost = d;
                NN->Alpha = a;
                if (Count < MaxCandidates) {
                    Count++;
                    NFrom++;
                }
                NFrom->To = 0;
            }
        }
        while ((To = To->Suc) != FirstNode);
    }
    while ((From = From->Suc) != FirstNode);
    // AddTourCandidates()函数会从用户指定的文件中添加候选边集，由于默认情况下不指定，所以AddTourCandidates()函数其实什么都没有做
    AddTourCandidates();
    if (Symmetric)
    // SymmetrizeCandidateSet()函数用来扩充候选集，目的是让每一条候选边都与它两端的节点关联
        SymmetrizeCandidateSet();
    if (TraceLevel >= 2)
        printff("done\n");
}
