#include "LKH.h"

/*
 * The ChooseInitialTour function generates a pseudo-random initial tour.
 * The algorithm constructs a tour as follows.
 *
 * First, a random node N is chosen.
 *
 * Then, as long as no all nodes have been chosen, choose the next node to
 * follow N in the tour, NextN, and set N equal to NextN.
 *
 * NextN is chosen as follows:
 *
 *  (A) If possible, choose NextN such that (N,NextN) is a fixed edge, or
 *      is common to two or more tours to be merged.
 *  (B) Otherwise, if possible, and Trial = 1, choose NextN such that
 *      (N,NextN) is an edge of a given initial tour.
 *  (C) Otherwise, if possible, choose NextN so that (N,NextN) is a
 *      candidate edge, the alpha-value of (N,NextN) is zero, and (N,NextN)
 *      belongs to the current best or next best tour.
 *  (D) Otherwise, if possible, choose NextN such that (N,NextN) is a
 *      candidate edge.
 *  (E) Otherwise, choose NextN at random among those nodes not already
 *      chosen.
 *
 *  When more than one node may be chosen, the node is chosen at random
 *  among the alternatives (a one-way list of nodes).
 *
 *  The sequence of chosen nodes constitutes the initial tour.
 */

/*
  ChooseInitialTour()函数会按照顺序生成一个伪随机的初始路径。
  这个函数使用迭代法生成初始解:
  步骤1：任意选择一个节点作为初始节点N
  步骤2：再选择一个节点NextN作为节点N的后继节点。
  步骤3：让节点N与节点NextN相等
  重复步骤2,3直到所有的节点都被包含进来
  按照如下优先级选择NextN节点
  (A)尽可能让(N,NextN)在一条fixed边上，或者让它同时被两条或两条以上即将合并的tours包含(LKH给的测试文件找不到满足该优先级的点)
  (B)否则，如果Tria=1,就尽可能让(N,NextN)是初始解上的一条边(LKH给的测试文件可以找到满足优先级的点)
  (C)尽可能让(N,NextN)是一条候选边，同时它的Alpha为零，同时这条边属于当前最优路径或者下一条最优路径(可以找到)
  (D)尽可能让(N,NextN)是一条候选边(可以找到)
  (E)只要NextN节点是一个还没有被选择过的节点即可(可以找到)
  这些节点都被选择以后，会组成一个初始解
 */

static int FixedOrCommonCandidates(Node * N);

void ChooseInitialTour()
{
    Node *N, *NextN, *FirstAlternative, *Last;
    Candidate *NN;
    int Alternatives, Count, i;
    //不会进入
    if (KickType > 0 && Kicks > 0 && Trial > 1) {
        for (Last = FirstNode; (N = Last->BestSuc) != FirstNode; Last = N)
            Follow(N, Last);
        for (i = 1; i <= Kicks; i++)
            KSwapKick(KickType);
        return;
    }
    // 无意义
    if (Trial == 1 && (!FirstNode->InitialSuc || InitialTourFraction < 1)) {
        // 不会进入
        if (InitialTourAlgorithm == BORUVKA ||
                InitialTourAlgorithm == GREEDY ||
                InitialTourAlgorithm == MOORE ||
                InitialTourAlgorithm == NEAREST_NEIGHBOR ||
                InitialTourAlgorithm == QUICK_BORUVKA ||
                InitialTourAlgorithm == SIERPINSKI) {
            GainType Cost = InitialTourAlgorithm == MOORE ||
                            InitialTourAlgorithm == SIERPINSKI ?
                            SFCTour(InitialTourAlgorithm) : GreedyTour();
            if (MaxTrials == 0) {
                BetterCost = Cost;
                RecordBetterTour();
            }
            if (!FirstNode->InitialSuc)
                return;
        }
    }

Start:
    /*
    把所有节点的V置为0
    N->V用来标记节点的访问状态:0----未访问    1----已访问
     */
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    Count = 0;

    //从所有节点中任意选择一个点作为FirstNode(要求这个点的fixed或common边小于2)
    do {
        //FixedOrCommonCandidates(N)函数会返回从N节点出发的fixed和common候选边的条数
        if (FixedOrCommonCandidates(N) < 2)
            break;
    }
    while ((N = N->Suc) != FirstNode);
    // 不会进入
    if (ProblemType == ATSP && N->Id <= DimensionSaved)
        FirstNode = N += DimensionSaved;

    //如果一个节点fixed和common候选边的条数等于2,就把这条边移动到FirstNode的前面,方便后面选择节点
    for (Last = FirstNode->Pred; N != Last; N = NextN) {
        NextN = N->Suc;
        if (FixedOrCommonCandidates(N) == 2)
            // 把N节点移动到Last节点的前面
            Follow(N, Last);
    }

    // 把FirstNode标记为已选择
    FirstNode->V = 1;
    N = FirstNode;

    //一直循环直到所有节点都被选中
    while (N->Suc != FirstNode) {
        // 第一个被选择的点
        FirstAlternative = 0;
        // 可以被选择的点
        Alternatives = 0;
        Count++;
        /*
        每个case都以if(Alternatives==0)开头，当有一个case满足时，就会把Alternatives++
        所以这5个case是互相排斥的，只要有一个满足，就不会进入其他的case，那么直接跳到case的末尾，将NextN节点作为N节点的后继
         */
        /*
        Case A：尽可能让(N,NextN)在一条fixed边上，或者让它同时被两条或两条以上即将合并的tours包含(无法找到)
        */
        // 没有意义
        for (NN = N->CandidateSet; (NextN = NN->To); NN++) {
            // 不会进入
            if (!NextN->V && Fixed(N, NextN)) {
                Alternatives++;
                NextN->Next = FirstAlternative;
                FirstAlternative = NextN;
            }
        }
        //不进入
        if (Alternatives == 0 && MergeTourFiles > 1) {

            for (NN = N->CandidateSet; (NextN = NN->To); NN++) {
                if (!NextN->V && IsCommonEdge(N, NextN)) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        /*
        Case B：如果Tria=1(这里的Trial会从1开始不断变大),就尽可能让(N,NextN)是初始解上的一条边(无法找到)
        */
        // 不会进入
        if (Alternatives == 0 && FirstNode->InitialSuc && Trial == 1 &&
                Count <= InitialTourFraction * Dimension) {
            for (NN = N->CandidateSet; (NextN = NN->To); NN++) {
                if (!NextN->V && InInitialTour(N, NextN)) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        /*
        Case C :尽可能让(N,NextN)是一条候选边，它的Alpha为零，同时这条边包含于当前的最优路径或者下一条最优路径(可以找到)
        */
        if (Alternatives == 0 && Trial > 1 &&
                ProblemType != HCP && ProblemType != HPP) {
            for (NN = N->CandidateSet; (NextN = NN->To); NN++) {
                if (!NextN->V && FixedOrCommonCandidates(NextN) < 2 &&
                        NN->Alpha == 0 && (InBestTour(N, NextN) ||
                                           InNextBestTour(N, NextN))) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        /*
        Case D：尽可能让(N,NextN)是一条候选边(可以找到)
        */
        if (Alternatives == 0) {
            for (NN = N->CandidateSet; (NextN = NN->To); NN++) {
                if (!NextN->V && FixedOrCommonCandidates(NextN) < 2) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        /*
        Case E：NextN节点是一个没有被选择过的节点。(满足优先级E的点是伪随机选择的)
        */
        if (Alternatives == 0) {
            NextN = N->Suc;
            //这个while语句不会运行
            while ((FixedOrCommonCandidates(NextN) == 2 || Forbidden(N, NextN))
                    && NextN->Suc != FirstNode)
                NextN = NextN->Suc;
            //这个if语句也不会运行
            if (FixedOrCommonCandidates(NextN) == 2 || Forbidden(N, NextN)) {
                FirstNode = FirstNode->Suc;
                goto Start;
            }
        }
        /*
        满足case X(X=A,B,C,D)的点不止一个，这个else用来从这些满足的点中随机选择一个
         */
        else {
            NextN = FirstAlternative;
            if (Alternatives > 1) {
                i = Random() % Alternatives;
                while (i--)
                    NextN = NextN->Next;
            }
        }
        // 把节点NextN作为节点N的后继节点
        Follow(NextN, N);
        N = NextN;
        N->V = 1;
    }
    /*
    Forbidden()函数用来判断在atsp问题中，(ta,tb)这条边是否是一条权重为无穷的边
    由于在默认情况下，图为tsp，所以Forbidden()函数永远返回0
     */
    // 不会进入
    if (Forbidden(N, N->Suc)) {
        FirstNode = FirstNode->Suc;
        goto Start;
    }
    //不会进入
    if (MaxTrials == 0) {
        GainType Cost = 0;
        N = FirstNode;
        do
            Cost += C(N, N->Suc) - N->Pi - N->Suc->Pi;
        while ((N = N->Suc) != FirstNode);
        Cost /= Precision;
        if (Cost < BetterCost) {
            BetterCost = Cost;
            RecordBetterTour();
        }
    }
}

//FixedOrCommonCandidates(N)函数会返回从N节点出发的fixed和common候选边的条数
static int FixedOrCommonCandidates(Node * N)
{
    int Count = 0;
    Candidate *NN;

    if (N->FixedTo2)
        return 2;
    if (!N->FixedTo1 && MergeTourFiles < 2)
        return 0;
    for (NN = N->CandidateSet; NN->To; NN++)
        if (FixedOrCommon(N, NN->To))
            Count++;
    return Count;
}
