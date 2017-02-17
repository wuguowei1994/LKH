#include "LKH.h"

/*
 * The CreateCandidateSet function determines for each node its set of incident
 * candidate edges.
 *
 * The Ascent function is called to determine a lower bound on the optimal tour
 * using subgradient optimization. But only if the penalties (the Pi-values) is
 * not available on file. In the latter case, the penalties is read from the
 * file, and the lower bound is computed from a minimum 1-tree.
 *
 * The function GenerateCandidates is called to compute the Alpha-values and to
 * associate to each node a set of incident candidate edges.
 *
 * The CreateCandidateSet function itself is called from LKHmain.
 */

/**
 * CreateCandidateSet()函数用来确定每个节点出度候选边
 * Ascent()函数使用次梯度优化算法计算最优解长度的下限
 * GenerateCandidates()函数用来计算每个节点的Alpha值,并把每个节点和一些候选边集合关联起来。
 * LKHmain这个文件会调用CreateCandidateSet()这个函数，接着CreateCandidateSet()会调用Ascent()和GenerateCandidates()函数
 */

void CreateCandidateSet()
{
    GainType Cost, MaxAlpha, A;
    Node *Na;
    int CandidatesRead = 0, i;
    double EntryTime = GetTime();

    Norm = 9999;
    // 这个if语句的作用把所有节点的cost乘以Precision(精度)
    if (C == C_EXPLICIT) {
        Na = FirstNode;
        do {
            for (i = 1; i < Na->Id; i++)

                Na->C[i] *= Precision;
        }
        while ((Na = Na->Suc) != FirstNode);
    }
    // 这个if循环在默认情况下永远不会进入
    if (Distance == Distance_1 ||
            (MaxTrials == 0 &&
             (FirstNode->InitialSuc || InitialTourAlgorithm == SIERPINSKI ||
              InitialTourAlgorithm == MOORE))) {
        ReadCandidates(MaxCandidates);
        AddTourCandidates();
        if (ProblemType == HCP || ProblemType == HPP)
            Ascent();
        goto End_CreateCandidateSet;
    }
    /*
     下面这个if语句不会进入，因为TraceLevel=1
     TRACE_LEVEL = <integer> 指定在计算解的过程中输出的细节的等级(0--1,0表示输出的数量最少;1表示输出的数量最大)。默认值：1
     */
    if (TraceLevel >= 2)
        printff("Creating candidates ...\n");
    /*
     这个if语句也不会进入，它适合已经指定了Pi值的情况
    */
    if (MaxCandidates > 0 &&
            (CandidateSetType == QUADRANT || CandidateSetType == NN)) {
        ReadPenalties();
        if (!(CandidatesRead = ReadCandidates(MaxCandidates)) &&
                MaxCandidates > 0) {
            if (CandidateSetType == QUADRANT)
                CreateQuadrantCandidateSet(MaxCandidates);
            else if (CandidateSetType == NN)
                CreateNearestNeighborCandidateSet(MaxCandidates);
        } else {
            AddTourCandidates();
            if (CandidateSetSymmetric)
                SymmetrizeCandidateSet();
        }
        goto End_CreateCandidateSet;
    }
    /*
    这个循环进入了，因为本例没有事先指定Pi值
     */
    if (!ReadPenalties()) {
        /* No PiFile specified or available */
        Na = FirstNode;
        // int count=0;
        // 一共有2392个节点，这个do while循环会把所有节点的Pi初始化为0
        do
        {
            Na->Pi = 0;
            // count++;
        }
        while ((Na = Na->Suc) != FirstNode);
        // printf("一共有%d个节点\n",count);
        // ReadCandidates()这个函数会从指定文件中读取候选集，成功则返回1，失败则返回0.
        // 由于没有指定，这里函数返回0,CandidatesRead(int类型)仍然是0。所以下面这条语句其实没有任何作用
        CandidatesRead = ReadCandidates(MaxCandidates);
        /*
         Ascent()函数使用次梯度优化算法计算最优解长度的下限。
         该函数运行时会计算每条边的Alpha值(在这里Alpha显示了这条边出现在最优路径中的可能性)
         该函数运行时会计算每个节点的Pi值，使得下限值 L(T(Pi)) - 2*PiSum 最大。
         */
        Cost = Ascent();
        /*
        Subgradient=1
        SubproblemSize=0
        这个if会进入
         */
        if (Subgradient && SubproblemSize == 0) {
            /*
            WritePenalties()函数会把Ascent()函数中计算出来的每个节点的Pi值写入到PiFileName文件中
            由于默认情况下没有指定PiFileName，所以这个函数其实什么都没有做。
             */
            WritePenalties();
            PiFile = 0;
        }
    }
    // 这个else if语句永远不会进入
    else if ((CandidatesRead = ReadCandidates(MaxCandidates)) ||
             MaxCandidates == 0) {
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        goto End_CreateCandidateSet;
    }
    // 这个else语句永远不会进入
    else {
        if (CandidateSetType != DELAUNAY && MaxCandidates > 0) {
            if (TraceLevel >= 2)
                printff("Computing lower bound ... ");
            Cost = Minimum1TreeCost(0);
            if (TraceLevel >= 2)
                printff("done\n");
        } else {
            CreateDelaunayCandidateSet();
            Na = FirstNode;
            do {
                Na->BestPi = Na->Pi;
                Na->Pi = 0;
            }
            while ((Na = Na->Suc) != FirstNode);
            if (TraceLevel >= 2)
                printff("Computing lower bound ... ");
            Cost = Minimum1TreeCost(1);
            if (TraceLevel >= 2)
                printff("done\n");
            Na = FirstNode;
            do {
                Na->Pi = Na->BestPi;
                Cost -= 2 * Na->Pi;
            }
            while ((Na = Na->Suc) != FirstNode);
        }
    }

    LowerBound = (double) Cost / Precision;
    /*
     TraceLevel=1
     会进入这个if语句,但是全部都是打印变量，没有什么实际作用
    */
    if (TraceLevel >= 1) {
        printff("Lower bound = %0.1f", LowerBound);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.2f%%",
                    100.0 * (Optimum - LowerBound) / Optimum);
        if (!PiFile)
            printff(", Ascent time = %0.2f sec.",
                    fabs(GetTime() - EntryTime));
        printff("\n");
    }

    MaxAlpha = (GainType) fabs(Excess * Cost);
    if ((A = Optimum * Precision - Cost) > 0 && A < MaxAlpha)
        MaxAlpha = A;
    //这个if语句永远不会进入
    if (CandidateSetType == DELAUNAY || MaxCandidates == 0)
        OrderCandidateSet(MaxCandidates, MaxAlpha, CandidateSetSymmetric);
    // 会进入这个else语句
    else
        /*
        MaxCandidates:5 节点候选边的最大个数
        MaxAlpha:15614  候选边Alpha的上限
        CandidateSetSymmetric:0  非零表示候选边的集合要被扩充，以至于让每一条候选边都被它的两个端点关联。0表示候选边集合不需要被扩充
        GenerateCandidates()函数会把每个节点和它的入度候选边集合关联起来。每个节点的候选边们会按照候选边的Alpha值排序(升序)
        GenerateCandidates()函数在前面已经被调用过一次，这里已经是第二次调用了
         */
        // 把每个节点和它的入度候选边集合关联起来。这里CandidateSetSymmetric=0,表示候选边集合不需要被扩充
        GenerateCandidates(MaxCandidates, MaxAlpha, CandidateSetSymmetric);
End_CreateCandidateSet:
    // ExtraCandidates:0 这个if不会进入
    if (ExtraCandidates > 0) {
        AddExtraCandidates(ExtraCandidates,
                           ExtraCandidateSetType,
                           ExtraCandidateSetSymmetric);
        AddTourCandidates();
    }
    // ResetCandidates()函数会移除候选边集合中正在使用的边和Alpha等于正无穷的边，然后将集合中的边重新排序
    ResetCandidateSet();
    Na = FirstNode;
    // 下面这个do while循环主要是用来打印的，当没有找到任何候选集的时候，打印错误信息
    // 一般情况下这个do while循环不会进入
    do {

        if (!Na->CandidateSet || !Na->CandidateSet[0].To) {
            if (MaxCandidates == 0)
                eprintf("MAX_CANDIDATES = 0: Node %d has no candidates",
                        Na->Id);
            else
                eprintf("Node %d has no candidates", Na->Id);
        }
    }
    while ((Na = Na->Suc) != FirstNode);
    /*
     CandidatesRead=0
     SubproblemSize=0
     */
    if (!CandidatesRead && SubproblemSize == 0)
    //WriteCandidates()函数会把候选边集合写入到CandidateFileName[0]文件中
    //由于默认情况下不指定输出文件，所以这个函数其实什么都没有做
        WriteCandidates();
    //默认情况下C = C_EXPLICIT 所以这个循环会进入
    if (C == C_EXPLICIT) {
        Na = FirstNode;
        do
            for (i = 1; i < Na->Id; i++)
                //令每个节点的C[i] = C[i]+Pi + NodeSet[i].Pi
                //C[i]是一个记录了权重的矩阵
                //NodeSet是一个一个包含所有节点的数组
                Na->C[i] += Na->Pi + NodeSet[i].Pi;
        while ((Na = Na->Suc) != FirstNode);
    }
    // TraceLevel=1
    if (TraceLevel >= 1) {
        // CandidateReport()函数会把一个节点的候选集中最小、最大和平均数打印出来
        CandidateReport();
        printff("Preprocessing time = %0.2f sec.\n",
                fabs(GetTime() - EntryTime));
    }
}
