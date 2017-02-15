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
 * Ascent()函数用来确定使用次梯度优化算法生成的最优路径的下限(前提条件是输入在文件中无法找到每个节点的Pi值)。
 * 如果在输入文件中已经包含了每个节点的Pi值，那么最优路径的下限会从最小1-tree树中计算出来。
 * GenerateCandidates()函数用来计算每个节点的Alpha值,并把每个节点和一些候选边集合关联起来。
 * LKHmain这个文件会调用CreateCandidateSet()这个函数，其余两个函数会被CreateCandidateSet()调用。
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


        
        Cost = Ascent();
        if (Subgradient && SubproblemSize == 0) {
            WritePenalties();
            PiFile = 0;
        }
    } else if ((CandidatesRead = ReadCandidates(MaxCandidates)) ||
               MaxCandidates == 0) {
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        goto End_CreateCandidateSet;
    } else {
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
    if (CandidateSetType == DELAUNAY || MaxCandidates == 0)
        OrderCandidateSet(MaxCandidates, MaxAlpha, CandidateSetSymmetric);
    else
        GenerateCandidates(MaxCandidates, MaxAlpha, CandidateSetSymmetric);

End_CreateCandidateSet:
    if (ExtraCandidates > 0) {
        AddExtraCandidates(ExtraCandidates,
                           ExtraCandidateSetType,
                           ExtraCandidateSetSymmetric);
        AddTourCandidates();
    }
    ResetCandidateSet();
    Na = FirstNode;
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
    if (!CandidatesRead && SubproblemSize == 0)
        WriteCandidates();
    if (C == C_EXPLICIT) {
        Na = FirstNode;
        do
            for (i = 1; i < Na->Id; i++)
                Na->C[i] += Na->Pi + NodeSet[i].Pi;
        while ((Na = Na->Suc) != FirstNode);
    }
    if (TraceLevel >= 1) {
        CandidateReport();
        printff("Preprocessing time = %0.2f sec.\n",
                fabs(GetTime() - EntryTime));
    }
}
