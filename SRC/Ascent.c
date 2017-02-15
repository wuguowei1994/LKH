#include "LKH.h"

/*
 * The Ascent function computes a lower bound on the optimal tour length
 * using subgradient optimization. The function also transforms the original
 * problem into a problem in which the Alpha-values reflect the likelihood
 * of edges being optimal.
 *
 * The function attempts to find penalties (Pi-values) that maximizes the
 * lower bound L(T(Pi)) - 2*PiSum, where L(T(Pi)) denotes the length of the
 * minimum spanning 1-tree computed from the transformed distances, and PiSum
 * denotes the sum of Pi-values. If C(i,j) denotes the length of an edge
 * (i,j), then the transformed distance D(i,j) of an edge is
 * C(i,j) + Pi(i) + Pi(j).
 *
 * The Minimum1TreeCost function is used to compute the cost of a minimum
 * 1-tree.The Generatecandidates function is called in order to generate
 * candidate sets. Minimum 1-trees are then computed in the corresponding
 * sparse graph.
 */
/*
  Ascent()函数计算了使用次梯度优化算法得到的最优解长度的下限。同时，这个函数还会把初始问题转化成一个具有Alpha值的问题。
  在这里，Alpha显示了这个节点出现在最优路径中的可能性。
  这个函数试图去寻找Pi的值，来让下限L(T(Pi)) - 2*PiSum最大化。在这里，L(T(Pi))表示将"距离转化"以后，最小1-tree生成树的长度
  PiSum表示了所有节点Pi值的和。前面说的"距离转化"公式:D(i,j)=C(i,j) + Pi(i) + Pi(j)。(这里C(i,j)表示节点i和j之间的距离)

  Ascent()函数会调用Minimum1TreeCost()函数，这个函数的作用是计算最小1-tree的cost。
  Ascent()函数会调用Generatecandidates()函数，这个函数的作用是生成候选集，每个点的候选集有了以后会得到一个稀疏图，最小1-tree树们就是
  从这里计算出来的。
 */
GainType Ascent()
{
    Node *t;
    GainType BestW, W, W0, Alpha, MaxAlpha = INT_MAX;
    int T, Period, P, InitialPhase, BestNorm;

Start:
    //初始化每个节点的Pi和BestPi为0
    t = FirstNode;
    do
        t->Pi = t->BestPi = 0;
    while ((t = t->Suc) != FirstNode);
    // 默认情况下CandidateSetType = ALPHA，所以下面这个if不会进入
    if (CandidateSetType == DELAUNAY)
        CreateDelaunayCandidateSet();
    // 默认情况下MaxCandidates=5，所以这个else if也不会进入
    else if (MaxCandidates == 0)
        AddTourCandidates();

    //计算最小1-tree的权重
    //第一次运行时，传入的实参为0;根据打印的情况来看，Minimum1TreeCost没有生成一个有效的环路
    W = Minimum1TreeCost(CandidateSetType == DELAUNAY
                         || MaxCandidates == 0);
    // printf("最小1-tree的权重为:%lld\n",W);
    /*
       如果下面的情况有一条满足，那么就会直接返回cost
       (1)用户不希望使用次梯度优化算法
       (2)Minimum1TreeCost()函数在运行过程中Norm为0(Norm为0表示1-tree里面已经包含了一个最优环路)
       很明显，第一次运行时，下面的if不会成立
     */
    if (!Subgradient || !Norm)
        return W;
    /*
      当前数值:
      Optimum=378032;
      Precision=100;
      W=34267000;
      Optimum * Precision - W=3536200
     */
    if (Optimum != MINUS_INFINITY && (Alpha = Optimum * Precision - W) >= 0)
        MaxAlpha = Alpha;
    // if执行以后:MaxAlpha = Alpha=3536200
    // MaxCandidates:每个节点能够关联的候选边的最大值(默认值：5)
    if (MaxCandidates > 0) {
        //给所有的节点生成对称候选集 CandidateSetType默认为ALPHA
        if (CandidateSetType != DELAUNAY)





            // AscentCandidates:50
            GenerateCandidates(AscentCandidates, MaxAlpha, 1);
        else {
            OrderCandidateSet(AscentCandidates, MaxAlpha, 1);
            W = Minimum1TreeCost(1);
            if (!Norm || W / Precision == Optimum)
                return W;
        }
    }
    if (ExtraCandidates > 0)
        AddExtraCandidates(ExtraCandidates, ExtraCandidateSetType,
                           ExtraCandidateSetSymmetric);
    if (TraceLevel >= 2) {
        CandidateReport();
        printff("Subgradient optimization ...\n");
    }

    /* Set LastV of every node to V (the node's degree in the 1-tree) */
    t = FirstNode;
    do
        t->LastV = t->V;
    while ((t = t->Suc) != FirstNode);

    BestW = W0 = W;
    BestNorm = Norm;
    InitialPhase = 1;
    /* Perform subradient optimization with decreasing period length
       and decreasing step size */
    for (Period = InitialPeriod, T = InitialStepSize * Precision;
            Period > 0 && T > 0 && Norm != 0; Period /= 2, T /= 2) {
        /* Period and step size are halved at each iteration */
        if (TraceLevel >= 2)
            printff
            ("  T = %d, Period = %d, BestW = %0.1f, Norm = %d\n",
             T, Period, (double) BestW / Precision, Norm);
        for (P = 1; T && P <= Period && Norm != 0; P++) {
            /* Adjust the Pi-values */
            t = FirstNode;
            do {
                if (t->V != 0) {
                    t->Pi += T * (7 * t->V + 3 * t->LastV) / 10;
                    if (t->Pi > INT_MAX / 4)
                        t->Pi = INT_MAX / 4;
                    else if (t->Pi < -INT_MAX / 4)
                        t->Pi = -INT_MAX / 4;
                }
                t->LastV = t->V;
            }
            while ((t = t->Suc) != FirstNode);
            /* Compute a minimum 1-tree in the sparse graph */
            W = Minimum1TreeCost(1);
            /* Test if an improvement has been found */
            if (W > BestW || (W == BestW && Norm < BestNorm)) {
                /* If the lower bound becomes greater than twice its
                   initial value it is taken as a sign that the graph might be
                   too sparse */
                if (W - W0 > (W0 >= 0 ? W0 : -W0) && AscentCandidates > 0
                        && AscentCandidates < Dimension) {
                    W = Minimum1TreeCost(CandidateSetType == DELAUNAY
                                         || MaxCandidates == 0);
                    if (W < W0) {
                        /* Double the number of candidate edges
                           and start all over again */
                        if (TraceLevel >= 2)
                            printff("Warning: AscentCandidates doubled\n");
                        if ((AscentCandidates *= 2) > Dimension)
                            AscentCandidates = Dimension;
                        goto Start;
                    }
                    W0 = W;
                }
                BestW = W;
                BestNorm = Norm;
                /* Update the BestPi-values */
                t = FirstNode;
                do
                    t->BestPi = t->Pi;
                while ((t = t->Suc) != FirstNode);
                if (TraceLevel >= 2)
                    printff
                    ("* T = %d, Period = %d, P = %d, BestW = %0.1f, Norm = %d\n",
                     T, Period, P, (double) BestW / Precision, Norm);
                /* If in the initial phase, the step size is doubled */
                if (InitialPhase && T * sqrt((double) Norm) > 0)
                    T *= 2;
                /* If the improvement was found at the last iteration of the
                   current period, then double the period */
                if (CandidateSetType != DELAUNAY && P == Period
                        && (Period *= 2) > InitialPeriod)
                    Period = InitialPeriod;
            } else {
                if (TraceLevel >= 3)
                    printff
                    ("  T = %d, Period = %d, P = %d, W = %0.1f, Norm = %d\n",
                     T, Period, P, (double) W / Precision, Norm);
                if (InitialPhase && P > Period / 2) {
                    /* Conclude the initial phase */
                    InitialPhase = 0;
                    P = 0;
                    T = 3 * T / 4;
                }
            }
        }
    }

    t = FirstNode;
    do {
        t->Pi = t->BestPi;
        t->BestPi = 0;
    } while ((t = t->Suc) != FirstNode);

    /* Compute a minimum 1-tree */
    W = BestW = Minimum1TreeCost(CandidateSetType == DELAUNAY
                                 || MaxCandidates == 0);

    if (MaxCandidates > 0) {
        FreeCandidateSets();
        if (CandidateSetType == DELAUNAY)
            CreateDelaunayCandidateSet();
    } else {
        Candidate *Nt;
        t = FirstNode;
        do {
            for (Nt = t->CandidateSet; Nt && Nt->To; Nt++)
                Nt->Cost += t->Pi + Nt->To->Pi;
        }
        while ((t = t->Suc) != FirstNode);
    }
    if (TraceLevel >= 2)
        printff("Ascent: BestW = %0.1f, Norm = %d\n",
                (double) BestW / Precision, Norm);
    return W;
}
