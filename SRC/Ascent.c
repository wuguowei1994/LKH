#include "LKH.h"

/*
  Ascent()函数使用次梯度优化算法计算最优解长度的下限。
  函数运行时会计算每条边的Alpha值(在这里Alpha显示了这条边出现在最优路径中的可能性)
  函数运行时会计算每个节点的Pi值，使得下限值 L(T(Pi)) - 2*PiSum 最大。
  上式中，
  L(T(Pi))表示使用公式"D(i,j)=C(i,j) + Pi(i) + Pi(j)"将每条边的长度转换后最小1-tree树的长度(C(i,j)表示节点i和j之间的距离)
  PiSum表示了所有节点Pi值的和。
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
            // GenerateCandidates()函数会把每个节点和它的入度候选边集合关联起来。
            // 这个函数运行完毕以后，每个节点的候选边们会按照边的Alpha值排序(升序)。
            // AscentCandidates:50
            GenerateCandidates(AscentCandidates, MaxAlpha, 1);
        // 这个else语句永远不会被执行
        else {
            OrderCandidateSet(AscentCandidates, MaxAlpha, 1);
            W = Minimum1TreeCost(1);
            if (!Norm || W / Precision == Optimum)
                return W;
        }
    }
    // ExtraCandidates的默认值为0。
    // 这个if不会进入
    if (ExtraCandidates > 0)
        AddExtraCandidates(ExtraCandidates, ExtraCandidateSetType,
                           ExtraCandidateSetSymmetric);
    if (TraceLevel >= 2) {
        CandidateReport();
        printff("Subgradient optimization ...\n");
    }

    /* Set LastV of every node to V (the node's degree in the 1-tree) */
    //把每一个节点的LastV设置为V,V代表的是每一个节点在1-tree中的度数
    t = FirstNode;
    do
        t->LastV = t->V;
    while ((t = t->Suc) != FirstNode);

    BestW = W0 = W;
    BestNorm = Norm;
    InitialPhase = 1;
    //使用次梯度优化算法,在运行过程中不断减少周期长度(period length)和步长因子(step size)
    for (Period = InitialPeriod, T = InitialStepSize * Precision;
            Period > 0 && T > 0 && Norm != 0; Period /= 2, T /= 2) {
        //在每次迭代的过程中，周期Period和步长因子step都会减半
        //这个if语句永远不会进入
        if (TraceLevel >= 2)
            printff
            ("  T = %d, Period = %d, BestW = %0.1f, Norm = %d\n",
             T, Period, (double) BestW / Precision, Norm);

        for (P = 1; T && P <= Period && Norm != 0; P++) {
            // 调整每个节点的Pi值
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
            /*
              计算出稀疏图中的一颗最小1-tree树
              Minimum1TreeCost()函数会返回最小1-tree的cost
             */ 
            W = Minimum1TreeCost(1);
            //判断是否找到了一条improvement路径
            //在初始情况下W和BestW一样，那么如果经过次梯度优化算法以后生成的1-tree权重变大或者
            //权重不变，但是Norm的值变小
            if (W > BestW || (W == BestW && Norm < BestNorm)) {
                /*
                如果下限值变成了原来的两倍，那么可以认为这是一张稀疏图
                AscentCandidates:50
                Dimension:2392(与节点个数一样)
                下面这个if语句从来没有运行过
                 */
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
                //BestPi用来临时保存次梯度优化算法运行过程中每个节点的Pi值，在次梯度完成以后，会令node->Pi=node->BestPi
                t = FirstNode;
                do
                    t->BestPi = t->Pi;
                while ((t = t->Suc) != FirstNode);
                if (TraceLevel >= 2)
                    printff
                    ("* T = %d, Period = %d, P = %d, BestW = %0.1f, Norm = %d\n",
                     T, Period, P, (double) BestW / Precision, Norm);
                //在初始相位上，step size会被加倍，这个if语句只会在程序的初始阶段运行
                if (InitialPhase && T * sqrt((double) Norm) > 0)
                    T *= 2;
                //如果improvement是在这个周期的最后一次迭代中被发现的，那么将周期加倍
                if (CandidateSetType != DELAUNAY && P == Period
                        && (Period *= 2) > InitialPeriod)
                    Period = InitialPeriod;
            } else {
                if (TraceLevel >= 3)
                    printff
                    ("  T = %d, Period = %d, P = %d, W = %0.1f, Norm = %d\n",
                     T, Period, P, (double) W / Precision, Norm);
                if (InitialPhase && P > Period / 2) {
                    //确定初始相位
                    InitialPhase = 0;
                    P = 0;
                    T = 3 * T / 4;
                }
            }
        }
    }
    //把在次梯度优化算法中得到的BestPi值置为节点最终的Pi值
    t = FirstNode;
    do {
        t->Pi = t->BestPi;
        t->BestPi = 0;
    } while ((t = t->Suc) != FirstNode);

    /*
      计算最小1-tree和它的cost
      MaxCandidates:5
      CandidateSetType!=DELAUNAY
      所以下面的这个Minimum1TreeCost()传入的实参为0
      补充:这里是最后一次调用Minimum1TreeCost()函数
     */
    W = BestW = Minimum1TreeCost(CandidateSetType == DELAUNAY
                                 || MaxCandidates == 0);

    if (MaxCandidates > 0) {
        //FreeCandidateSets()函数会free候选边集合
        FreeCandidateSets();
        //  CandidateSetType!=DELAUNAY
        if (CandidateSetType == DELAUNAY)
            CreateDelaunayCandidateSet();
    } 
    // 这个else不会进入
    else {
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
