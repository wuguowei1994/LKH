#include "Segment.h"
#include "LKH.h"
#include "Hashing.h"
#include "Sequence.h"

/*
  LinKernighan()函数通过opt交换来修正可行解
  这个函数会返回修正后的路径的权重
 */

GainType LinKernighan()
{
    Node *t1, *t2, *SUCt1;
    GainType Gain, G0, Cost;
    int X2, i, it = 0;
    Candidate *Nt1;
    Segment *S;
    SSegment *SS;
    double EntryTime = GetTime();

    Reversed = 0;
    S = FirstSegment;
    i = 0;
    do {
        S->Size = 0;
        S->Rank = ++i;
        S->Reversed = 0;
        S->First = S->Last = 0;
    }
    while ((S = S->Suc) != FirstSegment);
    SS = FirstSSegment;
    i = 0;
    do {
        SS->Size = 0;
        SS->Rank = ++i;
        SS->Reversed = 0;
        SS->First = SS->Last = 0;
    }
    while ((SS = SS->Suc) != FirstSSegment);

    FirstActive = LastActive = 0;
    Swaps = 0;

    /*
    1.计算出初始可行解的权重，并保存在变量Cost中。
    2.计算出对应的hash值，并保存在变量Hash中。
    3.初始化segment(分割)列表
    4.激活所有节点的访问状态
     */
    Cost = 0;
    Hash = 0;
    i = 0;
    t1 = FirstNode;
    do {
        t2 = t1->OldSuc = t1->Suc;
        t1->OldPred = t1->Pred;
        t1->Rank = ++i;
        Cost += (t1->SucCost = t2->PredCost = C(t1, t2)) - t1->Pi - t2->Pi;
        Hash ^= Rand[t1->Id] * Rand[t2->Id];
        t1->Cost = INT_MAX;
        for (Nt1 = t1->CandidateSet; (t2 = Nt1->To); Nt1++)
            if (t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
                t1->Cost = Nt1->Cost;
        t1->Parent = S;
        S->Size++;
        if (S->Size == 1)
            S->First = t1;
        S->Last = t1;
        if (SS->Size == 0)
            SS->First = S;
        S->Parent = SS;
        SS->Last = S;
        if (S->Size == GroupSize) {
            S = S->Suc;
            SS->Size++;
            if (SS->Size == SGroupSize)
                SS = SS->Suc;
        }
        t1->OldPredExcluded = t1->OldSucExcluded = 0;
        t1->Next = 0;
        if (Trial == 1 || KickType == 0 || Kicks == 0 ||
                !InBestTour(t1, t1->Pred) || !InBestTour(t1, t1->Suc))
            Activate(t1);
    }
    while ((t1 = t1->Suc) != FirstNode);
    // 会进入
    if (S->Size < GroupSize)
        SS->Size++;
    Cost /= Precision;
    //不会进入
    if (TraceLevel >= 3 || (TraceLevel == 2 && Cost < BetterCost)) {
        printff("Cost = " GainFormat, Cost);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.4f%%", 100.0 * (Cost - Optimum) / Optimum);
        printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
                Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
    }
    PredSucCostAvailable = 1;

    //循环终止条件:增益小于等于0(增益小于0说明无法找到更好的修正解)
    // do while循环会一直进行opt交换直到找不到新的正增益
    do {
        //取第一个被激活的节点为t1
        while ((t1 = RemoveFirstActive())) {
            //现在t1为非激活状态
            //取t1的下一个节点
            SUCt1 = SUC(t1);
            //不会进入
            if ((TraceLevel >= 3 || (TraceLevel == 2 && Trial == 1)) &&
                    ++it % (Dimension >= 100000 ? 10000 :
                            Dimension >= 10000 ? 1000 : 100) == 0)
                printff("#%d: Time = %0.2f sec.\n",
                        it, fabs(GetTime() - EntryTime));
            //取t2为t1在原始解中的邻居
            for (X2 = 1; X2 <= 2; X2++) {
                t2 = X2 == 1 ? PRED(t1) : SUCt1;
                if (FixedOrCommon(t1, t2) ||
                        (RestrictedSearch && Near(t1, t2) &&
                         (Trial == 1 ||
                          (Trial > BackboneTrials &&
                           (KickType == 0 || Kicks == 0)))))
                    continue;
                G0 = C(t1, t2);
                // 尝试能否找到一条修正解
                do
                    t2 = Swaps == 0 ? BestMove(t1, t2, &G0, &Gain) :
                         BestSubsequentMove(t1, t2, &G0, &Gain);
                while (t2);
                if (Gain > 0) {
                    //如果Gain>0证明修正解被找到
                    assert(Gain % Precision == 0);
                    Cost -= Gain / Precision;
                    //不会进入
                    if (TraceLevel >= 3 ||
                            (TraceLevel == 2 && Cost < BetterCost)) {
                        printff("Cost = " GainFormat, Cost);
                        if (Optimum != MINUS_INFINITY && Optimum != 0)
                            printff(", Gap = %0.4f%%",
                                    100.0 * (Cost - Optimum) / Optimum);
                        printff(", Time = %0.2f sec. %s\n",
                                fabs(GetTime() - EntryTime),
                                Cost < Optimum ? "<" : Cost ==
                                Optimum ? "=" : "");
                    }
                    // 这个函数会存储找到的修正解，令每个节点的OldPred=Pred,OldSuc=Suc,同时更新节点的Cost
                    StoreTour();
                    // HashSearch(HTable,Hash,Cost)：如果HTable中含有Hash和Cost,这个函数会返回1 否则返回0
                    if (HashSearch(HTable, Hash, Cost))
                        goto End_LinKernighan;
                    //重新激活t1节点
                    Activate(t1);
                    break;
                }
                // 如果Gain>0，前面会直接break;如果进入这里说明没有找到大于0的Gain。
                // RestoreTour()函数用来回退opt交换
                RestoreTour();
            }
        }
        // HashSearch(HTable,Hash,Cost)：如果HTable中含有Hash和Cost,这个函数会返回1 否则返回0
        if (HashSearch(HTable, Hash, Cost))
            goto End_LinKernighan;

        // 向哈希表HTable中插入Hash(key)和Cost(value)
        HashInsert(HTable, Hash, Cost);
        /*
          使用4/5-opt交换修正可行解
         */
        Gain = 0;
        if (Gain23Used && (Gain = Gain23()) > 0) {
            // 如果Gain除以Precisio的余数为0,说明找到了修正解。(Gain一般都是100的倍数，Precision=100)
            assert(Gain % Precision == 0);
            Cost -= Gain / Precision;
            // 存储找到的修正解
            StoreTour();
            // TraceLevel == 1，不会进入if语句
            if (TraceLevel >= 3 || (TraceLevel == 2 && Cost < BetterCost)) {
                printff("Cost = " GainFormat, Cost);
                if (Optimum != MINUS_INFINITY && Optimum != 0)
                    printff(", Gap = %0.4f%%",
                            100.0 * (Cost - Optimum) / Optimum);
                printff(", Time = %0.2f sec. + %s\n",
                        fabs(GetTime() - EntryTime),
                        Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
            }
            //不会进入
            if (HashSearch(HTable, Hash, Cost))
                goto End_LinKernighan;
        }
    }
    while (Gain > 0);

End_LinKernighan:
    PredSucCostAvailable = 0;
    // NormalizeNodeList()函数用来交换节点的Suc和Pred，然后得到一个循环的双向链表
    NormalizeNodeList();    
    NormalizeSegmentList();
    return Cost;
}
