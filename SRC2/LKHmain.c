#include "LKH.h"
#include "Genetic.h"

/*
 * This file contains the main function of the program.
 */

int main(int argc, char *argv[])
{
    // GainType是long long 类型的
    GainType Cost, OldOptimum;
    double Time, LastTime = GetTime();

    /* Read the specification of the problem */
    // 获取输入文件路径
    if (argc >= 2)
        ParameterFileName = argv[1];
    //读取默认参数，包括问题的类型，计算精度等等
    ReadParameters();
    MaxMatrixDimension = 10000;
    // 读取问题
    ReadProblem();
    // SubproblemSize默认值为0,表示不会对原问题进行分割
    if (SubproblemSize > 0) {
        if (DelaunayPartitioning)
            SolveDelaunaySubproblems();
        else if (KarpPartitioning)
            SolveKarpSubproblems();
        else if (KCenterPartitioning)
            SolveKCenterSubproblems();
        else if (KMeansPartitioning)
            SolveKMeansSubproblems();
        else if (RohePartitioning)
            SolveRoheSubproblems();
        else if (MoorePartitioning || SierpinskiPartitioning)
            SolveSFCSubproblems();
        else
            SolveTourSegmentSubproblems();
        return EXIT_SUCCESS;
    }
	// 分配所有除了节点和候选集以外的内存结构
    AllocateStructures();
	// CreateCandidateSet()函数用来确定每个节点出度候选边	
    CreateCandidateSet();
    // 初始化一些统计用的变量
    InitializeStatistics();
    // Norm为178
    if (Norm != 0)
        // Norm!=0说明在前面调用的ascent()函数没有获得最优解，此时就让BestCost为正无穷
        BestCost = PLUS_INFINITY;
    else {
        // 如果进入这个else语句，说明在前面调用的ascent()函数中已经获得了最优解(判断是否获得了最优解的条件就是Norm=0)
        Optimum = BestCost = (GainType) LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        WriteTour(OutputTourFileName, BestTour, BestCost);
        WriteTour(TourFileName, BestTour, BestCost);
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */
    // 默认情况下Runs=10
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        // FindTour()函数会使用LKH算法(里面又调用了opt交换)修正可行解，返回最优解的权重
        Cost = FindTour();    
        // MaxPopulationSize=0
        // 不会进入
        if (MaxPopulationSize > 1) {
            /* Genetic algorithm */
            int i;
            for (i = 0; i < PopulationSize; i++) {
                GainType OldCost = Cost;
                Cost = MergeTourWithIndividual(i);
                if (TraceLevel >= 1 && Cost < OldCost) {
                    printff("  Merged with %d: Cost = " GainFormat, i + 1,
                            Cost);
                    if (Optimum != MINUS_INFINITY && Optimum != 0)
                        printff(", Gap = %0.4f%%",
                                100.0 * (Cost - Optimum) / Optimum);
                    printff("\n");
                }
            }
            if (!HasFitness(Cost)) {
                if (PopulationSize < MaxPopulationSize) {
                    AddToPopulation(Cost);
                    if (TraceLevel >= 1)
                        PrintPopulation();
                } else if (Cost < Fitness[PopulationSize - 1]) {
                    i = ReplacementIndividual(Cost);
                    ReplaceIndividualWithTour(i, Cost);
                    if (TraceLevel >= 1)
                        PrintPopulation();
                }
            }
        } else if (Run > 1)
        // MergeTourWithBestTour()函数会把当前的路径和BestTour[]数组中的储存的路径合并起来得到一个新的路径
            Cost = MergeTourWithBestTour();
        // 进入这个循环证明前面的MergeTourWithBestTour()函数找到了更好的路径
        if (Cost < BestCost) {
            BestCost = Cost;
        // RecordBetterTour()函数会把这个更好的解记录在BetterTour[]数组中.
            RecordBetterTour();
        //RecordBestTour()函数会把当前的最优解记录到BestTour[]数组中
            RecordBestTour();
            WriteTour(OutputTourFileName, BestTour, BestCost);
            WriteTour(TourFileName, BestTour, BestCost);
        }
        OldOptimum = Optimum;
        // 不会进入
        if (Cost < Optimum) {
            if (FirstNode->InputSuc) {
                Node *N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode);
            }
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime);
        // 更新数据
        UpdateStatistics(Cost, Time);
        // 打印
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%",
                        100.0 * (Cost - Optimum) / Optimum);
            printff(", Time = %0.2f sec. %s\n\n", Time,
                    Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
        }
        //不会进入
        if (StopAtOptimum && Cost == OldOptimum && MaxPopulationSize >= 1) {
            Runs = Run;
            break;
        }
        //不会进入
        if (PopulationSize >= 2 &&
            (PopulationSize == MaxPopulationSize ||
             Run >= 2 * MaxPopulationSize) && Run < Runs) {
            Node *N;
            int Parent1, Parent2;
            Parent1 = LinearSelection(PopulationSize, 1.25);
            do
                Parent2 = LinearSelection(PopulationSize, 1.25);
            while (Parent2 == Parent1);
            ApplyCrossover(Parent1, Parent2);
            N = FirstNode;
            do {
                if (ProblemType != HCP && ProblemType != HPP) {
                    int d = C(N, N->Suc);
                    AddCandidate(N, N->Suc, d, INT_MAX);
                    AddCandidate(N->Suc, N, d, INT_MAX);
                }
                N = N->InitialSuc = N->Suc;
            }
            while (N != FirstNode);
        }
        // SRandom(Seed)函数使用给定的seed生成一系列的伪随机数
        SRandom(++Seed);
    }
    PrintStatistics();
    return EXIT_SUCCESS;
}
