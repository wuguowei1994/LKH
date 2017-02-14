#include "LKH.h"
/**
 * 必须指定的参数:
 * PROBLEM_FILE = <string> 指定问题描述文件
 * ASCENT_CANDIDATES = <integer> 指定在ascent阶段，每个节点求几条候选边。默认值 50
 * BACKBONE_TRIALS = <integer> 指定在每次运行时backbone trials的数量。默认值 0
 * BACKTRACKING = { YES | NO } 指定在第一次K-opt交换时是否使用回溯k-opt。默认值 NO
 * CANDIDATE_FILE = <string> 指定要把候选集写在哪个文件里.如果这个文件已经存在，那么直接读取它作为候选边集合。
 * 这个文件的第一行表示了节点的数量。
 * 以后的每一行是这种格式:节点的索引和它在最小生成树中父节点的索引(如果为0证明没有父节点);接着就是从这个节点出发的候选边的数量，然后是具体的候选边索引。
 * 对于每一个候选边，均给出它的终点和alpha值。补充:你可以指定不止一个CANDIDATE_FILE
 * CANDIDATE_SET_TYPE = { ALPHA | DELAUNAY [ PURE ] | NEAREST-NEIGHBOR | QUADRANT } 指定candidate set的类型。默认值:ALPHA
 * EOF 终止输入数据
 * (说明:real数据类型用于存储单精度浮点数)
 * EXCESS = <real> 每一个候选边的最大alpha值。默认值: 1.0/DIMENSION
 * EXTRA_CANDIDATES = <integer> [ SYMMETRIC ] : 每一个节点的额外候选边集合(它的候选类型需要另外指定) 默认值 0
 * integer的后面可以跟一个关键词 SYMMETRIC ,表示这些候选边需要补全来让每一条边都和它的两个终点关联起来。
 * EXTRA_CANDIDATE_SET_TYPE = { NEAREST-NEIGHBOR | QUADRANT } 额外候选边集合的类型。默认值:QUADRANT
 * GAIN23 = { YES | NO } 是否使用Gain23函数。默认值:YES
 * GAIN_CRITERION = { YES | NO }  是否使用LKH的增益规范。默认值:YES
 * INITIAL_PERIOD = <integer> 在ascent中第一周期的长度。默认值:DIMENSION/2 (不得小于100)
 * INITIAL_STEP_SIZE = <integer> 在ascent中初始步长因子。默认值:1
 * INITIAL_TOUR_ALGORITHM = { BORUVKA | GREEDY | MOORE | NEAREST-NEIGHBOR | QUICK-BORUVKA | SIERPINSKI | WALK } 指定获得初始解的算法。默认值:WALK
 * INITIAL_TOUR_FILE = <string> 指定包含初始解的文件名，这个文件里的初始解会作为LKH算法的输入。这个初始解是一个节点的列表,这个文件以-1作为结束标志
 * INITIAL_TOUR_FRACTION = <real in [0;1]> 指定创建初始解用到的小数 默认值:1.0
 * INPUT_TOUR_FILE = <string> 指定含有tour的输入文件。这条路径是用来限制搜索的，另外，每条边的Alpha会被设置为0。和前面一样，这个tour也是节点的list，以-1作为结束标志
 * KICK_TYPE = <integer> 指定k-swap中kick(踢，反冲)的值。如果KICK_TYPE 为0,那么LKH就会使用特殊的kicking策略-----WALK.默认值：0
 * KICKS = <integer> 指定kick LKH路径的次数。如果KICKS为0,那么就是用WALK策略。默认值：1
 * MAX_BREADTH = <integer> 在LKH每一次move时，考虑的候选边的最大值。默认值：INT_MAX(正无穷)
 * MAX_CANDIDATES = <integer> [ SYMMETRIC ] 每个节点能够关联的候选边的最大值 。 默认值：5
 * MAX_SWAPS = <integer>指定在每一次路径优化过程中，swaps的最大值。默认值：与DIMENSION相等
 * MAX_TRIALS = <integer> 在每次运行时trials的最大值。默认值：节点的数量(也就是problem file中的DIMENSION)
 * MERGE_TOUR_FILE = <string> 指定要被合并的路径名。这个路径中的边将会被添加到candidate sets中。MERGE_TOUR_FILE 可以指定多个。
 * MOVE_TYPE = <integer> 指定LKH中sequential move的类型。默认值：5
 * NONSEQUENTIAL_MOVE_TYPE = <integer> 指定LKH中nonsequential move的类型。默认值: (MOVE_TYPE + PATCHING_C + PATCHING_A - 1).
 * OUTPUT_TOUR_FILE = <string> 指定最优解的输出文件路径
 * OPTIMUM = <integer> 已知的最优路径长度。如果STOP_AT_OPTIMUM为YES，那么程序会在路径长度等于OPTIMUM的情况下停止。默认值：负无穷
 * PATCHING_A = <integer> [ RESTRICTED | EXTENDED ] 在修正路径的过程中，不相交的交换环路的最大值。默认值：1
 * PATCHING_C = <integer> [ RESTRICTED | EXTENDED ] 在寻找可行、有增益的move过程中，不相交的环路的最大值。默认值：0
 * PI_FILE = <string> 指定Pi-values的输出路径。这个文件的第一行是节点的数量，接下来的每一行都是这种格式: <integer> <integer>。第一个integer是节点的索引，第二个integer是这个节点的Pi-value。
 * 如果你指定	file name为0表示所有节点的Pi是0
 * POPULATION_SIZE = <integer> 指定遗传算法中人口的最大值。默认值：0
 * PRECISION = <integer> 在转换距离( d[i][j] = PRECISION*c[i][j] + pi[i] + pi[j])的过程中，使用的精度。在刚才的公式中,d[i][j],c[i][j],pi[i]和pi[j]都是整数。默认值:100
 * RESTRICTED_SEARCH = { YES | NO } 指定优化的技巧：如果已经有了一个解，那么第一条被删除的边必须不属于当前解。如果当前没有解，那么第一条被删除的边必须不属于最小1-tree生成树。默认值：YES
 * RUNS = <integer> 运行的总次数：10
 * SEED = <integer> 指定遗传算法中的初始seed。默认值:1
 * STOP_AT_OPTIMUM = { YES | NO } 当路径长度和OPTIMUM相等时，终止计算。默认值:YES
 * SUBGRADIENT = { YES | NO } 指定是否使用次梯度优化算法计算Pi的值。默认值：YES
 * SUBPROBLEM_SIZE = <integer> [ DELAUNAY | KARP | K-CENTER | K-MEANS | MOORE |
 *                               ROHE | SIERPINSKI ] [ BORDERS ] [ COMPRESSED ]
 * 在原始问题中被分配到subproblems中的节点数量 。默认值:0(表示不会对原始问题进行分割)
 * SUBPROBLEM_TOUR_FILE = <string> 由于前面SUBPROBLEM_SIZE为0,这个不会被用到
 * SUBSEQUENT_MOVE_TYPE = <integer> 默认值：0;表示所有的opt-move都是同类型的。
 * SUBSEQUENT_PATCHING = { YES | NO } 指定在第一次move以后是否需要修正路径。默认值：YES
 * TIME_LIMIT = <real> 指定每次运行的时间限制(单位:秒)。默认值：DBL_MAX
 * TOUR_FILE = <string> 指定在运行过程中每次有新的最优解被找到时，best tour的输出路径
 * TRACE_LEVEL = <integer> 指定在计算解的过程中输出的细节的等级(0--1,0表示输出的数量最少;1表示输出的数量最大)。默认值：1
 *
 */
/**
 * [ReadParameters description]
 * 这个函数会读取输入文件，并设置所有LKH需要用到的变量为默认值
 */
static char Delimiters[] = "= \n\t\r\f\v\xef\xbb\xbf";
static char *GetFileName(char *Line);
static char *ReadYesOrNo(int *V);
#undef max
static size_t max(size_t a, size_t b);
void ReadParameters()
{
	// 把上面提到的变量设置为默认值
	char *Line, *Keyword, *Token, *Name;
	unsigned int i;
	ProblemFileName = PiFileName = InputTourFileName =
	                                   OutputTourFileName = TourFileName = 0;
	CandidateFiles = MergeTourFiles = 0;
	AscentCandidates = 50;
	BackboneTrials = 0;
	Backtracking = 0;
	CandidateSetSymmetric = 0;
	CandidateSetType = ALPHA;
	Crossover = ERXT;
	DelaunayPartitioning = 0;
	DelaunayPure = 0;
	Excess = -1;
	ExtraCandidates = 0;
	ExtraCandidateSetSymmetric = 0;
	ExtraCandidateSetType = QUADRANT;
	Gain23Used = 1;
	GainCriterionUsed = 1;
	InitialPeriod = -1;
	InitialStepSize = 0;
	InitialTourAlgorithm = WALK;
	InitialTourFraction = 1.0;
	KarpPartitioning = 0;
	KCenterPartitioning = 0;
	KMeansPartitioning = 0;
	Kicks = 1;
	KickType = 0;
	MaxBreadth = INT_MAX;
	MaxCandidates = 5;
	MaxPopulationSize = 0;
	MaxSwaps = -1;
	MaxTrials = -1;
	MoorePartitioning = 0;
	MoveType = 5;
	NonsequentialMoveType = -1;
	Optimum = MINUS_INFINITY;
	PatchingA = 1;
	PatchingC = 0;
	PatchingAExtended = 0;
	PatchingARestricted = 0;
	PatchingCExtended = 0;
	PatchingCRestricted = 0;
	Precision = 100;
	RestrictedSearch = 1;
	RohePartitioning = 0;
	Runs = 0;
	Seed = 1;
	SierpinskiPartitioning = 0;
	StopAtOptimum = 1;
	Subgradient = 1;
	SubproblemBorders = 0;
	SubproblemsCompressed = 0;
	SubproblemSize = 0;
	SubsequentMoveType = 0;
	SubsequentPatching = 1;
	TimeLimit = DBL_MAX;
	TraceLevel = 1;

	/**
	 * 读取输入文件
	 */
	if (ParameterFileName) {
		if (!(ParameterFile = fopen(ParameterFileName, "r")))
			eprintf("Cannot open PARAMETER_FILE: \"%s\"",
			        ParameterFileName);
		printff("PARAMETER_FILE = %s\n", ParameterFileName);
	} else {
		while (1) {
			printff("PARAMETER_FILE = ");
			if (!(ParameterFileName = GetFileName(ReadLine(stdin)))) {
				do {
					printff("PROBLEM_FILE = ");
					ProblemFileName = GetFileName(ReadLine(stdin));
				} while (!ProblemFileName);
				return;
			} else if (!(ParameterFile = fopen(ParameterFileName, "r")))
				printff("Cannot open \"%s\". Please try again.\n",
				        ParameterFileName);
			else
				break;
		}
	}
	while ((Line = ReadLine(ParameterFile))) {
		if (!(Keyword = strtok(Line, Delimiters)))
			continue;
		if (Keyword[0] == '#')
			continue;
		for (i = 0; i < strlen(Keyword); i++)
			Keyword[i] = (char) toupper(Keyword[i]);
		if (!strcmp(Keyword, "ASCENT_CANDIDATES")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &AscentCandidates))
				eprintf("ASCENT_CANDIDATES: integer expected");
			if (AscentCandidates < 2)
				eprintf("ASCENT_CANDIDATES: >= 2 expected");
		} else if (!strcmp(Keyword, "BACKBONE_TRIALS")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &BackboneTrials))
				eprintf("BACKBONE_TRIALS: integer expected");
			if (BackboneTrials < 0)
				eprintf("BACKBONE_TRIALS: non-negative integer expected");
		} else if (!strcmp(Keyword, "BACKTRACKING")) {
			if (!ReadYesOrNo(&Backtracking))
				eprintf("BACKTRACKING: YES or NO expected");
		} else if (!strcmp(Keyword, "CANDIDATE_FILE")) {
			if (!(Name = GetFileName(0)))
				eprintf("CANDIDATE_FILE: string expected");
			if (CandidateFiles == 0) {
				assert(CandidateFileName =
				           (char **) malloc(sizeof(char *)));
				CandidateFileName[CandidateFiles++] = Name;
			} else {
				int i;
				for (i = 0; i < CandidateFiles; i++)
					if (!strcmp(Name, CandidateFileName[i]))
						break;
				if (i == CandidateFiles) {
					assert(CandidateFileName =
					           (char **) realloc(CandidateFileName,
					                             (CandidateFiles +
					                              1) * sizeof(char *)));
					CandidateFileName[CandidateFiles++] = Name;
				}
			}
		} else if (!strcmp(Keyword, "CANDIDATE_SET_TYPE")) {
			if (!(Token = strtok(0, Delimiters)))
				eprintf("%s", "CANDIDATE_SET_TYPE: "
				        "ALPHA, DELAUNAY, NEAREST-NEIGHBOR, "
				        "or QUADRANT expected");
			for (i = 0; i < strlen(Token); i++)
				Token[i] = (char) toupper(Token[i]);
			if (!strncmp(Token, "ALPHA", strlen(Token)))
				CandidateSetType = ALPHA;
			else if (!strncmp(Token, "DELAUNAY", strlen(Token))) {
				CandidateSetType = DELAUNAY;
			} else if (!strncmp(Token, "NEAREST-NEIGHBOR", strlen(Token)))
				CandidateSetType = NN;
			else if (!strncmp(Token, "QUADRANT", strlen(Token)))
				CandidateSetType = QUADRANT;
			else
				eprintf("%s", "CANDIDATE_SET_TYPE: "
				        "ALPHA, DELAUNAY, NEAREST-NEIGHBOR, "
				        "or QUADRANT expected");
			if (CandidateSetType == DELAUNAY) {
				if ((Token = strtok(0, Delimiters))) {
					for (i = 0; i < strlen(Token); i++)
						Token[i] = (char) toupper(Token[i]);
					if (strncmp(Token, "PURE", strlen(Token)))
						eprintf("%s", "CANDIDATE_SET_TYPE (DELAUNAY): "
						        "PURE or no token expected");
					DelaunayPure = 1;
				}
			}
		} else if (!strcmp(Keyword, "COMMENT"))
			continue;
		else if (!strcmp(Keyword, "EOF"))
			break;
		else if (!strcmp(Keyword, "EXCESS")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%lf", &Excess))
				eprintf("EXCESS: real expected");
			if (Excess < 0)
				eprintf("EXCESS: non-negeative real expected");
		} else if (!strcmp(Keyword, "EXTRA_CANDIDATES")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &ExtraCandidates))
				eprintf("EXTRA_CANDIDATES: integer expected");
			if (ExtraCandidates < 0)
				eprintf("EXTRA_CANDIDATES: non-negative integer expected");
			if ((Token = strtok(0, Delimiters))) {
				for (i = 0; i < strlen(Token); i++)
					Token[i] = (char) toupper(Token[i]);
				if (strncmp(Token, "SYMMETRIC", strlen(Token)))
					eprintf
					("(EXTRA_CANDIDATES) Illegal SYMMETRIC specification");
				ExtraCandidateSetSymmetric = 1;
			}
		} else if (!strcmp(Keyword, "EXTRA_CANDIDATE_SET_TYPE")) {
			if (!(Token = strtok(0, Delimiters)))
				eprintf("%s", "EXTRA_CANDIDATE_SET_TYPE: "
				        "NEAREST-NEIGHBOR, or QUADRANT expected");
			for (i = 0; i < strlen(Token); i++)
				Token[i] = (char) toupper(Token[i]);
			if (!strncmp(Token, "NEAREST-NEIGHBOR", strlen(Token)))
				ExtraCandidateSetType = NN;
			else if (!strncmp(Token, "QUADRANT", strlen(Token)))
				ExtraCandidateSetType = QUADRANT;
			else
				eprintf("%s", "EXTRA_CANDIDATE_SET_TYPE: "
				        "NEAREST-NEIGHBOR or QUADRANT expected");
		} else if (!strcmp(Keyword, "GAIN23")) {
			if (!ReadYesOrNo(&Gain23Used))
				eprintf("GAIN23: YES or NO expected");
		} else if (!strcmp(Keyword, "GAIN_CRITERION")) {
			if (!ReadYesOrNo(&GainCriterionUsed))
				eprintf("GAIN_CRITERION: YES or NO expected");
		} else if (!strcmp(Keyword, "INITIAL_PERIOD")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &InitialPeriod))
				eprintf("INITIAL_PERIOD: integer expected");
			if (InitialPeriod < 0)
				eprintf("INITIAL_PERIOD: non-negative integer expected");
		} else if (!strcmp(Keyword, "INITIAL_STEP_SIZE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &InitialStepSize))
				eprintf("INITIAL_STEP_SIZE: integer expected");
			if (InitialStepSize <= 0)
				eprintf("INITIAL_STEP_SIZE: positive integer expected");
		} else if (!strcmp(Keyword, "INITIAL_TOUR_ALGORITHM")) {
			if (!(Token = strtok(0, Delimiters)))
				eprintf("INITIAL_TOUR_ALGORITHM: "
				        "BORUVKA, GREEDY, MOORE, NEAREST-NEIGHBOR,\n"
				        "QUICK-BORUVKA, SIERPINSKI, or WALK expected");
			for (i = 0; i < strlen(Token); i++)
				Token[i] = (char) toupper(Token[i]);
			if (!strncmp(Token, "BORUVKA", strlen(Token)))
				InitialTourAlgorithm = BORUVKA;
			else if (!strncmp(Token, "GREEDY", strlen(Token)))
				InitialTourAlgorithm = GREEDY;
			else if (!strncmp(Token, "MOORE", strlen(Token)))
				InitialTourAlgorithm = MOORE;
			else if (!strncmp(Token, "NEAREST-NEIGHBOR", strlen(Token)))
				InitialTourAlgorithm = NEAREST_NEIGHBOR;
			else if (!strncmp(Token, "QUICK-BORUVKA", strlen(Token)))
				InitialTourAlgorithm = QUICK_BORUVKA;
			else if (!strncmp(Token, "SIERPINSKI", strlen(Token)))
				InitialTourAlgorithm = SIERPINSKI;
			else if (!strncmp(Token, "WALK", strlen(Token)))
				InitialTourAlgorithm = WALK;
			else
				eprintf("INITIAL_TOUR_ALGORITHM: "
				        "BORUVKA, GREEDY, MOORE, NEAREST-NEIGHBOR,\n"
				        "QUICK-BORUVKA, SIERPINSKI or WALK expected");
		} else if (!strcmp(Keyword, "INITIAL_TOUR_FILE")) {
			if (!(InitialTourFileName = GetFileName(0)))
				eprintf("INITIAL_TOUR_FILE: string expected");
		} else if (!strcmp(Keyword, "INITIAL_TOUR_FRACTION")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%lf", &InitialTourFraction))
				eprintf("INITIAL_TOUR_FRACTION: real expected");
			if (InitialTourFraction < 0 || InitialTourFraction > 1)
				eprintf("INITIAL_TOUR_FRACTION: >= 0 or <= 1 expected");
		} else if (!strcmp(Keyword, "INPUT_TOUR_FILE")) {
			if (!(InputTourFileName = GetFileName(0)))
				eprintf("INPUT_TOUR_FILE: string expected");
		} else if (!strcmp(Keyword, "KICK_TYPE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &KickType))
				eprintf("KICK_TYPE: integer expected");
			if (KickType != 0 && KickType < 4)
				eprintf("KICK_TYPE: integer >= 4 expected");
		} else if (!strcmp(Keyword, "KICKS")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &Kicks))
				eprintf("KICKS: integer expected");
			if (Kicks < 0)
				eprintf("KICKS: non-negative integer expected");
		} else if (!strcmp(Keyword, "MAX_BREADTH")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &MaxBreadth))
				eprintf("MAX_BREADTH: integer expected");
			if (MaxBreadth < 0)
				eprintf("MAX_BREADTH: non-negative integer expected");
		} else if (!strcmp(Keyword, "MAX_CANDIDATES")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &MaxCandidates))
				eprintf("MAX_CANDIDATES: integer expected");
			if (MaxCandidates < 0)
				eprintf("MAX_CANDIDATES: non-negative integer expected");
			if ((Token = strtok(0, Delimiters))) {
				for (i = 0; i < strlen(Token); i++)
					Token[i] = (char) toupper(Token[i]);
				if (!strncmp(Token, "SYMMETRIC", strlen(Token)))
					CandidateSetSymmetric = 1;
				else
					eprintf
					("(MAX_CANDIDATES) Illegal SYMMETRIC specification");
			}
		} else if (!strcmp(Keyword, "MAX_SWAPS")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &MaxSwaps))
				eprintf("MAX_SWAPS: integer expected");
			if (MaxSwaps < 0)
				eprintf("MAX_SWAPS: non-negative integer expected");
		} else if (!strcmp(Keyword, "MAX_TRIALS")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &MaxTrials))
				eprintf("MAX_TRIALS: integer expected");
			if (MaxTrials < 0)
				eprintf("MAX_TRIALS: non-negative integer expected");
		} else if (!strcmp(Keyword, "MERGE_TOUR_FILE")) {
			if (!(Name = GetFileName(0)))
				eprintf("MERGE_TOUR_FILE: string expected");
			if (MergeTourFiles == 0) {
				assert(MergeTourFileName =
				           (char **) malloc(sizeof(char *)));
				MergeTourFileName[MergeTourFiles++] = Name;
			} else {
				int i;
				for (i = 0; i < MergeTourFiles; i++)
					if (!strcmp(Name, MergeTourFileName[i]))
						break;
				if (i == MergeTourFiles) {
					assert(MergeTourFileName =
					           (char **) realloc(MergeTourFileName,
					                             (MergeTourFiles +
					                              1) * sizeof(char *)));
					MergeTourFileName[MergeTourFiles++] = Name;
				}
			}
		} else if (!strcmp(Keyword, "MOVE_TYPE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &MoveType))
				eprintf("MOVE_TYPE: integer expected");
			if (MoveType < 2)
				eprintf("MOVE_TYPE: >= 2 expected");
		} else if (!strcmp(Keyword, "NONSEQUENTIAL_MOVE_TYPE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &NonsequentialMoveType))
				eprintf("NONSEQUENTIAL_MOVE_TYPE: integer expected");
			if (NonsequentialMoveType < 4)
				eprintf("NONSEQUENTIAL_MOVE_TYPE: >= 4 expected");
		} else if (!strcmp(Keyword, "OPTIMUM")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, GainInputFormat, &Optimum))
				eprintf("OPTIMUM: integer expected");
		} else if (!strcmp(Keyword, "OUTPUT_TOUR_FILE")) {
			if (!(OutputTourFileName = GetFileName(0)))
				eprintf("OUTPUT_TOUR_FILE: string expected");
		} else if (!strcmp(Keyword, "PATCHING_A")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &PatchingA))
				eprintf("PATCHING_A: integer expected");
			if (PatchingA < 0)
				eprintf("PATCHING_A: non-negative integer expected");
			if ((Token = strtok(0, Delimiters))) {
				for (i = 0; i < strlen(Token); i++)
					Token[i] = (char) toupper(Token[i]);
				if (!strncmp(Token, "RESTRICTED", strlen(Token)))
					PatchingARestricted = 1;
				else if (!strncmp(Token, "EXTENDED", strlen(Token)))
					PatchingAExtended = 1;
				else
					eprintf("%s", "(PATCHING_A) "
					        "Illegal RESTRICTED or EXTENDED specification");
			}
		} else if (!strcmp(Keyword, "PATCHING_C")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &PatchingC))
				eprintf("PATCHING_C: integer expected");
			if (PatchingC < 0)
				eprintf("PATCHING_C: non-negative integer expected");
			if ((Token = strtok(0, Delimiters))) {
				for (i = 0; i < strlen(Token); i++)
					Token[i] = (char) toupper(Token[i]);
				if (!strncmp(Token, "RESTRICTED", strlen(Token)))
					PatchingCRestricted = 1;
				else if (!strncmp(Token, "EXTENDED", strlen(Token)))
					PatchingCExtended = 1;
				else
					eprintf("%s", "(PATCHING_C) ",
					        "Illegal RESTRICTED or EXTENDED specification");
			}
		} else if (!strcmp(Keyword, "PI_FILE")) {
			if (!(PiFileName = GetFileName(0)))
				eprintf("PI_FILE: string expected");
		} else if (!strcmp(Keyword, "POPULATION_SIZE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &MaxPopulationSize))
				eprintf("POPULATION_SIZE: integer expected");
		} else if (!strcmp(Keyword, "PRECISION")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &Precision))
				eprintf("PRECISION: integer expected");
		} else if (!strcmp(Keyword, "PROBLEM_FILE")) {
			if (!(ProblemFileName = GetFileName(0)))
				eprintf("PROBLEM_FILE: string expected");
		} else if (!strcmp(Keyword, "RESTRICTED_SEARCH")) {
			if (!ReadYesOrNo(&RestrictedSearch))
				eprintf("RESTRICTED_SEARCH: YES or NO expected");
		} else if (!strcmp(Keyword, "RUNS")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &Runs))
				eprintf("RUNS: integer expected");
			if (Runs <= 0)
				eprintf("RUNS: positive integer expected");
		} else if (!strcmp(Keyword, "SEED")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%u", &Seed))
				eprintf("SEED: integer expected");
		} else if (!strcmp(Keyword, "STOP_AT_OPTIMUM")) {
			if (!ReadYesOrNo(&StopAtOptimum))
				eprintf("STOP_AT_OPTIMUM: YES or NO expected");
		} else if (!strcmp(Keyword, "SUBGRADIENT")) {
			if (!ReadYesOrNo(&Subgradient))
				eprintf("SUBGRADIENT: YES or NO expected");
		} else if (!strcmp(Keyword, "SUBPROBLEM_TOUR_FILE")) {
			if (!(SubproblemTourFileName = GetFileName(0)))
				eprintf("SUBPROBLEM_TOUR_FILE: string expected");
		} else if (!strcmp(Keyword, "SUBPROBLEM_SIZE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &SubproblemSize))
				eprintf("SUBPROBLEM_SIZE: integer expected");
			if (SubproblemSize < 3)
				eprintf("SUBPROBLEM_SIZE: >= 3 expected");
			if ((Token = strtok(0, Delimiters))) {
				for (i = 0; i < strlen(Token); i++)
					Token[i] = (char) toupper(Token[i]);
				if (!strncmp(Token, "DELAUNAY", strlen(Token)))
					DelaunayPartitioning = 1;
				else if (!strncmp(Token, "KARP", max(strlen(Token), 2)))
					KarpPartitioning = 1;
				else if (!strncmp(Token, "K-CENTER", max(strlen(Token), 3)))
					KCenterPartitioning = 1;
				else if (!strncmp(Token, "K-MEANS", max(strlen(Token), 3)))
					KMeansPartitioning = 1;
				else if (!strncmp(Token, "MOORE", strlen(Token)))
					MoorePartitioning = 1;
				else if (!strncmp(Token, "ROHE", strlen(Token)))
					RohePartitioning = 1;
				else if (!strncmp(Token, "SIERPINSKI", strlen(Token)))
					SierpinskiPartitioning = 1;
				else if (!strncmp(Token, "BORDERS", strlen(Token)))
					SubproblemBorders = 1;
				else if (!strncmp(Token, "COMPRESSED", strlen(Token)))
					SubproblemsCompressed = 1;
				else
					eprintf
					("(SUBPROBLEM_SIZE) Illegal DELAUNAY, KARP, K-CENTER, "
					 "K-MEANS, MOORE, ROHE,\n SIERPINSKI, "
					 "BORDERS or COMPRESSED specification");
				while ((Token = strtok(0, Delimiters))) {
					for (i = 0; i < strlen(Token); i++)
						Token[i] = (char) toupper(Token[i]);
					if (!strncmp(Token, "BORDERS", strlen(Token)))
						SubproblemBorders = 1;
					else if (!strncmp(Token, "COMPRESSED", strlen(Token)))
						SubproblemsCompressed = 1;
					else
						eprintf
						("(SUBPROBLEM_SIZE) Illegal BORDERS or "
						 "COMPRESSED specification");
				}
			}
		} else if (!strcmp(Keyword, "SUBSEQUENT_MOVE_TYPE")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &SubsequentMoveType))
				eprintf("SUBSEQUENT_MOVE_TYPE: integer expected");
			if (SubsequentMoveType != 0 && SubsequentMoveType < 2)
				eprintf("SUBSEQUENT_MOVE_TYPE: 0 or >= 2 expected");
		} else if (!strcmp(Keyword, "SUBSEQUENT_PATCHING")) {
			if (!ReadYesOrNo(&SubsequentPatching))
				eprintf("SUBSEQUENT_PATCHING: YES or NO expected");
		} else if (!strcmp(Keyword, "TIME_LIMIT")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%lf", &TimeLimit))
				eprintf("TIME_LIMIT: real expected");
			if (TimeLimit < 0)
				eprintf("TIME_LIMIT: >= 0 expected");
		} else if (!strcmp(Keyword, "TOUR_FILE")) {
			if (!(TourFileName = GetFileName(0)))
				eprintf("TOUR_FILE: string expected");
		} else if (!strcmp(Keyword, "TRACE_LEVEL")) {
			if (!(Token = strtok(0, Delimiters)) ||
			        !sscanf(Token, "%d", &TraceLevel))
				eprintf("TRACE_LEVEL: integer expected");
		} else
			eprintf("Unknown keyword: %s", Keyword);
		if ((Token = strtok(0, Delimiters)) && Token[0] != '#')
			eprintf("Junk at end of line: %s", Token);
	}
	if (!ProblemFileName)
		eprintf("Problem file name is missing");
	if (SubproblemSize == 0 && SubproblemTourFileName != 0)
		eprintf("SUBPROBLEM_SIZE specification is missing");
	if (SubproblemSize > 0 && SubproblemTourFileName == 0)
		eprintf("SUBPROBLEM_TOUR_FILE specification is missing");
	fclose(ParameterFile);
	free(LastLine);
	LastLine = 0;
}

static char *GetFileName(char *Line)
{
	char *Rest = strtok(Line, "\n\t\r\f"), *t;

	if (!Rest)
		return 0;
	while (isspace(*Rest))
		Rest++;
	if (!Line) {
		if (*Rest == '=')
			Rest++;
	}
	while (isspace(*Rest))
		Rest++;
	for (t = Rest + strlen(Rest) - 1; isspace(*t); t--)
		*t = '\0';
	if (!strlen(Rest))
		return 0;
	assert(t = (char *) malloc(strlen(Rest) + 1));
	strcpy(t, Rest);
	return t;
}

static char *ReadYesOrNo(int *V)
{
	char *Token = strtok(0, Delimiters);

	if (Token) {
		unsigned int i;
		for (i = 0; i < strlen(Token); i++)
			Token[i] = (char) toupper(Token[i]);
		if (!strncmp(Token, "YES", strlen(Token)))
			*V = 1;
		else if (!strncmp(Token, "NO", strlen(Token)))
			*V = 0;
		else
			Token = 0;
	}
	return Token;
}

static size_t max(size_t a, size_t b)
{
	return a > b ? a : b;
}