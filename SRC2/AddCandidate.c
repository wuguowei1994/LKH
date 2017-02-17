#include "LKH.h"

/*
 * The AddCandidate function adds a given edge (From, To) to the set
 * of candidate edges associated with the node From. The cost and
 * alpha-value of the edge are passed as parameters to the function.
 *
 * The function has no effect if the edge is already in the candidate
 * set.
 *
 * If the edge was added, the function returns 1; otherwise 0.
 *
 * The function is called from the functions CreateDelaunaySet and
 * OrderCandidateSet.
 */
/*
  AddCandidate()函数会把(from,to)这条边放到节点from的候选边集合中。行参Cost和Alpha表示
  (from,to)这条边的cost和Alpha值
  如果这条边已经在from节点的候选集中，该函数不会造成任何影响
  如果添加成功，函数返回1;否则返回0
 */

int AddCandidate(Node * From, Node * To, int Cost, int Alpha)
{
    int Count;
    Candidate *NFrom;
    // 这个if不会进入
    if (From->Subproblem != FirstNode->Subproblem)
        return 0;
    // 这个if语句不会进入
    if (From->CandidateSet == 0)
        assert(From->CandidateSet =
                   (Candidate *) calloc(3, sizeof(Candidate)));

    //这个if语句不会进入
    if (From == To || To->Subproblem != FirstNode->Subproblem ||
            !IsPossibleCandidate(From, To))
        return 0;
    // 函数从下面才真正开始
    Count = 0;
    for (NFrom = From->CandidateSet; NFrom->To && NFrom->To != To; NFrom++)
        Count++;
    if (NFrom->To) {
        if (NFrom->Alpha == INT_MAX)
            NFrom->Alpha = Alpha;
        return 0;
    }
    NFrom->Cost = Cost;
    NFrom->Alpha = Alpha;
    NFrom->To = To;
    assert(From->CandidateSet =
               (Candidate *) realloc(From->CandidateSet,
                                     (Count + 2) * sizeof(Candidate)));
    From->CandidateSet[Count + 1].To = 0;
    return 1;
}
