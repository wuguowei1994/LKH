#include "LKH.h"

/*
 * Each time a trial has resulted in a shorter tour the candidate set is
 * adjusted (by AdjustCandidateSet). The ResetCandidates function resets
 * the candidate set. The original order is re-established (using, and 
 * edges with Alpha == INT_MAX are excluded.
 *
 * The function is called from FindTour and OrderCandidates.
 */
/*
 每当程序找到更短的路径时，就会调整候选边集合(通过AdjustCandidateSet()这个函数)。
 ResetCandidates()函数会移除候选边集合中正在使用的边和Alpha等于正无穷的边，然后将集合中的边重新排序
 这个函数会被FindTour()和OrderCandidates()函数调用
 */
void ResetCandidateSet()
{
    Node *From;
    Candidate *NFrom, *NN, Temp;

    From = FirstNode;
    /* Loop for all nodes */
    do {
        if (!From->CandidateSet)
            continue;
        /* Reorder the candidate array of From */
        for (NFrom = From->CandidateSet; NFrom->To; NFrom++) {
            Temp = *NFrom;
            for (NN = NFrom - 1;
                 NN >= From->CandidateSet &&
                 (Temp.Alpha < NN->Alpha ||
                  (Temp.Alpha == NN->Alpha && Temp.Cost < NN->Cost)); NN--)
                *(NN + 1) = *NN;
            *(NN + 1) = Temp;
        }
        NFrom--;
        /* Remove included edges */
        while (NFrom >= From->CandidateSet + 2 && NFrom->Alpha == INT_MAX)
            NFrom--;
        NFrom++;
        NFrom->To = 0;
        /* Remove impossible candidates */
        for (NFrom = From->CandidateSet; NFrom->To; NFrom++) {
            if (!IsPossibleCandidate(From, NFrom->To)) {
                for (NN = NFrom; NN->To; NN++)
                    *NN = *(NN + 1);
                NFrom--;
            }
        }
    }
    while ((From = From->Suc) != FirstNode);
}
