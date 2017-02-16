#include "LKH.h"


/*
  AddTourCandidates()函数会从用户指定的文件中添加候选边集，由于默认情况下不指定，所以AddTourCandidates()函数其实什么都没有做
  这个函数会被GenerateCandidateSet()函数和OrderCandidateSet()函数调用
 */

void AddTourCandidates()
{
    Node *Na, *Nb;
    int i, d, Subproblem = FirstNode->Subproblem;

    //添加fixed edges
    Na = FirstNode;
    do {
        if (Na->FixedTo1)
            AddCandidate(Na, Na->FixedTo1, D(Na, Na->FixedTo1), 0);
        if (Na->FixedTo2)
            AddCandidate(Na, Na->FixedTo2, D(Na, Na->FixedTo2), 0);
    }
    while ((Na = Na->Suc) != FirstNode);

    //添加MERGE_TOUR_FILE这个文件中的边
    //默认情况下,MergeTourFiles为零，所以这个for循环一次都不会进入
    for (i = 0; i < MergeTourFiles; i++) {
        Na = FirstNode;
        do {
            Nb = Na->MergeSuc[i];
            if (!Nb)
                break;
            if (Na->Subproblem == Subproblem &&
                    Nb->Subproblem == Subproblem) {
                d = D(Na, Nb);
                AddCandidate(Na, Nb, d, 1);
                AddCandidate(Nb, Na, d, 1);
            }
        }
        while ((Na = Nb) != FirstNode);
    }

    // 把INITIAL_TOUR_FILE这个文件中的边添加进来，这个do while循环内部的第二个if语句永远不会进入。其实什么都没有做
    Na = FirstNode;
    do {
        Nb = Na->InitialSuc;
        if (!Nb)
            break;
        if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
            d = D(Na, Nb);
            AddCandidate(Na, Nb, d, 1);
            AddCandidate(Nb, Na, d, 1);
        }
    }
    while ((Na = Nb) != FirstNode);

    /* Add INPUT_TOUR_FILE edges */
    Na = FirstNode;
    do {
        Nb = Na->InputSuc;
        if (!Nb)
            break;
        if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
            d = D(Na, Nb);
            AddCandidate(Na, Nb, d, 1);
            AddCandidate(Nb, Na, d, 1);
        }
    }
    while ((Na = Nb) != FirstNode);

    /* Add SUBPROBLEM_TOUR_FILE edges */
    Na = FirstNode;
    do {
        Nb = Na->SubproblemSuc;
        if (!Nb)
            break;
        if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
            d = D(Na, Nb);
            AddCandidate(Na, Nb, d, 1);
            AddCandidate(Nb, Na, d, 1);
        }
    } while ((Na = Nb) != FirstNode);
}
