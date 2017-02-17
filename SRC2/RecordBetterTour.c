#include "LKH.h"

/*
 * The RecordBetterTour function is called by FindTour each time
 * the LinKernighan function has returned a better tour.
 *
 * The function records the tour in the BetterTour array and in the
 * BestSuc field of each node. Furthermore, for each node the previous 
 * value of BestSuc is saved in the NextBestSuc field.
 *
 * Recording a better tour in the BetterTour array when the problem is 
 * asymmetric requires special treatment since the number of nodes has
 * been doubled.  
 */
/*
    每当LinKernighan()函数找到一个更好的解就会调用RecordBetter()函数
    RecordBetterTour()函数会把这个更好的解记录在BetterTour[]数组中。如果这个数组已经有值，就把原来的
    值存在NextBestSuc[]数组中，然后才更新。
 */

void RecordBetterTour()
{
    Node *N = FirstNode, *Stop = N;

    if (ProblemType != ATSP) {
        int i = 1;
        do
            BetterTour[i++] = N->Id;
        while ((N = N->Suc) != Stop);
    } else {
        if (Stop->Id > DimensionSaved)
            Stop = N = Stop->Suc;
        if (N->Suc->Id != DimensionSaved + N->Id) {
            int i = 1;
            do
                if (N->Id <= DimensionSaved)
                    BetterTour[i++] = N->Id;
            while ((N = N->Suc) != Stop);
        } else {
            int i = DimensionSaved;
            do
                if (N->Id <= DimensionSaved)
                    BetterTour[i--] = N->Id;
            while ((N = N->Suc) != Stop);
        }
    }
    BetterTour[0] = BetterTour[ProblemType != ATSP ?
                               Dimension : DimensionSaved];
    N = FirstNode;
    do {
        N->NextBestSuc = N->BestSuc;
        N->BestSuc = N->Suc;
    }
    while ((N = N->Suc) != FirstNode);
}
