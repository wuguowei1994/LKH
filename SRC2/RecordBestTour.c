#include "LKH.h"

/*
 * The RecordBestTour function records the current best tour in the BestTour 
 * array. 
 *
 * The function is called by LKHmain each time a run has resulted in a
 * shorter tour. Thus, when the predetermined number of runs have been
 * completed, BestTour contains an array representation of the best tour
 * found.    
 */
/*
  RecordBestTour()函数会把当前的最优解记录到BestTour[]数组中
 */

void RecordBestTour()
{
    int i, Dim = ProblemType != ATSP ? Dimension : Dimension / 2;

    for (i = 0; i <= Dim; i++)
        BestTour[i] = BetterTour[i];
}
