#include "LKH.h"

/*
 * The WritePenalties function writes node penalties (Pi-values)
 * to file PiFileName.
 *
 * The first line of the file contains the number of nodes.
 *
 * Each of the following lines is of the form
 *       <integer> <integer>
 * where the first integer is a node number, and the second integer
 * is the Pi-value associated with the node.
 *
 * The function is called from the CreateCandidateSet function.
 */
/*
  WritePenalties()函数会把Ascent()函数中计算出来的每个节点的Pi值写入到PiFileName文件中
  PiFileName第一行表示节点的数量，以后的每一行包含两个integer，第一个表示节点的索引，第二个表示
  这个节点的Pi值。
 */

void WritePenalties()
{
    Node *N;
    if (PiFileName == 0 || !(PiFile = fopen(PiFileName, "w")))
        return;
    if (TraceLevel >= 1)
        printff("Writing PI_FILE: \"%s\" ... ", PiFileName);
    fprintf(PiFile, "%d\n", Dimension);
    N = FirstNode;
    do
        fprintf(PiFile, "%d %d\n", N->Id, N->Pi);
    while ((N = N->Suc) != FirstNode);
    fprintf(PiFile, "-1\nEOF\n");
    fclose(PiFile);
    if (TraceLevel >= 1)
        printff("done\n", PiFileName);
}
