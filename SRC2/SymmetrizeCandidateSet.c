#include "LKH.h"

/*
 * The SymmetrizeCandidateSet function complements the candidate set such
 * that every candidate edge is associated with both its two end nodes.
*/
/*
 SymmetrizeCandidateSet()函数用来扩充候选集，目的是让每一条候选边都与它两端的节点关联
 */
void SymmetrizeCandidateSet()
{
	Node *From, *To;
	Candidate *NFrom;

	From = FirstNode;
	do {
		for (NFrom = From->CandidateSet; NFrom && (To = NFrom->To); NFrom++)
			/*
			 AddCandidate()函数会把(To,From)这条边放到节点To的候选边集合中。
			 实参NFrom->Cost和NFrom->Alpha表示(To,From)这条边的cost和Alpha值
			*/
			AddCandidate(To, From, NFrom->Cost, NFrom->Alpha);
	}
	while ((From = From->Suc) != FirstNode);
	// ResetCandidates()函数会移除候选边集合中正在使用的边和Alpha等于正无穷的边，然后将集合中的边重新排序
	ResetCandidateSet();
}
