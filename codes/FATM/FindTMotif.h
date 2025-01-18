#pragma once
#include "stdafx.h"
#include "TGraph.h"
#include "TMotif.h"
#include "DynamicConnectivity.h"
#include "LinkedList.h"
#include "Util.h"
#include "DisjointSet.h"

/*implement the algorithm for exact label matches*/
class FindTMotif {
public:
	
	/*
		static algorithm for exact label matches
		@parameter:
		graph: input temporal graph
		result: output result
		choiceStartT£¬choiceEndT: start timestamp and end timestamp of the algorithm
	*/
	static void FRTMExact(TGraph*& graph,
		 vec(TMotifI*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT);//precise
	
	#pragma region function of printing

	static bool cmp(TMotifI*& a, TMotifI*& b) {
		sort(a->getMotifEdge()->begin(), a->getMotifEdge()->end());
		sort(b->getMotifEdge()->begin(), b->getMotifEdge()->end());
		return a->getMotifEdge()->begin()->id < b->getMotifEdge()->begin()->id;
	}
	static bool cmp2(TMotifII*& a, TMotifII*& b) {
		sort(a->getMotifEdge()->begin(), a->getMotifEdge()->end());
		sort(b->getMotifEdge()->begin(), b->getMotifEdge()->end());
		return *(a->getMotifEdge()->begin()) < *(b->getMotifEdge()->begin());
	}

	/*print the list of motifs*/
	//static void print(TGraph*& temporal_graph, vec(TMotifI*)*& res, int pos, bool outputTime) {
	//	vec(TMotifI*)& lis = res[pos];
	//	if (lis.size() == 0) {
	//		cout << "empty result\n" << endl;
	//		return;
	//	}
	//	//sort(lis.begin(), lis.end(), cmp);//testing
	//	veciter(TMotifI*) listIter = lis.begin(),
	//		listEnd = lis.end();
	//	int i = 1;
	//	TMotifI* motif;
	//	int motifStartT, motifEndT, intvLen;
	//	int size;
	//	if (listIter != listEnd) {
	//		motifStartT = (*listIter)->getStartT();
	//		motifEndT = (*listIter)->getEndT();
	//		intvLen = motifEndT - motifStartT + 1;
	//		Test::maxIntvLen = Test::maxIntvLen < intvLen ? intvLen : Test::maxIntvLen;
	//		Test::sumIntvLen += lis.size() * intvLen;
	//	}
	//	int counter = 0;
	//	for (; listIter != listEnd; ++listIter) {
	//		motif = *listIter;
	//		size = (int)motif->getSize();
	//		motifSum += size;
	//		motifMaxNum = size <= motifMaxNum ? motifMaxNum : size;
	//		
	//		//if (size < filterEdgeNumber) continue;
	//		if (output >= 1) {
	//			if (outputTime) {
	//				cout << "startT: " << motifStartT
	//					<< "\tendT: " << motifEndT << endl;
	//			}
	//			if (output >= 2) {
	//				temporal_graph->printMotif(res, motif, i, k/*, fixednode1, fixednode2*/);
	//				i++;
	//			}
	//			else if (output == 1) {
	//				cout << MOTIF_ID << i++ << endl;
	//				cout << EDGE_NUM << size << endl;
	//			}
	//		}
	//	}
	//	if (output != 0) cout << "\n";
	//}

	static void print(TGraph*& temporal_graph, vec(TMotifII*)*& res, int pos, bool outputTime) {
		vec(TMotifII*)& lis = res[pos];
		if (lis.size() == 0) {
			cout << "empty result\n" << endl;
			return;
		}
		sort(lis.begin(), lis.end(), cmp2);//testing

		veciter(TMotifII*) listIter = lis.begin(),
			listEnd = lis.end();
		int i = 1;
		TMotifII* motif;
		int motifStartT, motifEndT, intvLen;
		int size;
		if (listIter != listEnd) {
			motifStartT = (*listIter)->getStartT();
			motifEndT = (*listIter)->getEndT();
			intvLen = motifEndT - motifStartT + 1;
			Test::maxIntvLen = Test::maxIntvLen < intvLen ? intvLen : Test::maxIntvLen;
			Test::sumIntvLen += lis.size() * intvLen;
		}
		int counter = 0;
		for (; listIter != listEnd; ++listIter) {
			motif = *listIter;
			size = (int)motif->getSize();
			motifSum += size;
			motifMaxNum = size <= motifMaxNum ? motifMaxNum : size;
			
			//if (size < filterEdgeNumber) continue;

			if (output >= 1) {
				if (outputTime) {
					cout << "startT: " << motifStartT
						<< "\tendT: " << motifEndT << endl;
				}
				if (output >= 3) {
					temporal_graph->printMotif2(motif, i/*, fixednode1, fixednode2*/);
					i++;
				}
				else if (output == 2) {
					temporal_graph->printMotif(motif, i/*, fixednode1, fixednode2*/);
					i++;
				}
				else if (output == 1) {
					cout << MOTIF_ID << i++ << endl;
					cout << EDGE_NUM << size << endl;
				}
			}
		}
		if (output != 0) cout << "\n";
	}
	

	#pragma endregion

	static int k;//frequency condition

	static long long motifNumber;//the number of motifs

	static int output;//output level
						   
	static bool isEdgeTypeFixed;//whether fix the edge's label of motifs(default:false)

	static char outputSrc[FILE_NAME_LENGTH];//output file name

	static long long  motifMaxNum, motifSum;

	//static int filterIntvSize;//output motifs with interval length > filterIntvSize, default filterIntvSize = 0
	//static int filterEdgeNumber;//output motifs with edge number > filterIntvSize, default filterEdgeNumber = 0
	//static int* edgeBef;
	//static int fixednode1, fixednode2;//output motifs with a specific edge
};


