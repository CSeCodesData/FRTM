#include "FindTMotif.h"


void FindTMotif::FRTMExact(TGraph*& graph, vec(TMotifI*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT) {
	int nEdge = graph->getNEdge();

#pragma region initialize
	int maxIntervalLength = graph->getMaxIntervalLength();

	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(CComponents*) tempComponents;//temporary component list CC
	SAVEINFO_Vec* edgeSetsR = DBG_NEW SAVEINFO_Vec[maxIntervalLength];// R edge sets
	DynamicConnectivity* disjointSet;//disjoint set

	SAVEINFO_VIter iterEnd, iterStart;//iterator for edgeSetsR
	veciter(CComponents*) listEnd;//iterator for tempComponents
	int selectedNum;//the number of edges in edgeSetsR

	vec(int) saveMotifPos;//record the position of components which need to be saved
	/*record the position of components which need to combine with other components*/
	vec(int) combineMotifPos;

	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]

	int vertexNum;//the number of vertexs
	//int realMotifNum;//the number of generated connected component
	int lastTime = choiceEndT - k + 1;//last start time 
	int edgeEndT;//record the currently considering end time of R edge sets
	int Te;//the end time of interval
#pragma endregion
	TGraph::posUsedForEdgeFilter = (choiceStartT + k - 1)*nEdge;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		// select edges which keep fixed in [Ts,Ts+k-1]
		auto begin = clock();
		selectedNum = 0;
		Te = Ts + k - 1;
		graph->edgeFilter(Ts, Te,
			edgeSetsR, selectedNum,
			fixLabel, isEdgeTypeFixed);
		//initalize for every row
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;
		disjointSet = DBG_NEW DisjointSet((int)((selectedNum << 1) + 1));
		vertex2Pos.clear();
		root2Comp.clear();
		vertexNum = 0;
		edgeEndT = choiceEndT;
		//cout << clock() - begin << endl;
		Test::compr += clock() - begin;

		//Test S(m,T)
		if (maxIntervalLength <= choiceEndT - Te) smt = 0;
		else smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;

		//R2L
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT >= maxIntervalLength) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			if (iterStart == iterEnd) continue;
			//generateMaxTM and generatedExpTM
			graph->genMotifInOneIntv(iterStart, iterEnd,
				vertex2Pos, disjointSet, vertexNum,
				combineMotifPos, /*realMotifNum,*/ root2Comp, tempComponents,
				saveMotifPos, Ts, edgeEndT,
				result, k, motifNumber);
			edgeSetsR[tempT].clear();
		}

		//testing
		Test::updateMemoryUse();

		//release
		delete disjointSet;
		listEnd = tempComponents.end();
		for (veciter(CComponents*) listIter = tempComponents.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != nullptr)
				delete *listIter;
		}
		tempComponents.clear();
	}
	delete[] edgeSetsR;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}
