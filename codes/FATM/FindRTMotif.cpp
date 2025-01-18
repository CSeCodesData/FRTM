#include "FindRTMotif.h"

#pragma region FRTM

void FindRTMotif::FRTM(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel) {
#pragma region initialize
	int choiceStartT = graph->getStartT(), choiceEndT = graph->getEndT(),
		currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
	i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	int realSize = currNTimestamp - k + 1;
	vec(int)* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR
	
	vec(CComponentsII*) setCC;//cc list
	veciter(CComponentsII*) listEnd;//iterator for setCC
	DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

	i2iHMap subCCMap(nNode), root2Id(nNode);//used for recomputed ccs
	DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);//disjoint set used for recomputed ccs
	int* subCCId = DBG_NEW int[nEdge];

	vec(int) maxCC;//record the position of ccs which need to be checked
	
	
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map the position of cc in setCC to previous node in checkedCC
	int lastTime = choiceEndT - k + 1;//last start time 
	int Te;//the end time of interval
	int edgeEndT, vertexNum;
	int rightEndpoint;//max position of R set
	bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue; 
	//i2tupHMap expCheck;
	TGraph::posUsedForEdgeFilter = 0;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		BEGIN_TIMER(b)
		selectedNum = 0;
		rightEndpoint = 0;
		Te = Ts + k - 1;
		graph->edgeFilterFRTM(Ts, Te,
			edgeSetsR, selectedNum, rightEndpoint,
			k, fixLabel, isEdgeTypeFixed);
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;
		if (Ts != choiceStartT) {
			disjointSet->reset();
		}
		vertex2Pos.clear();
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		vertexNum = 0;
		edgeEndT = choiceEndT;

		//Test S(m,T)
		smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		//R2L
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT > rightEndpoint) continue;

			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();

			BEGIN_TIMER(c)
#pragma region maxCheck
				if (iterStart != iterEnd) {//new edges fetched
					/*maintain connectivity*/
					graph->maintainConnectivity(iterStart, iterEnd,
						vertex2Pos, disjointSet, vertexNum,
						combineCCPos);
					/*combine ccs*/
					graph->combineComponentsFRTM(setCC,
						vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
						/*tempComponentsSize, realMotifNum,*/
						combineCCPos);
				}
			/*update ccs and generate maximal RTMs*/
			graph->updateNewEdgeInfoFRTM(iterStart, iterEnd,
				setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				maxCC, Ts, edgeEndT, k, &TGraph::maintainTempCCForUF, &TGraph::updateCCFRTM, ComponentsD5Type::CCFRTM);//O(¦¤Es)
			Test::gm += END_TIMER(c);

#pragma endregion

			BEGIN_TIMER(d)
				graph->expCheckFRTM(maxCC, setCC,
					resultTF, k, Ts,
					edgeEndT, expandMask,/* expCheck, */motifNumber, &TGraph::checkExpandableFRTM);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();
		}

		//testing
		Test::updateMemoryUse();

		//release
		listEnd = setCC.end();
		for (veciter(CComponentsII*) listIter = setCC.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != nullptr)
				delete *listIter;
		}
		setCC.clear();
	}
	delete disjointSet;
	delete[] expandMask;
	delete[] subCCId;
	delete[] edgeSetsR;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}

void FindRTMotif::FRTMDYN(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int choiceStartT, int choiceEndT) {
#pragma region initialize
	int currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(int)* edgeSetsR = DBG_NEW vec(int)[choiceEndT - choiceStartT - k + 2];// R edge sets
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR

	vec(CComponentsII*) setCC;//cc list
	veciter(CComponentsII*) listEnd;//iterator for setCC
	DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode + 1);//disjoint set
	vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

	i2iHMap subCCMap, root2Id;//used for recomputed ccs
	DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);//disjoint set used for recomputed ccs
	int* subCCId = DBG_NEW int[nEdge];//recompute ccs

	vec(int) maxCC;//record the position of ccs which need to be checked

	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked; // map the position of cc in setCC to previous node in checkedCC
	int vertexNum;//the number of vertexs
	int lastTime = choiceEndT - k + 1;//last start time 
	int edgeEndT;//record the currently considering end time of R edge sets
	int Te;//the end time of interval
	int rightEndpoint;//max position of R set

	bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue; 
	//i2tupHMap expCheck;
	//BEGIN_TIMER(h)
	TGraph::posUsedForEdgeFilter = choiceStartT * nEdge;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		BEGIN_TIMER(b)
			selectedNum = 0;
		rightEndpoint = 0;
		Te = Ts + k - 1;
		graph->edgeFilterFRTMMidR(Ts, Te,
			edgeSetsR, selectedNum, rightEndpoint,
			k, fixLabel, isEdgeTypeFixed);
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;

		if (Ts != choiceStartT) {
			disjointSet->reset();
		}
		vertex2Pos.clear();
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		vertexNum = 0;
		edgeEndT = choiceEndT;
		//Test S(m,T)
		smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		//R2L
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT > rightEndpoint) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();

			BEGIN_TIMER(c)
#pragma region maxCheck
				if (iterStart != iterEnd) {
					/*maintain connectivity*/
					graph->maintainConnectivity(iterStart, iterEnd,
						vertex2Pos, disjointSet, vertexNum,
						combineCCPos);
					/*combine ccs*/
					graph->combineComponentsFRTM(setCC,
						vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
						/*tempComponentsSize, realMotifNum,*/
						combineCCPos);
				}

			/*update ccs and generate maximal RTMs*/
			graph->updateNewEdgeInfoFRTM(iterStart, iterEnd,
				setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				maxCC, Ts, edgeEndT, k, &TGraph::maintainTempCCForUF, &TGraph::updateCCFRTM, ComponentsD5Type::CCFRTM);//O(¦¤Es)
			Test::gm += END_TIMER(c);

#pragma endregion

			BEGIN_TIMER(d)
				graph->expCheckFRTM(maxCC, setCC,
					resultTF, k, Ts,
					edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableFRTMMidR);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();
		}
		//testing
		Test::updateMemoryUse();

		//release
		listEnd = setCC.end();
		for (veciter(CComponentsII*) listIter = setCC.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != nullptr) {
				CComponentsFRTM* listIt = (CComponentsFRTM*)*listIter;
				delete *listIter;
			}
		}
		setCC.clear();

	}
	graph->updateMidResult();
	delete disjointSet;
	delete tempDisjointSet;
	delete[] expandMask;
	delete[] subCCId;
	delete[] edgeSetsR;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}

void FindRTMotif::DFRTM(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int graphEndT) {
	int startT = graph->getStartT();
	int endT = graph->getEndT();
	int lastRow = graphEndT - k + 1;
#pragma region copy from original result
	for (int i = lastRow; i >= startT; i--) {
		int first = i + FindTMotif::k - 1;
		int tempPos = resultPos(i, graphEndT, startT, graphEndT, FindTMotif::k);
		if (i != startT) {
			int nowPos = resultPos(i, graphEndT, startT, endT, FindTMotif::k);
			for (int j = graphEndT; j >= first; j--, tempPos--, nowPos--) {
				int resultSize = (int)resultTF[tempPos].size();
				for (int s = 0; s < resultSize; s++) {
					resultTF[nowPos].emplace_back(resultTF[tempPos][s]);
				}
				resultTF[tempPos].clear();
				motifNumber += resultSize;
			}
		}
		else {
			for (int j = graphEndT; j >= first; j--, tempPos--) {
				int resultSize = (int)resultTF[tempPos].size();
				motifNumber += resultSize;
			}
		}

	}
#pragma endregion
	if (graphEndT == endT) return;
#pragma region update result
	//row number<=T-k+1
	graph->findRTMotifsDynamic(k, resultTF, graphEndT, motifNumber, fixLabel, isEdgeTypeFixed);
	//row number>T-k+1 
	FRTMDYN(graph, resultTF, fixLabel, lastRow + 1, endT);
#pragma endregion
}

#pragma endregion 


#pragma region FRTMOpt1

void FindRTMotif::FRTMOpt1(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel) {
#pragma region initialize
	int choiceStartT = graph->getStartT(), choiceEndT = graph->getEndT(),
		currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
	i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	int realSize = currNTimestamp - k + 1;
	vec(int)* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR
	vec(int)* edgeSetsRAdd = DBG_NEW vec(int)[realSize];// R+ edge sets for each interval
	unordered_set<int> edgesAdd(nEdge << 1);// R+ edge sets

	vec(CComponentsII*) setCC;//cc list
	veciter(CComponentsII*) listEnd;//iterator for setCC
	DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

	i2iHMap subCCMap(nNode), root2Id(nNode);//used for recomputed ccs
	DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);//disjoint set used for recomputed ccs
	int* subCCId = DBG_NEW int[nEdge];

	vec(int) maxCC;//record the position of ccs which need to be checked
	
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); //  map the position of cc in setCC to previous node in checkedCC
	int lastTime = choiceEndT - k + 1;//last start time 
	int Te;//the end time of interval
	int edgeEndT, vertexNum = 0;
	int rightEndpoint;//max position of R set
	bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)
	
	i2iHMap haveNewEdgeCC;//have edge in R+
	int* remainEdges = DBG_NEW int[nEdge];//record remained edges after removing edges to be removed
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue; 
	//i2tupHMap expCheck;
	TGraph::posUsedForEdgeFilter = 0;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		BEGIN_TIMER(b)
		selectedNum = 0;
		rightEndpoint = 0;
		Te = Ts + k - 1;
		graph->edgeFilterOpt1(Ts, Te, 
			edgeSetsRAdd, edgeSetsR, selectedNum, rightEndpoint, k, fixLabel, isEdgeTypeFixed);
		Test::compr += END_TIMER(b); 
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;

		if (Ts != choiceStartT) {
			disjointSet->reset();
		}
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		edgeEndT = choiceEndT;

		//Test S(m,T)
		smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		//smk = edgeSetsR[0].size();
		//maxSmk = max(maxSmk, smk);
		//allSmk += smk;
		
		//R2L
		edgesAdd.clear(); 
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT > rightEndpoint) continue;
			//maintain R+ sets
			for (auto add : edgeSetsRAdd[tempT]) {
				edgesAdd.emplace(add);
			}
			edgeSetsRAdd[tempT].clear();
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			BEGIN_TIMER(c)
#pragma region maxCheck
				if (iterStart != iterEnd) {
					/*maintain connectivity*/
					graph->maintainConnectivity(iterStart, iterEnd,
						vertex2Pos, disjointSet, vertexNum,
						combineCCPos);
					/*combine ccs*/
					graph->combineComponentsOpt1(setCC,
						vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
						combineCCPos);
				}
			/*update ccs and generate maximal RTMs*/
			graph->updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd, 
				setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				 maxCC, Ts, edgeEndT, k, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
			Test::gm += END_TIMER(c); 
#pragma endregion
			BEGIN_TIMER(d)
				graph->expCheckFRTM(maxCC, setCC,
					resultTF, k, Ts,
					edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();
		}
		//testing
		Test::updateMemoryUse();
		//release
		listEnd = setCC.end();
		for (veciter(CComponentsII*) listIter = setCC.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != nullptr)
				delete *listIter;
		}
		setCC.clear();
	}
	delete disjointSet;
	delete[] subCCId;
	delete tempDisjointSet;
	delete[] expandMask;
	delete[] edgeSetsR;
	delete[] edgeSetsRAdd;
	delete[] remainEdges;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}

void FindRTMotif::FRTMOpt1DYN(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int choiceStartT, int choiceEndT) {
#pragma region initialize
	int currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
	i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	int realSize = choiceEndT - choiceStartT - k + 2;
	
	vec(int)* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR
	vec(int)*edgeSetsRAdd = DBG_NEW vec(int)[realSize];// R+ edge sets for each interval
	unordered_set<int> edgesAdd(nEdge << 1);// R+ edge sets

	veciter(CComponentsII*) listEnd;//iterator for setCC
	vec(CComponentsII*) setCC;//cc list
	DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

	i2iHMap subCCMap(nNode), root2Id(nNode);//used for recomputed ccs
	DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);//disjoint set used for recomputed ccs
	int* subCCId = DBG_NEW int[nEdge];

	vec(int) maxCC;//record the position of ccs which need to be checked
	
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map ccpos in setCC to previous node in checkedCC
	int lastTime = choiceEndT - k + 1;//last start time 
	int Te;//the end time of interval
	int edgeEndT, vertexNum = 0;
	int rightEndpoint;//max position of R set
	bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)

	i2iHMap haveNewEdgeCC;//have edge in R+
	int* remainEdges = DBG_NEW int[nEdge];//record remained edges after removing edges to be removed
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue;
	//i2tupHMap expCheck;
	TGraph::posUsedForEdgeFilter = choiceStartT * nEdge;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		BEGIN_TIMER(b)
			selectedNum = 0;
		rightEndpoint = 0;
		Te = Ts + k - 1;
		graph->edgeFilterOpt1MidR(Ts, Te,
			edgeSetsRAdd, edgeSetsR, selectedNum, rightEndpoint,
			k, fixLabel, isEdgeTypeFixed);
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;

		if (Ts != choiceStartT) {
			disjointSet->reset();
		}
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		edgeEndT = choiceEndT;

		//Test S(m,T)
		smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		
		//R2L
		edgesAdd.clear();
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT > rightEndpoint) continue;

			//maintain R+ sets 
			for (auto add : edgeSetsRAdd[tempT]) {
				edgesAdd.emplace(add);
			}
			edgeSetsRAdd[tempT].clear();

			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();

			BEGIN_TIMER(c)
#pragma region maxCheck
				if (iterStart != iterEnd) {//new edges fetched
					/*maintain connectivity*/
					graph->maintainConnectivity(iterStart, iterEnd,
						vertex2Pos, disjointSet, vertexNum,
						combineCCPos);
					/*combine ccs*/
					graph->combineComponentsOpt1(setCC,
						vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
						combineCCPos);
				}

			/*update ccs and generate maximal RTMs*/
			graph->updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd, 
				setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				maxCC, Ts, edgeEndT, k, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
			Test::gm += END_TIMER(c);

#pragma endregion
			BEGIN_TIMER(d)
				graph->expCheckFRTM(maxCC, setCC,
					resultTF, k, Ts,
					edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1MidR);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();
		}

		//testing
		Test::updateMemoryUse();

		//release
		listEnd = setCC.end();
		for (veciter(CComponentsII*) listIter = setCC.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != nullptr)
				delete *listIter;
		}
		setCC.clear();
	}
	graph->updateMidResult();
	delete disjointSet;
	delete[] subCCId;
	delete tempDisjointSet;
	delete[] expandMask;
	delete[] remainEdges;
	delete[] edgeSetsR;
	delete[] edgeSetsRAdd;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}

void FindRTMotif::DFRTMOpt1(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int graphEndT) {
	int startT = graph->getStartT();
	int endT = graph->getEndT();
	int lastRow = graphEndT - k + 1;
#pragma region copy from original result
	for (int i = lastRow; i >= startT; i--) {
		int first = i + FindTMotif::k - 1;
		int tempPos = resultPos(i, graphEndT, startT, graphEndT, FindTMotif::k);
		if (i != startT) {
			int nowPos = resultPos(i, graphEndT, startT, endT, FindTMotif::k);
			for (int j = graphEndT; j >= first; j--, tempPos--, nowPos--) {
				int resultSize = (int)resultTF[tempPos].size();
				for (int s = 0; s < resultSize; s++) {
					resultTF[nowPos].emplace_back(resultTF[tempPos][s]);
				}
				resultTF[tempPos].clear();
				motifNumber += resultSize;
			}
		}
		else {
			for (int j = graphEndT; j >= first; j--, tempPos--) {
				int resultSize = (int)resultTF[tempPos].size();
				motifNumber += resultSize;
			}
		}
	}
#pragma endregion
	if (graphEndT == endT) return;
#pragma region update result
	//row number<=T-k+1
	graph->findRTMotifsDynamicOpt1(k, resultTF, graphEndT, motifNumber, fixLabel, isEdgeTypeFixed);
	//row number>T-k+1 
	FRTMOpt1DYN(graph, resultTF, fixLabel, lastRow + 1, endT);
#pragma endregion
}
#pragma endregion

#pragma region FRTMPlus
void FindRTMotif::FRTMPlus(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel) {
	//for short intervals
	if (Setting::delta == 0 || Setting::c == 0) {
		FRTMOpt1(graph, resultTF, fixLabel);
		return;
	}
	int newk = (int)(ceil(1 / Setting::delta) + 0.5);
	int maxLengthForShortIntv = newk - k;
	if (maxLengthForShortIntv >= 1) {

#pragma region initialize
		int choiceStartT = graph->getStartT(), choiceEndT = graph->getEndT(),
			currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
		i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
		/*map the root in disjoint set to the position of its corresponding
		connected component in the component list*/
		i2iHMap root2Comp(nNode);

		int realSize = currNTimestamp - k + 1;
		vec(int) * edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
		veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
		int selectedNum;//the number of edges in edgeSetsR
		vec(int) * edgeSetsRAdd = DBG_NEW vec(int)[realSize];// R+ edge sets for each interval
		unordered_set<int> edgesAdd(nEdge << 1);// R+ edge sets

		vec(CComponentsII*) setCC;//temporary cc list
		veciter(CComponentsII*) listEnd;//iterator for setCC
		DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
		vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/
		
		DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);//disjoint set used for recomputed ccs
		i2iHMap subCCMap(nNode), root2Id(nNode);//used for recomputed ccs
		int* subCCId = DBG_NEW int[nEdge];

		vec(int) maxCC;//record the position of ccs which need to be checked
		
		long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
		long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
		long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
		LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
		unordered_map<int, LinkedNode<int>*> hasChecked(nEdge);  //map the position of cc in setCC to previous node in checkedCC
		int lastTime = choiceEndT - k + 1;//last start time 
		int Te, TeShortInterval;//the end time of interval
		int edgeEndT, vertexNum = 0;
		int rightEndpoint;//max position of R set
		bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)

		i2iHMap haveNewEdgeCC;//have edge in R+
		int* remainEdges = DBG_NEW int[nEdge];//record remained edges after removing edges to be removed
#pragma endregion
		vector<pair<int, int>> tempRecordFromQueue;

		vec(int) * edgeSetsRShortIntv = DBG_NEW vec(int)[maxLengthForShortIntv];// R edge sets for short intervals
		veciter(int) iterEndShortIntv, iterStartShortIntv;//iterator for edgeSetsRShortIntv
		veciter(CComponentsShortIntv*) listEndShortIntv;//iterator for setCCShortIntv
		i2iHMap root2CompShortIntv(nNode);/*root2Comp for short intervals*/
		vec(CComponentsShortIntv*) setCCShortIntv;//cc list for short intervals
		DynamicConnectivity* disjointSetShortIntv = DBG_NEW DisjointSet(nNode);//disjoint set for short intervals
		vec(int) maxCCShortIntv;//record the position of ccs which need to be checked for short intervals
		vec(int) combineCCPosShortIntv;//record the position of ccs which need to combine with other ccs for short intervals
		int selectedNumShortIntv;
		int lastTimeShortIntv = choiceEndT - newk + 1;
		bool* hasE = DBG_NEW bool[nEdge];//whether the edge is in R[m,m+newk-1] ... R[m,T]
		int* newE = DBG_NEW int[nEdge];//edges not in R[m,m+newk-1] ... R[m,T]
		int newENum;//the number of edges in newE
		int* addE = DBG_NEW int[nEdge];//edges to add in set setCCShortIntv
		int addENum;//the number of edges in addE
		bool* saveR = DBG_NEW bool[nNode];//record ccs computed in large intervals to be saved
		//i2tupHMap expCheck;
		TGraph::posUsedForEdgeFilter = choiceStartT * nEdge;
		TGraph::posUsedForEdgeFilterShortIntv = (choiceStartT + k - 1)*nEdge;
		//Test::counter9 = -1;
		for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge, TGraph::posUsedForEdgeFilterShortIntv += nEdge) {//O(T)
			BEGIN_TIMER(b)
			selectedNum = 0;
			selectedNumShortIntv = 0;
			rightEndpoint = -1;
			TeShortInterval = Ts + k - 1;
			Te = Ts + newk - 1;
			if (Ts <= lastTimeShortIntv + 1) {
				CLEARALL(hasE, 0, nEdge, bool);
				CLEARALL(saveR, 0, nNode, bool);
				newENum = 0;
			}
			if (Ts <= lastTimeShortIntv) {//have large intervals and short intervals for [m,m+k-1]...[m,T]
				if (Ts != choiceStartT) {
					root2Comp.clear();
					//expCheck.clear();
					disjointSet->reset(vertexNum);
					hasChecked.clear();
					checkedCC->release();
				}
				graph->edgeFilterPlus(Ts, Te, TeShortInterval, maxLengthForShortIntv - 1, hasE, newE, newENum, 
					edgeSetsRAdd,
					edgeSetsR, selectedNum, edgeSetsRShortIntv, selectedNumShortIntv, rightEndpoint, newk, fixLabel, isEdgeTypeFixed);
			}
			else {//only short intervals for [m,m+k-1]...[m,T]
				graph->edgeFilterShortIntv(Ts, TeShortInterval, maxLengthForShortIntv - 1, edgeSetsRShortIntv, selectedNumShortIntv, fixLabel, isEdgeTypeFixed);
			}
			Test::compr += END_TIMER(b);
			maxSelect = max(maxSelect, selectedNum + selectedNumShortIntv);
			selectedSum += selectedNum + selectedNumShortIntv;

			if (Ts != choiceStartT) {
				disjointSetShortIntv->reset(vertexNum);
				root2CompShortIntv.clear();
			}
			//if (Ts <= 3000) {
			//	cout << Test::compr <<" "<< Test::gne << "#"<< endl;
			//}
			//if (Ts == 2001) return;
			edgeEndT = choiceEndT;
			//Test S(m,T)
			if (Ts <= lastTimeShortIntv) {
				smt = edgeSetsR[choiceEndT - Te].size();
			}
			else {
				smt = edgeSetsRShortIntv[choiceEndT - TeShortInterval].size();
			}
			maxSmt = max(maxSmt, smt);
			allSmt += smt;
			//R2L
			edgesAdd.clear();
			int tempT = choiceEndT - Te;
			for (int tempTAll = choiceEndT - TeShortInterval; tempTAll >= 0; tempTAll--, tempT--, edgeEndT--) {
				if (tempTAll >= maxLengthForShortIntv) {//large intervals
					if (tempT > rightEndpoint) continue;

					for (auto add : edgeSetsRAdd[tempT]) {
						edgesAdd.emplace(add);
					}
					edgeSetsRAdd[tempT].clear();
					/*if (Ts < 21001) {
						edgeSetsR[tempT].clear(); continue;
					}*/
					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();
					BEGIN_TIMER(c)
#pragma region maxCheck
						if (iterStart != iterEnd) {
							graph->maintainConnectivity(iterStart, iterEnd,
								vertex2Pos, disjointSet, vertexNum,
								combineCCPos);
							graph->combineComponentsOpt1(setCC,
								vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
								combineCCPos);
						}

					graph->updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,
						setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
						maxCC, Ts, edgeEndT, newk, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						graph->expCheckFRTM(maxCC, setCC,
							resultTF, k, Ts, 
							edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();
				}
				else {//short intervals
					/*if (Ts < 21001) {
						edgeSetsRShortIntv[tempTAll].clear(); continue;
					}*/
					if (Ts <= lastTimeShortIntv && tempTAll == maxLengthForShortIntv - 1) {//record root of disjointSet to be saved
						graph->recordSavedRoot(setCC, hasChecked, newE, newENum, saveR, disjointSet, vertex2Pos);
					}
					iterStartShortIntv = edgeSetsRShortIntv[tempTAll].begin();
					iterEndShortIntv = edgeSetsRShortIntv[tempTAll].end();

					if (iterStartShortIntv == iterEndShortIntv) continue;
					BEGIN_TIMER(c)
#pragma region maxCheck
						addENum = 0;
					graph->maintainConnectivityPlus(iterStartShortIntv, iterEndShortIntv, addE, addENum,
						vertex2Pos, disjointSetShortIntv, vertexNum, hasE, saveR, disjointSet, combineCCPosShortIntv);
					graph->combineComponentsPlus(setCCShortIntv,
						vertex2Pos, disjointSetShortIntv, root2CompShortIntv, combineCCPosShortIntv);
					graph->updateNewEdgeInfoPlus(iterStartShortIntv, addE, addENum,
						setCCShortIntv, vertex2Pos, disjointSetShortIntv,
						root2CompShortIntv, hasE, maxCCShortIntv, Ts, edgeEndT);
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						graph->expCheckFRTMPlus(maxCCShortIntv, setCCShortIntv,
							resultTF, k, Ts, edgeEndT, expandMask, /*expCheck,*/ motifNumber, disjointSet, hasChecked, root2Comp, vertex2Pos, setCC);
					Test::gne += END_TIMER(d);
					edgeSetsRShortIntv[tempTAll].clear();
				}
			}

			//testing
			Test::updateMemoryUse();
			//release
			listEnd = setCC.end();
			for (veciter(CComponentsII*) listIter = setCC.begin();
				listIter != listEnd; ++listIter) {
				if (*listIter != nullptr)
					delete *listIter;
			}
			setCC.clear();
			listEndShortIntv = setCCShortIntv.end();
			for (veciter(CComponentsShortIntv*) listIterShortIntv = setCCShortIntv.begin();
				listIterShortIntv != listEndShortIntv; ++listIterShortIntv) {
				if (*listIterShortIntv != nullptr)
					delete *listIterShortIntv;
			}
			setCCShortIntv.clear();
		}
		delete[] hasE;
		delete[] newE;
		delete[] addE;
		delete[] saveR;
		delete[] edgeSetsRShortIntv;
		delete disjointSetShortIntv;
		delete disjointSet;
		delete[] subCCId;
		delete tempDisjointSet;
		delete[] expandMask;
		delete[] edgeSetsR;
		delete[] edgeSetsRAdd;
		delete[] remainEdges;
		delete checkedCC;
		cout << SELECT_EDGE << maxSelect << endl;
		cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
		cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
	}
	else {
		FRTMOpt1(graph, resultTF, fixLabel);
	}
}


void FindRTMotif::FRTMPlusDYN(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int choiceStartT, int choiceEndT) {
	//for short intervals
	if (Setting::delta == 0 || Setting::c == 0) {
		FRTMOpt1DYN(graph, resultTF, fixLabel, choiceStartT, choiceEndT);
		return;
	}
	int newk = (int)(ceil(1 / Setting::delta) + 0.5);
	int maxLengthForShortIntv = newk - k;
	if (maxLengthForShortIntv >= 1) {
#pragma region initialize
		int currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
		i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
		/*map the root in disjoint set to the position of its corresponding
		connected component in the component list*/
		i2iHMap root2Comp(nNode);

		int realSize = currNTimestamp - k + 1;
		vec(int)* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
		veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
		int selectedNum;//the number of edges in edgeSetsR
		vec(int) * edgeSetsRAdd = DBG_NEW vec(int)[realSize];// R+ edge sets for each interval
		unordered_set<int> edgesAdd(nEdge << 1);// R+ edge sets
		
		vec(CComponentsII*) setCC;//cc list
		veciter(CComponentsII*) listEnd;//iterator for setCC
		DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
		vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

		DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);//disjoint set used for recomputed ccs
		i2iHMap subCCMap(nNode), root2Id(nNode);//used for recomputed ccs
		int* subCCId = DBG_NEW int[nEdge];

		vec(int) maxCC;//record the position of ccs which need to be checked

		long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
		long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
		long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
		LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
		unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); //map the position of cc in setCC to previous node in checkedCC
		int lastTime = choiceEndT - k + 1;//last start time 
		int Te, TeShortInterval;//the end time of interval
		int edgeEndT, vertexNum = 0;
		int rightEndpoint;//max position of R set
		bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)

		i2iHMap haveNewEdgeCC;//have edge in R+
		int* remainEdges = DBG_NEW int[nEdge];//record remained edges after removing edges to be removed
#pragma endregion
		vector<pair<int, int>> tempRecordFromQueue;

		vec(int) * edgeSetsRShortIntv = DBG_NEW vec(int)[maxLengthForShortIntv];// R edge sets for short intervals
		veciter(int) iterEndShortIntv, iterStartShortIntv;//iterator for edgeSetsRShortIntv
		veciter(CComponentsShortIntv*) listEndShortIntv;//iterator for setCCShortIntv
		i2iHMap root2CompShortIntv(nNode);/*root2Comp for short intervals*/
		vec(CComponentsShortIntv*) setCCShortIntv;//cc list for short intervals
		DynamicConnectivity* disjointSetShortIntv = DBG_NEW DisjointSet(nNode);//disjoint set for short intervals
		vec(int) maxCCShortIntv;//record the position of ccs which need to be checked for short intervals
		vec(int) combineCCPosShortIntv;//record the position of ccs which need to combine with other ccs for short intervals
		int selectedNumShortIntv;
		int lastTimeShortIntv = choiceEndT - newk + 1;
		bool* hasE = DBG_NEW bool[nEdge];//whether the edge is in R[m,m+newk-1] ... R[m,T]
		int* newE = DBG_NEW int[nEdge];//edges not in R[m,m+newk-1] ... R[m,T]
		int newENum;//the number of edges in newE
		int* addE = DBG_NEW int[nEdge];//edges to add in set setCCShortIntv
		int addENum;//the number of edges in addE
		bool* saveR = DBG_NEW bool[nNode];//record ccs computed in large intervals to be saved
		//i2tupHMap expCheck;

		TGraph::posUsedForEdgeFilter = choiceStartT * nEdge;
		TGraph::posUsedForEdgeFilterShortIntv = (choiceStartT + k - 1)*nEdge;
		for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge, TGraph::posUsedForEdgeFilterShortIntv += nEdge) {//O(T)
			BEGIN_TIMER(b)
			selectedNum = 0;
			selectedNumShortIntv = 0;
			rightEndpoint = -1;
			TeShortInterval = Ts + k - 1;
			Te = Ts + newk - 1;
			if (Ts == choiceStartT || Ts <= lastTimeShortIntv + 1) {
				CLEARALL(hasE, 0, nEdge, bool);
				CLEARALL(saveR, 0, nNode, bool);
				newENum = 0;
			}
			
			if (Ts <= lastTimeShortIntv) {//have large intervals and short intervals for [m,m+k-1]...[m,T]
				if (Ts != choiceStartT) {
					root2Comp.clear();
					disjointSet->reset(vertexNum);
					hasChecked.clear();
					checkedCC->release();
				}
				graph->edgeFilterPlusMidR(Ts, Te, TeShortInterval, maxLengthForShortIntv - 1, hasE, newE, newENum,
					edgeSetsRAdd, edgeSetsR, selectedNum, edgeSetsRShortIntv, selectedNumShortIntv, rightEndpoint, newk, fixLabel, isEdgeTypeFixed);
			}
			else {//only short intervals for [m,m+k-1]...[m,T]
				graph->edgeFilterShortIntvMidR(Ts, TeShortInterval, maxLengthForShortIntv - 1, lastTimeShortIntv, edgeSetsRShortIntv, selectedNumShortIntv, choiceEndT, fixLabel, isEdgeTypeFixed);
			}
			//initalize for every row
			Test::compr += END_TIMER(b);
			maxSelect = max(maxSelect, selectedNum + selectedNumShortIntv);
			selectedSum += selectedNum + selectedNumShortIntv;

			if (Ts != choiceStartT) {
				disjointSetShortIntv->reset(vertexNum);
				root2CompShortIntv.clear();
			}
			edgeEndT = choiceEndT;

			//Test S(m,T)
			if (Ts <= lastTimeShortIntv) {
				smt = edgeSetsR[choiceEndT - Te].size();
			}
			else {
				smt = edgeSetsRShortIntv[choiceEndT - TeShortInterval].size();
			}
			maxSmt = max(maxSmt, smt);
			allSmt += smt;
			//R2L
			
			edgesAdd.clear();
			int tempT = choiceEndT - Te;//for having noise
			for (int tempTAll = choiceEndT - TeShortInterval; tempTAll >= 0; tempTAll--, tempT--, edgeEndT--) {
				if (tempTAll >= maxLengthForShortIntv) {//large intervals
					if (tempT > rightEndpoint) continue;

					for (auto add : edgeSetsRAdd[tempT]) {
						edgesAdd.emplace(add);
					}
					edgeSetsRAdd[tempT].clear();

					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();

					BEGIN_TIMER(c)
#pragma region generateMaxTM
						if (iterStart != iterEnd) {
							graph->maintainConnectivity(iterStart, iterEnd,
								vertex2Pos, disjointSet, vertexNum,
								combineCCPos);
							graph->combineComponentsOpt1(setCC,
								vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
								combineCCPos);
						}

					graph->updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,/**/
						setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
						maxCC, Ts, edgeEndT, newk, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						graph->expCheckFRTM(maxCC, setCC,
							resultTF, k, Ts,
							edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1MidR);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();

				}
				else {//short intervals
					
					if (Ts <= lastTimeShortIntv && tempTAll == maxLengthForShortIntv - 1) {//record root of disjointSet to be saved
						graph->recordSavedRoot(setCC, hasChecked, newE, newENum, saveR, disjointSet, vertex2Pos);
					}
					iterStartShortIntv = edgeSetsRShortIntv[tempTAll].begin();
					iterEndShortIntv = edgeSetsRShortIntv[tempTAll].end();
					if (iterStartShortIntv == iterEndShortIntv) continue;
					BEGIN_TIMER(c)
#pragma region generateMaxTM
					addENum = 0;
					graph->maintainConnectivityPlus(iterStartShortIntv, iterEndShortIntv, addE, addENum,
							vertex2Pos, disjointSetShortIntv, vertexNum, hasE, saveR, disjointSet, combineCCPosShortIntv);
					graph->combineComponentsPlus(setCCShortIntv,
						vertex2Pos, disjointSetShortIntv, root2CompShortIntv, combineCCPosShortIntv);
					graph->updateNewEdgeInfoPlus(iterStartShortIntv, addE, addENum,
						setCCShortIntv, vertex2Pos, disjointSetShortIntv,
						root2CompShortIntv, hasE, maxCCShortIntv, Ts, edgeEndT);
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						graph->expCheckFRTMPlusMidR(maxCCShortIntv, setCCShortIntv,
							resultTF, k, Ts, edgeEndT, expandMask, /*expCheck,*/ motifNumber, disjointSet, hasChecked, root2Comp, vertex2Pos, setCC);
					Test::gne += END_TIMER(d);
					edgeSetsRShortIntv[tempTAll].clear();
				}
			}
			//testing
			Test::updateMemoryUse();
			//release
			listEnd = setCC.end();
			for (veciter(CComponentsII*) listIter = setCC.begin();
				listIter != listEnd; ++listIter) {
				if (*listIter != nullptr)
					delete *listIter;
			}
			setCC.clear();
			listEndShortIntv = setCCShortIntv.end();
			for (veciter(CComponentsShortIntv*) listIterShortIntv = setCCShortIntv.begin();
				listIterShortIntv != listEndShortIntv; ++listIterShortIntv) {
				if (*listIterShortIntv != nullptr)
					delete *listIterShortIntv;
			}
			setCCShortIntv.clear();
		}
		graph->updateMidResult();
		delete[] hasE;
		delete[] newE;
		delete[] addE;
		delete[] saveR;
		delete[] edgeSetsRShortIntv;
		delete disjointSetShortIntv;

		delete disjointSet;
		delete[] subCCId;
		delete tempDisjointSet;
		delete[] expandMask;
		delete[] edgeSetsR;
		delete[] edgeSetsRAdd;
		delete[] remainEdges;
		delete checkedCC;
		cout << SELECT_EDGE << maxSelect << endl;
		cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
		cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
	}
	else {
		FRTMOpt1DYN(graph, resultTF, fixLabel, choiceStartT, choiceEndT);
	}
}

void FindRTMotif::DFRTMPLus(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int graphEndT) {
	int startT = graph->getStartT();
	int endT = graph->getEndT();
	int lastRow = graphEndT - k + 1;
#pragma region copy from original result
	for (int i = lastRow; i >= startT; i--) {
		int first = i + FindTMotif::k - 1;
		int tempPos = resultPos(i, graphEndT, startT, graphEndT, FindTMotif::k);
		if (i != startT) {
			int nowPos = resultPos(i, graphEndT, startT, endT, FindTMotif::k);
			for (int j = graphEndT; j >= first; j--, tempPos--, nowPos--) {
				int resultSize = (int)resultTF[tempPos].size();
				for (int s = 0; s < resultSize; s++) {
					resultTF[nowPos].emplace_back(resultTF[tempPos][s]);
				}
				resultTF[tempPos].clear();
				motifNumber += resultSize;
			}
		}
		else {
			for (int j = graphEndT; j >= first; j--, tempPos--) {
				int resultSize = (int)resultTF[tempPos].size();
				motifNumber += resultSize;
			}
		}

	}
#pragma endregion
	if (graphEndT == endT) return;
#pragma region update result
	//row number<=T-k+1
	graph->findRTMotifsDynamicPlus(k, resultTF, graphEndT, motifNumber, fixLabel, isEdgeTypeFixed);
	//row number>T-k+1 
	FRTMPlusDYN(graph, resultTF, fixLabel, lastRow + 1, endT);
#pragma endregion
}

#pragma endregion

