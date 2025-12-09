#include "FindRTMotif.h"

void FindRTMotif::reorganizeResult(TGraph*& graph, vec(TMotifII*)*& resultTF, int graphEndT, int choiceEndT) {
	int startT = graph->getStartT();
	int lastRow = graphEndT - k + 1;
#pragma region copy from original result
	for (int i = lastRow; i >= startT; i--) {
		int first = i + FindTMotif::k - 1;
		int tempPos = resultPos(i, graphEndT, startT, graphEndT, FindTMotif::k);
		if (i != startT) {
			int nowPos = resultPos(i, graphEndT, startT, choiceEndT, FindTMotif::k);
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
}
#pragma region FRTM

void FindRTMotif::FRTM(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int choiceStartT, int choiceEndT, bool dynamicMode, ComponentsType cctype) {
#pragma region initialize
	int currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
	//i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	int realSize = currNTimestamp - k + 1;
	vec(int)* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR
	
	CComponentsFRTM* emptyCC = nullptr;
	vec(CComponentsII*) setCC;//cc list
	veciter(CComponentsII*) listEnd;//iterator for setCC
	DisjointSet* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

	i2iHMap root2Id(nNode);//used for recomputed ccs
	DisjointSetWithTime* tempDisjointSet = DBG_NEW DisjointSetWithTime(nNode);//disjoint set used for recomputed ccs
	int* subCCId = DBG_NEW int[nEdge];

	vec(int) maxCC;//record the position of ccs which need to be checked
	
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map the position of cc in setCC to previous node in checkedCC
	int lastTime = choiceEndT - k + 1;//last start time 
	int Te;//the end time of interval
	int edgeEndT/*, vertexNum*/;
	int rightEndpoint;//max position of R set
	bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue; 
	//i2tupHMap expCheck;
	TGraph::posUsedForEdgeFilter = choiceStartT * nEdge;

	auto checkExpandableFun = dynamicMode?
		&TGraph::checkExpandableFRTMMidR<CComponentsFRTM> : &TGraph::checkExpandableFRTM<CComponentsFRTM>;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		BEGIN_TIMER(b)
		selectedNum = 0;
		rightEndpoint = 0;
		Te = Ts + k - 1;
		graph->edgeFilterFRTM(Ts, Te,
			edgeSetsR, selectedNum, rightEndpoint,
			k, fixLabel, isEdgeTypeFixed, dynamicMode);
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;
		if (Ts != choiceStartT) {
			disjointSet->reset();
		}
		//vertex2Pos.clear();
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		//vertexNum = 0;
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
						/*vertex2Pos,*/ disjointSet, /*vertexNum,*/
						combineCCPos);
					/*combine ccs*/
					graph->combineComponentsFRTM(setCC,
						/*vertex2Pos,*/ disjointSet, root2Comp, checkedCC, hasChecked,
						combineCCPos);
				}
			/*update ccs and generate maximal RTMs*/
			graph->updateNewEdgeInfoFRTM(iterStart, iterEnd,
				setCC,/* vertex2Pos,*/ disjointSet, tempDisjointSet, root2Comp, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				maxCC, Ts, edgeEndT, k, cctype);//O(¦¤Es)
			Test::gm += END_TIMER(c);

#pragma endregion

			BEGIN_TIMER(d)
			graph->expCheckFRTM(maxCC, setCC, resultTF, k, Ts, edgeEndT, expandMask, motifNumber, 
				checkExpandableFun, emptyCC);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();

			//cout << Ts << " " << edgeEndT << " " << Test::counter << " " << Test::counter2 << " " << Test::counter3 << 
			// " " << Test::counter4 << " " << Test::counter5<< endl;
		}

		//testing
		Test::updateMemoryUse();

		//release
		graph->releaseCCs(setCC);
	}
	if (dynamicMode) { 
		graph->updateMidResult(); 
	}
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
	bool*& fixLabel, int graphEndT, ComponentsType cctype) {
	int endT = graph->getEndT();
	int lastRow = graphEndT - k + 1;
	
	reorganizeResult(graph, resultTF, graphEndT, endT);
	if (graphEndT == endT) return;

#pragma region update result
	//row number<=T-k+1
	graph->findRTMotifsDynamic(k, resultTF, graphEndT, motifNumber, fixLabel, isEdgeTypeFixed, cctype);
	//row number>T-k+1 
	FRTM(graph, resultTF, fixLabel, lastRow + 1, endT, true, cctype);
#pragma endregion
}

#pragma endregion 


#pragma region FRTMOpt1

void FindRTMotif::FRTMOpt1(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int choiceStartT, int choiceEndT, bool dynamicMode, ComponentsType cctype) {
#pragma region initialize
	int currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
	//i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	int realSize = currNTimestamp - k + 1;
	vec(int)* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR
	int* edgesAdd = DBG_NEW int[nEdge];// R+ edge sets
	CLEARALL(edgesAdd, 0, nEdge, int);

	CCTYPEOPT1* emptyCC = nullptr;
	vec(CComponentsII*) setCC;//cc list
	DisjointSet* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	
	vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/

	i2iHMap root2Id(nNode);//used for recomputed ccs
	
	DisjointSetWithTime* tempDisjointSet = nullptr;
	int* subCCId = nullptr;
	bool* expandMask = nullptr;
	int* tempHaveNewEdgeCC = nullptr;
	if (Setting::delta != 0 && Setting::c != 0) {
		tempDisjointSet = DBG_NEW DisjointSetWithTime(nNode);//disjoint set used for recomputed ccs
		subCCId = DBG_NEW int[nEdge];
		expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)
		tempHaveNewEdgeCC = DBG_NEW int[nNode];//have edge in R+
	}
	
	vec(int) maxCC;//record the position of ccs which need to be checked
	
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); //  map the position of cc in setCC to previous node in checkedCC
	int lastTime = choiceEndT - k + 1;//last start time 
	int Te;//the end time of interval
	int edgeEndT/*, vertexNum = 0*/;
	int rightEndpoint;//max position of R set
	
	int* haveNewEdgeCC = DBG_NEW int[nEdge];//have edge in R+
	CLEARALL(haveNewEdgeCC,-1,nEdge,int);
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue; 
	TGraph::posUsedForEdgeFilter = choiceStartT * nEdge;
	auto checkExpandableFun = dynamicMode ? 
		&TGraph::checkExpandableFRTMMidR<CCTYPEOPT1> : &TGraph::checkExpandableFRTM<CCTYPEOPT1>;
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge) {//O(T)
		//BEGIN_TIMER(b)
		selectedNum = 0;
		rightEndpoint = 0;
		Te = Ts + k - 1;
		graph->edgeFilterFRTM(Ts, Te, 
			edgeSetsR, selectedNum, rightEndpoint, k, fixLabel, isEdgeTypeFixed, dynamicMode, edgesAdd);
		//Test::compr += END_TIMER(b);
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
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT > rightEndpoint) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			//BEGIN_TIMER(c)
#pragma region maxCheck
				if (iterStart != iterEnd) {
					/*maintain connectivity*/
					graph->maintainConnectivity(iterStart, iterEnd,
						/*vertex2Pos,*/ disjointSet, /*vertexNum,*/
						combineCCPos);
					/*combine ccs*/
					graph->combineComponentsFRTMOPT1(setCC, /*vertex2Pos,*/ disjointSet, root2Comp, checkedCC, 
						hasChecked, combineCCPos, haveNewEdgeCC);
				}
			/*update ccs and generate maximal RTMs*/
			graph->updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd, 
				setCC, /*vertex2Pos,*/ disjointSet, tempDisjointSet, root2Comp, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				 maxCC, Ts, edgeEndT, k, haveNewEdgeCC, tempHaveNewEdgeCC, /*remainEdges,*/ cctype);//O(¦¤Es)
			//Test::gm += END_TIMER(c); 
#pragma endregion
			//BEGIN_TIMER(d)
			graph->expCheckFRTM(maxCC, setCC, resultTF, k, Ts, edgeEndT, expandMask, motifNumber, checkExpandableFun, emptyCC);
			//Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();

			//cout << Ts<<" "<<edgeEndT<<" "<<Test::counter << " " << Test::counter2 << " " << 
			// Test::counter3 << " " << Test::counter4 << " " << Test::counter5<< endl;
		}
		//testing
		Test::updateMemoryUse();
			//release
		graph->releaseCCs(setCC);
	}
	if (dynamicMode) {
		graph->updateMidResult();
	}
	delete disjointSet;
	if (Setting::delta != 0 && Setting::c != 0) {
		delete[] subCCId;
		delete tempDisjointSet;
		delete[] expandMask;
		delete[] tempHaveNewEdgeCC;
	}
	delete[] edgeSetsR;
	delete[] edgesAdd;
	delete[] haveNewEdgeCC;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}


void FindRTMotif::DFRTMOpt1(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int graphEndT, ComponentsType cctype) {
	int endT = graph->getEndT();
	int lastRow = graphEndT - k + 1;

	reorganizeResult(graph, resultTF, graphEndT, endT);
	if (graphEndT == endT) return;

#pragma region update result
	//row number<=T-k+1
	graph->findRTMotifsDynamicOpt1(k, resultTF, graphEndT, motifNumber, fixLabel, isEdgeTypeFixed, cctype);
	//row number>T-k+1 
	FRTMOpt1(graph, resultTF, fixLabel, lastRow + 1, endT, true, cctype);
#pragma endregion
}
#pragma endregion

#pragma region FRTMPlus
void FindRTMotif::FRTMPlus(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int choiceStartT, int choiceEndT, bool dynamicMode, ComponentsType cctype) {
	//for short intervals
	if (Setting::delta == 0 || Setting::c == 0) {
		FRTMOpt1(graph, resultTF, fixLabel, choiceStartT, choiceEndT, dynamicMode, cctype);
		return;
	}
	int newk = (int)(ceil(1 / Setting::delta) + 0.5);
	int maxLengthForShortIntv = newk - k;
	if (maxLengthForShortIntv >= 1) {

#pragma region initialize
		int currNTimestamp = graph->getCurrNTimestamp(), nNode = graph->getNNode(), nEdge = graph->getNEdge();
		//i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
		/*map the root in disjoint set to the position of its corresponding
		connected component in the component list*/
		i2iHMap root2Comp(nNode);

		int realSize = currNTimestamp - k + 1;
		vec(int) * edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
		veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
		int selectedNum;//the number of edges in edgeSetsR
		//unordered_set<int> edgesAdd((nEdge << 1));// R+ edge sets
		int* edgesAdd = DBG_NEW int[nEdge];// R+ edge sets
		CLEARALL(edgesAdd, 0, nEdge, int);

		CCTYPEOPT1* emptyCC = nullptr;
		vec(CComponentsII*) setCC;//temporary cc list
		DisjointSet* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
		vec(int) combineCCPos;/*record the position of ccs which need to combine with other ccs*/
		
		DisjointSetWithTime* tempDisjointSet = DBG_NEW DisjointSetWithTime(nNode);//disjoint set used for recomputed ccs
		i2iHMap /*subCCMap(nNode),*/ root2Id(nNode);//used for recomputed ccs
		int* subCCId = DBG_NEW int[nEdge];

		vec(int) maxCC;//record the position of ccs which need to be checked
		
		long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
		long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
		long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
		LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
		unordered_map<int, LinkedNode<int>*> hasChecked(nEdge);  //map the position of cc in setCC to previous node in checkedCC
		int lastTime = choiceEndT - k + 1;//last start time 
		int Te, TeShortInterval;//the end time of interval
		int edgeEndT/*, vertexNum = 0*/;
		int rightEndpoint;//max position of R set
		bool* expandMask = DBG_NEW bool[currNTimestamp];//mark certain intervals where RTM can appear (for expandable checking)

		int* tempHaveNewEdgeCC = DBG_NEW int[nNode];//have edge in R+
		int* haveNewEdgeCC = DBG_NEW int[nEdge];//have edge in R+
		CLEARALL(haveNewEdgeCC, -1, nEdge, int); 
		//int* remainEdges = DBG_NEW int[nEdge];//record remained edges after removing edges to be removed
#pragma endregion
		vector<pair<int, int>> tempRecordFromQueue;

		vec(int) * edgeSetsRShortIntv = DBG_NEW vec(int)[maxLengthForShortIntv];// R edge sets for short intervals
		veciter(int) iterEndShortIntv, iterStartShortIntv;//iterator for edgeSetsRShortIntv
		veciter(CComponentsShortIntv*) listEndShortIntv;//iterator for setCCShortIntv
		i2iHMap root2CompShortIntv(nNode);/*root2Comp for short intervals*/
		vec(CComponentsShortIntv*) setCCShortIntv;//cc list for short intervals
		DisjointSet* disjointSetShortIntv = DBG_NEW DisjointSet(nNode);//disjoint set for short intervals
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
		int saveRSize = nNode;
		
		auto checkExpandableFunc = dynamicMode ?
			&TGraph::checkExpandableFRTMMidR<CCTYPEOPT1> : &TGraph::checkExpandableFRTM<CCTYPEOPT1>;
		auto expCheckFRTMFunc = dynamicMode ? &TGraph::expCheckFRTMPlusMidR<CCTYPEOPT1> : &TGraph::expCheckFRTMPlus<CCTYPEOPT1>;

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
				CLEARALL(saveR, 0, saveRSize, bool);
				newENum = 0;
			}
			if (Ts <= lastTimeShortIntv) {//have large intervals and short intervals for [m,m+k-1]...[m,T]
				if (Ts != choiceStartT) {
					root2Comp.clear();
					disjointSet->reset();
					hasChecked.clear();
					checkedCC->release();
				}
				graph->edgeFilterPlus(Ts, Te, TeShortInterval, maxLengthForShortIntv - 1, hasE, newE, newENum, edgesAdd,
					edgeSetsR, selectedNum, edgeSetsRShortIntv, selectedNumShortIntv, rightEndpoint, newk, fixLabel, isEdgeTypeFixed, dynamicMode);
			}
			else {//only short intervals for [m,m+k-1]...[m,T]
				if (dynamicMode) {
					graph->edgeFilterShortIntvMidR(Ts, TeShortInterval, maxLengthForShortIntv - 1, lastTimeShortIntv, edgeSetsRShortIntv, selectedNumShortIntv, choiceEndT, fixLabel, isEdgeTypeFixed);
				}
				else {
					graph->edgeFilterShortIntv(Ts, TeShortInterval, maxLengthForShortIntv - 1, edgeSetsRShortIntv, selectedNumShortIntv, fixLabel, isEdgeTypeFixed);
				}
			}
			Test::compr += END_TIMER(b);
			maxSelect = max(maxSelect, selectedNum + selectedNumShortIntv);
			selectedSum += selectedNum + selectedNumShortIntv;

			if (Ts != choiceStartT) {
				disjointSetShortIntv->reset();
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
			//edgesAdd.clear();
			int tempT = choiceEndT - Te;
			for (int tempTAll = choiceEndT - TeShortInterval; tempTAll >= 0; tempTAll--, tempT--, edgeEndT--) {
				if (tempTAll >= maxLengthForShortIntv) {//large intervals
					if (tempT > rightEndpoint) continue;

					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();
					BEGIN_TIMER(c)
#pragma region maxCheck
						if (iterStart != iterEnd) {
							//BEGIN_TIMER(cc)
							graph->maintainConnectivity(iterStart, iterEnd,  disjointSet,
								combineCCPos);
							graph->combineComponentsFRTMOPT1(setCC, disjointSet, root2Comp, checkedCC, hasChecked,
								combineCCPos, haveNewEdgeCC);
							//Test::counter3 += END_TIMER(cc);
						}

					graph->updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,
						setCC, disjointSet, tempDisjointSet, root2Comp,/* &subCCMap,*/ subCCId, root2Id, checkedCC, hasChecked, 
						tempRecordFromQueue, maxCC, Ts, edgeEndT, newk, haveNewEdgeCC, tempHaveNewEdgeCC, cctype);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						graph->expCheckFRTM(maxCC, setCC, resultTF, k, Ts, 
							edgeEndT, expandMask, motifNumber, checkExpandableFunc, emptyCC);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();


				}
				else {//short intervals
					if (Ts <= lastTimeShortIntv && tempTAll == maxLengthForShortIntv - 1) {//record root of disjointSet to be saved
						graph->recordSavedRoot(setCC, hasChecked, newE, newENum, saveR, disjointSet);
					}
					iterStartShortIntv = edgeSetsRShortIntv[tempTAll].begin();
					iterEndShortIntv = edgeSetsRShortIntv[tempTAll].end();

					if (iterStartShortIntv == iterEndShortIntv) continue;
					BEGIN_TIMER(c)
#pragma region maxCheck
						addENum = 0;
					graph->maintainConnectivityPlus(iterStartShortIntv, iterEndShortIntv, addE, addENum,
						/*vertex2Pos,*/ disjointSetShortIntv, /*vertexNum,*/ hasE, saveR, disjointSet, combineCCPosShortIntv);
					graph->combineComponentsPlus(setCCShortIntv,
						disjointSetShortIntv, root2CompShortIntv, combineCCPosShortIntv);
					graph->updateNewEdgeInfoPlus(iterStartShortIntv, addE, addENum,
						setCCShortIntv, /*vertex2Pos,*/ disjointSetShortIntv,
						root2CompShortIntv, hasE, maxCCShortIntv, Ts, edgeEndT, cctype);
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
					(graph->*expCheckFRTMFunc)(maxCCShortIntv, setCCShortIntv,
						resultTF, k, Ts, edgeEndT, expandMask, motifNumber, disjointSet, hasChecked, root2Comp, 
						/*vertex2Pos,*/ setCC, emptyCC);
					Test::gne += END_TIMER(d);
					edgeSetsRShortIntv[tempTAll].clear();
				}
			}

			
			//testing
			Test::updateMemoryUse();
			//release
			graph->releaseCCs(setCC);
			graph->releaseCCs(setCCShortIntv);
		}
		if (dynamicMode) graph->updateMidResult();
		delete[] hasE;
		delete[] newE;
		delete[] addE;
		delete[] saveR;
		delete[] edgeSetsRShortIntv;
		delete disjointSetShortIntv;
		delete disjointSet;
		delete[] subCCId;
		delete[] edgesAdd;
		delete[] haveNewEdgeCC;
		delete tempDisjointSet;
		delete[] expandMask;
		delete[] edgeSetsR;
		delete[] tempHaveNewEdgeCC;
		delete checkedCC;
		cout << SELECT_EDGE << maxSelect << endl;
		cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
		cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
	}
	else {
		FRTMOpt1(graph, resultTF, fixLabel, choiceStartT, choiceEndT, dynamicMode, cctype);
	}
}


void FindRTMotif::DFRTMPLus(TGraph*& graph,
	vec(TMotifII*)*& resultTF,
	bool*& fixLabel, int graphEndT, ComponentsType cctype) {
	int endT = graph->getEndT();
	int lastRow = graphEndT - k + 1;

	reorganizeResult(graph, resultTF, graphEndT, endT);
	if (graphEndT == endT) return;

#pragma region update result
	//row number<=T-k+1
	graph->findRTMotifsDynamicPlus(k, resultTF, graphEndT, motifNumber, fixLabel, isEdgeTypeFixed, cctype);
	//row number>T-k+1 
	FRTMPlus(graph, resultTF, fixLabel, lastRow + 1, endT, true, cctype);
#pragma endregion
}

#pragma endregion


