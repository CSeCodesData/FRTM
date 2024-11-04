#include"TGraph.h"
#include"stdafx.h"

TGraph::TGraph(const TGraph& ances):
	startT(ances.startT), endT(ances.endT),
	currNTimestamp(ances.currNTimestamp),
	allNTimestamp(ances.allNTimestamp),
	nNode(ances.nNode),nEdge(ances.nEdge){
	edgeList = DBG_NEW NodePair[nEdge];
	for (int i = 0; i < nEdge; i++) {
		edgeList[i] = ances.edgeList[i];
	}
}

#pragma region construct and update the temporal graph

/* construct the temporal graph
1) graph file src: each line of the file describe an edge in a timestamp.
each line has four number: u,v,t,w (separated by ',' with no other space)
for weight w of edge (u,v) in time t. Ids of node u and v are not guaranteed to be
continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
first, and followed by edges in time 1, 2...
*/
void TGraph::constructGraph(const char* src, int k, int fixedE) {
	loadInfomation(src, fixedE, k);
}
#pragma endregion


void TGraph::createStructForDefType() {
	int size = nEdge * allNTimestamp;

	//method2
	scope = DBG_NEW pair<int, int>[nEdge];

	int eLabelSize = nEdge * numOfLabel;
	scanT = DBG_NEW int[eLabelSize];
	CLEARALL(scanT, 0, eLabelSize, int);

	vioT = DBG_NEW ForbidPairList[eLabelSize];

	CircularQueue<NVIntv>::queueSize = numOfLabel * (int)(ceil(Setting::delta*currNTimestamp) + 0.5);
	if (Setting::choice == AlgorithmType::FRTMPLUS || Setting::choice == AlgorithmType::DFRTMPLUS || Setting::choice == AlgorithmType::FRTMPLUSDYN) {
		maxIntvShortIntv = DBG_NEW pair<int, int>[nEdge];
		for (int j = 0; j < nEdge; j++) {
			maxIntvShortIntv[j].first = -1;
			maxIntvShortIntv[j].second = -1;
			scope[j].first = -1;
			scope[j].second = -1;
		}
	}

	preMaxIntv = DBG_NEW CircularQueue<NVIntv>[nEdge];

	if (Setting::choice == AlgorithmType::DFRTM || Setting::choice == AlgorithmType::FRTMFORDYN
		|| Setting::choice == AlgorithmType::DFRTMOPT1 || Setting::choice == AlgorithmType::FRTMOPT1DYN
		|| Setting::choice == AlgorithmType::DFRTMPLUS || Setting::choice == AlgorithmType::FRTMPLUSDYN) {
		newMIntR = DBG_NEW vec(MotifPos);
		newEIntR = DBG_NEW vec(MidResult*);
		newPosInEIntR = DBG_NEW int[eLabelSize];
		for (int i = 0; i < eLabelSize; i++) newPosInEIntR[i] = EMaxIntvlChange::INTVINIT;//initialize

		edgesFetched = DBG_NEW bool[nEdge];
		CLEARALL(edgesFetched, false, nEdge, bool);
		edgesInEIntR = DBG_NEW vec(int)();
	}
}

void TGraph::releaseStructForDefType() {
	
	delete[] scope;
	delete[] preMaxIntv;
	delete[] scanT;
	delete[] vioT;
		
	if (Setting::choice == AlgorithmType::FRTMPLUS|| Setting::choice == AlgorithmType::DFRTMPLUS || Setting::choice == AlgorithmType::FRTMPLUSDYN) {
		delete[] maxIntvShortIntv;
	}

	if (Setting::choice == AlgorithmType::DFRTM || Setting::choice == AlgorithmType::FRTMFORDYN
		|| Setting::choice == AlgorithmType::DFRTMOPT1 || Setting::choice == AlgorithmType::FRTMOPT1DYN
		|| Setting::choice == AlgorithmType::DFRTMPLUS || Setting::choice == AlgorithmType::FRTMPLUSDYN) {
		if (posInEIntR != nullptr) {
			delete[] posInEIntR;
		}
		delete[] newPosInEIntR;

		for (auto record : *newEIntR) {
			delete record;
		}
		delete newEIntR;

		if (MIntR != nullptr) {
			delete MIntR;
		}
		delete newMIntR;

		delete[] edgesFetched;
		delete edgesInEIntR;
	}
}

void TGraph::createCommonStruct() {
	maxIntv = DBG_NEW pair<int, int>[nEdge];
	for (int j = 0; j < nEdge; j++) {
		maxIntv[j].first = -1;
		maxIntv[j].second = -1;
	}
}

void TGraph::releaseCommonStruct() {
	delete[] maxIntv;
}



void TGraph::resetStruct() {
	for (int j = 0; j < nEdge; j++) {
		maxIntv[j].first = -1;
		maxIntv[j].second = -1;
		if (Setting::choice == AlgorithmType::DFRTMPLUS) {
			maxIntvShortIntv[j].first = -1;
			maxIntvShortIntv[j].second = -1;
		}
	}
	CLEARALL(scanT, 0, nEdge*numOfLabel, int);
	for (int j = 0; j < nEdge; j++) {
		for (int i = 0; i < numOfLabel; i++) {
			vioT[j*numOfLabel + i].release();
		}
		preMaxIntv[j].reset();
	}
}


#pragma region FRTMExact

/*generate motifs in one interval for FRTMExact*/
void TGraph::genMotifInOneIntv(SAVEINFO_VIter& iterStart, SAVEINFO_VIter& iterEnd,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, int& vertexNum,
	vec(int)& combineCCPos, /*int& realMotifNum,*/
	i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
	vec(int)& saveCCPos, int motifStartT, int motifEndT,
	vec(TMotifI*)*& result, int k, long long& motifNumber) {
	BEGIN_TIMER(gmtimer)

	#pragma region generateMaxTM
	/*maintain connectivity*/
	maintainConnectivity(iterStart, iterEnd,
		vertex2Pos, connectivity, vertexNum,
		combineCCPos);
	/*combine ccs*/
	combineComponents(tempComponents,
		vertex2Pos, connectivity, root2Comp,
		/*tempComponentsSize, realMotifNum,*/
		combineCCPos);
	/*update ccs and generate new motif*/
	updateNewEdgeInfo(iterStart, iterEnd,
		tempComponents, vertex2Pos, connectivity,
		root2Comp, saveCCPos, motifStartT, motifEndT);//O(¦¤Es)
	Test::gm += END_TIMER(gmtimer);

	#pragma endregion

	BEGIN_TIMER(gnetimer)
	expCheck(saveCCPos, tempComponents,
		result, k, motifStartT,
		motifEndT, motifNumber);
	Test::gne += END_TIMER(gnetimer);
}

#pragma region generate maximal motifs

/*fetch new edges from one R edge set and insert into the structure which maintains connectivity*/
void TGraph::maintainConnectivity(SAVEINFO_VIter& infoBegin,
	SAVEINFO_VIter& infoEnd, i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, int&vertexNum,
	vec(int)& combineCCPos) {
	NodePair* temp;
	int motifPos;
	i2iHMap_Iter vertexIt;
	int sId, tId, vertex;//vertexs'id in the disjoint set
	for (SAVEINFO_VIter iter = infoBegin;
		iter != infoEnd; ++iter) {
		temp = &edgeList[iter->edgeId];//new edge from one R edge set
		vertex = temp->first;
		vertexIt = vertex2Pos.find(vertex);//O(1)
		if (vertexIt == vertex2Pos.end()) {
			sId = vertexNum;
			vertex2Pos[vertex] = vertexNum++;
		}
		else sId = vertexIt->second;
		vertex = temp->second;
		vertexIt = vertex2Pos.find(vertex);//O(1)
		if (vertexIt == vertex2Pos.end()) {
			tId = vertexNum;
			vertex2Pos[vertex] = vertexNum++;
		}
		else tId = vertexIt->second;
		/*union-find operation means that sId and tId are connected*/ 
		//if (sId < tId) motifPos = connectivity->addE(iter->edgeId, sId, tId);
		//else motifPos = connectivity->addE(iter->edgeId, tId, sId);
		
		motifPos = connectivity->addE(iter->edgeId, sId, tId);
		//Util::output(motifPos, sId, tId);
		//cout << endl;

		if (motifPos != -1)
			combineCCPos.emplace_back(motifPos);
	}
}

inline void TGraph::combineCC(CComponents*& origin, CComponents*& now) {
	vec(TEdge)* tempEdges = &origin->edges;//edges of currentCComponents
	veciter(TEdge) edgeEnd = tempEdges->end();
	//insert edges of currentCComponents into newCComponent 
	for (auto edgeIter = tempEdges->begin();
		edgeIter != edgeEnd; ++edgeIter) {
		now->edges.emplace_back(*edgeIter);
	}

	vec(SaveCCInfo)* tempSaveInfo = &origin->saveInfo;
	veciter(SaveCCInfo) saveInfoEnd = tempSaveInfo->end();
	//insert saveInfo of currentCComponents into newCComponent 
	for (auto saveInfoIter = tempSaveInfo->begin();
		saveInfoIter != saveInfoEnd; ++saveInfoIter) {
		now->saveInfo.emplace_back(*saveInfoIter);
	}

	if (now->startT < origin->startT) {
		now->startT = origin->startT;
	}
}

/*combine ccs*/
void TGraph::combineComponents(vec(CComponents*)& tempComponents,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity,
	i2iHMap& root2Comp, vec(int)& combineCCPos) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in setCC
	CComponents* currentCComponents;//move currentCComponents's edges to newCComponent
	veciter(int) listEnd = combineCCPos.end();//setCC's iterator
	i2iHMap_Iter combineIter, mapIter;//root2Comp's iterator
	int combinePos;//the position of currentCComponents
#pragma endregion

#pragma region combination
	for (auto listIter = combineCCPos.begin();
		listIter != listEnd; ++listIter) {
		combineIter = root2Comp.find(*listIter);//cc position
		if (combineIter != root2Comp.end()) {//cc is combined to another cc
			combinePos = combineIter->second;
			currentCComponents = tempComponents[combinePos];
			
			oldRoot = currentCComponents->root;
			newRoot = connectivity->findRoot(oldRoot);
			root2Comp.erase(oldRoot);
			mapIter = root2Comp.find(newRoot);
			if (mapIter == root2Comp.end()) {//new edge will be added
				currentCComponents->root = newRoot;
				root2Comp[newRoot] = combinePos;
				continue;
			}

			combineCC(currentCComponents, tempComponents[mapIter->second]);
			
			delete currentCComponents;//currentCComponents need to be deleted
			tempComponents[combinePos] = nullptr;//this position will be deleted
		}
	}
#pragma endregion
	combineCCPos.clear();
}
/*add edges into ccs*/
void TGraph::updateNewEdgeInfo(
	SAVEINFO_VIter& infoBegin, SAVEINFO_VIter& infoEnd,
	vec(CComponents*)& tempComponents,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity,
	i2iHMap& root2Comp, vec(int)& saveCCPos, int startTime, int endTime) {

#pragma region initialize
	i2bHMap hasSaved;//means whether the motif has been saved 
	i2iHMap_Iter mapIter;//root2Comp's iterator 
	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int startT;//the edge e's start time of maxIntv[e]
	int label;//the edge e's label
	int ccPos;//cc's position in setCC
	/*the cc which has already generated/will generate*/
	CComponents* generatedCC, *newCC;
	int tempComponentsSize;//size of setCC
#pragma endregion

	for (SAVEINFO_VIter infoIter = infoBegin;
		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets

		id = infoIter->edgeId;//edge's id

		/*the root of the edge's vertex in the disjoint set*/
		root = connectivity->findRoot(vertex2Pos[edgeList[id].first]);

		mapIter = root2Comp.find(root);//check which cc the edge belongs to 
		startT = infoIter->startT;//start time of the edge
		label = infoIter->labelid; //label of the edge

		if (mapIter == root2Comp.end()) {//generateMaxTM case 1
			newCC = DBG_NEW CComponents(startT, root);
			newCC->edges.emplace_back(id, label);
			tempComponentsSize = (int)tempComponents.size();
			root2Comp[root] = tempComponentsSize;//update root2Comp
			if (startT == startTime) {//check left unexpandable
				saveCCPos.emplace_back(tempComponentsSize);
				hasSaved[tempComponentsSize] = true;
			}
			tempComponents.emplace_back(newCC);
		}
		else {//generateMaxTM case 2
			ccPos = mapIter->second;
			generatedCC = tempComponents[ccPos];

			generatedCC->edges.emplace_back(id, label);
			if (generatedCC->startT < startT) {
				generatedCC->startT = startT;
			}
			/* check whether this cc has already been inserted into save list */
			if (hasSaved.find(ccPos) == hasSaved.end()) {
				if (generatedCC->startT == startTime) {//check left unexpandable
					saveCCPos.emplace_back(ccPos);
					hasSaved[ccPos] = true;
				}
			}
		}
	}
}

void TGraph::expCheck(vec(int)&saveCCPos,
	vec(CComponents*)& tempComponents,
	vec(TMotifI*)*& result, int k, int motifStartT, int motifEndT,
	long long& motifNumber) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k), tempSavePos;
	veciter(int) saveMotifEnd = saveCCPos.end();
	TMotifI* motif;
	veciter(SaveCCInfo) saveInfoEnd;
	vec(SaveCCInfo)* saveInfoPtr;
	CComponents* tempCC;
	TMotifI* linkMotif;
	int saveInfoP;
	//newly generated connected components
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		tempCC = tempComponents[*saveMotifIter];

		motif = DBG_NEW TMotifI(motifStartT, motifEndT);
		motif->setMotifEdge(tempCC->edges);

		saveInfoPtr = &tempCC->saveInfo;
		saveInfoEnd = saveInfoPtr->end();
		SaveCCInfo saveinfo((int)result[savePos].size(), motifEndT);
		for (auto saveInfoIter = saveInfoPtr->begin();
			saveInfoIter != saveInfoEnd; ++saveInfoIter) {// O(motif number) <= O(original saved edge number)
			//reuse edges from original motifs
			tempSavePos = resultPos(motifStartT, saveInfoIter->saveEndTime, startT, endT, k);
			linkMotif = (TMotifI*)result[tempSavePos][saveInfoIter->savePos];
			saveInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
			if (saveInfoIter->saveEndTime == endT) {
				saveinfo.saveCCInfoPos = saveInfoP;
				linkMotif->tempLinkToMotifs(saveinfo);
			}
		}

		tempCC->edges.clear();
		tempCC->saveInfo.clear();
		tempCC->saveInfo.emplace_back(saveinfo);

		result[savePos].emplace_back(motif);
	}
	motifNumber += saveCCPos.size();
	saveCCPos.clear();
}
#pragma endregion 

#pragma endregion


void TGraph::updateMidResult() {
	/*auto it = newValidMidResult.begin();
	for (auto be = newEIntR->begin(); be != newEIntR->end(); be++,it++) {
		if(*it)
		cout << (*be)->row << " " << (*be)->edgeId << endl;
	}
	exit(0);*/
	
	if (posInEIntR != nullptr) {
		delete[] posInEIntR;
	}
	posInEIntR = newPosInEIntR;
		
	int eLabelSize = nEdge* numOfLabel;
	newPosInEIntR = DBG_NEW int[eLabelSize];
	for (int i = 0; i < eLabelSize; i++) newPosInEIntR[i] = EMaxIntvlChange::INTVINIT;
		
	if (EIntR != nullptr) {
		delete EIntR;
	}
	EIntR = newEIntR;

	if (MIntR != nullptr) {
		delete MIntR;
	}
	MIntR = newMIntR;

	newMIntR = DBG_NEW vec(MotifPos);
	newEIntR = DBG_NEW vec(MidResult*);
	validMidResult.swap(newValidMidResult);
	newValidMidResult.clear();
	cout << "EIntR: " << EIntR->size() << " MIntR: " << MIntR->size()<<endl;
	
}

void TGraph::showMidResult(vec(TMotifII*)*& result) {
	int edgeN = 0;
	unordered_set<int> s;
	for (auto item : *EIntR) {
		if (s.find(item->edgeId) == s.end()) {
			s.emplace(item->edgeId);
			edgeN++;
		}
	}
	cout << " Edge in EIntR: " << edgeN;
	edgeN = 0;
	for (auto item : *MIntR) {
		edgeN += result[item.first][item.second]->getSize();
	}
	cout << " Edge in MIntR: " << edgeN << endl;
	//exit(0);
}


/*print motif*/
void TGraph::printMotif(TMotifII*& motif, int motifId, int node1, int node2) {
	vec(int)* motifEdge = motif->getMotifEdge();
	veciter(int) listEnd = motifEdge->end();
	sort(motifEdge->begin(), listEnd);
	if (node1 == -1 || node2 == -1) {
		cout << MOTIF_ID << motifId << endl;
		cout << EDGE_NUM << motif->getSize() << endl;
		//vec(int)* maskEdge = motif->getMaskEdge();

		NodePair *edge;
		//bool finded;
		for (auto iter = motifEdge->begin();
			iter != listEnd; ++iter) {
			//finded = Util::findItem(maskEdge, iter->id);
			//if (!finded) {
			edge = &edgeList[*iter];
			cout << edge->first << "," << edge->second <<
				"," << idToLabel[getEdgeLabel(*iter, motif->getStartT())] << endl;
			//}
		}
	}
	else {
		bool output = false;
		NodePair *edge;
		for (auto iter = motifEdge->begin();
			iter != listEnd; ++iter) {
			edge = &edgeList[*iter];
			if ((edge->first == node1 && edge->second == node2) || (edge->first == node2 && edge->second == node1)) {
				output = true;
				break;
			}
		}
		if (output) {
			cout << "***************************************" << endl;
			cout << MOTIF_ID << motifId << endl;
			cout << EDGE_NUM << motif->getSize() << endl;
			for (auto iter = motifEdge->begin();
				iter != listEnd; ++iter) {
				edge = &edgeList[*iter];
				cout << edge->first << "," << edge->second <<
					"," << idToLabel[getEdgeLabel(*iter, motif->getStartT())] << endl;
			}
			cout << "***************************************" << endl;
		}
	}
}

