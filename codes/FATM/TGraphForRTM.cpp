#pragma once
#include"TGraph.h"

#pragma region maintain connectivity in maxCheck for FRTM and FRTMOPT1
void TGraph::connected(NodePair* temp, int edgeId, DisjointSet*& connectivity, /*i2iHMap& vertex2Pos,*/ vec(int)& combineCCPos/*, int&vertexNum*/) {
	int motifPos = connectivity->addE(temp->first/*sId*/, temp->second/*tId*/);

	if (motifPos != -1)
		combineCCPos.emplace_back(motifPos);
}

/*fetch new edges from one R edge set and insert into the structure which maintains connectivity*/
void TGraph::maintainConnectivity(veciter(int)& infoBegin,
	veciter(int)& infoEnd, /*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivity, /*int&vertexNum,*/
	vec(int)& combineCCPos) {
	for (auto iter = infoBegin;
		iter != infoEnd; ++iter) {
		connected(&edgeList[*iter], *iter, connectivity, /*vertex2Pos,*/ combineCCPos/*, vertexNum*/);
	}
}

#pragma endregion

#pragma region maxCheck for FRTM

void TGraph::combineComponentsFRTM(vec(CComponentsII*)& setCC,
	DisjointSet*& connectivity,
	i2iHMap& root2Comp, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
	vec(int)& combineCCPos) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in setCC
	CComponentsII* currentCComponents;//move currentCComponents's edges to newCComponent

	veciter(TEdge) edgeEnd;//tempEdges' iterator
	veciter(int) vecIntEnd;
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
			currentCComponents = setCC[combinePos];
			oldRoot = currentCComponents->root;
			newRoot = connectivity->DisjointSet::findRoot(oldRoot);

			root2Comp.erase(oldRoot);
			mapIter = root2Comp.find(newRoot);
			if (mapIter == root2Comp.end()) {//new edge will be added
				currentCComponents->root = newRoot;
				root2Comp[newRoot] = combinePos;
				continue;
			}

			//update haschecked
			auto preNode = hasChecked.find(combinePos);
			if (preNode != hasChecked.end()) {
				if (preNode->second == nullptr) {
					checkedCC->deleteFirstNode();
					if (checkedCC->first != nullptr) hasChecked[checkedCC->first->item] = nullptr;
				}
				else {
					checkedCC->deleteNextNode(preNode->second);
					auto nextNode = preNode->second->next;
					if (nextNode != nullptr) hasChecked[nextNode->item] = preNode->second;
				}
				hasChecked.erase(combinePos);
			}

			setCC[mapIter->second]->combineFrom(currentCComponents);

			delete currentCComponents;//currentCComponents need to be deleted
			setCC[combinePos] = nullptr;//this position will be deleted
		}
	}
#pragma endregion
	combineCCPos.clear();
}

void TGraph::combineComponentsFRTMOPT1(vec(CComponentsII*)& setCC,
	DisjointSet*& connectivity,
	i2iHMap& root2Comp, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
	vec(int)& combineCCPos, int*& haveNewEdgeCC) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in setCC
	CComponentsII* currentCComponents;//move currentCComponents's edges to newCComponent

	veciter(TEdge) edgeEnd;//tempEdges' iterator
	veciter(int) vecIntEnd;
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
			currentCComponents = setCC[combinePos];
			oldRoot = currentCComponents->root;
			newRoot = connectivity->DisjointSet::findRoot(oldRoot);

			root2Comp.erase(oldRoot);
			mapIter = root2Comp.find(newRoot);
			if (mapIter == root2Comp.end()) {//new edge will be added
				currentCComponents->root = newRoot;
				root2Comp[newRoot] = combinePos;
				continue;
			}

			//update haschecked
			auto preNode = hasChecked.find(combinePos);
			if (preNode != hasChecked.end()) {
				if (preNode->second == nullptr) {
					checkedCC->deleteFirstNode();
					if (checkedCC->first != nullptr) hasChecked[checkedCC->first->item] = nullptr;
				}
				else {
					checkedCC->deleteNextNode(preNode->second);
					auto nextNode = preNode->second->next;
					if (nextNode != nullptr) hasChecked[nextNode->item] = preNode->second;
				}
				hasChecked.erase(combinePos);
			}

			setCC[mapIter->second]->combineFrom(currentCComponents);
			haveNewEdgeCC[mapIter->second] = max(haveNewEdgeCC[mapIter->second], haveNewEdgeCC[combinePos]);
			
			delete currentCComponents;//currentCComponents need to be deleted
			setCC[combinePos] = nullptr;//this position will be deleted
		}
	}
#pragma endregion
	combineCCPos.clear();
}

void TGraph::updateCCFRTM(vec(CComponentsII*)& setCC, i2iHMap&  root2Comp,
	LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,/* int& realMotifNum,*/
	veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsType type) {
	i2iHMap_Iter mapIter = root2Comp.find(root);//check which cc the edge belongs to 
	ForbidPairNode* now;
	int label = getEdgeLabel(id, startTime - startT);
	if (mapIter == root2Comp.end()) {//generateMaxTM case 1
		auto newCC = (CComponentsFRTM*)CComponentsIIFactory::instance(root, type);
		newCC->addEdge(id);

		scope[id] = maxIntv[id];
		newCC->minEMaxIntvlEndT = scope[id].second;
		newCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		if (Setting::delta != 0 && Setting::c != 0) {
			if (!preMaxIntv[id].empty()) {
				auto intvf = preMaxIntv[id].front;
				while (intvf != preMaxIntv[id].rear) {
					auto& intv = preMaxIntv[id].q[intvf];
					if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
						newCC->preEMaxIntvlChangePos->emplace(intv.first, intv.second, id);
					}
					intvf++;
				}
			}

			now = vioT[numOfLabel * id + label].first;
			while (now != nullptr) {//O(delta log delta)
				int intvStartT, intvEndT;
				tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
				if (intvStartT > endTime) break;
				else if (intvStartT >= filterTime && intvStartT <= endTime) {
					newCC->noisePosQueue->emplace(
						min(intvEndT, endTime), // noise start position
						NoisePos(
							intvStartT, // noise end position
							0 // edge position in newCC->edges
						)

					);
				}
				now = now->next;
			}
		}
		
		newCC->tabuTChangePos = -1;

		int setCCSize = (int)setCC.size();
		setCC.emplace_back(newCC);
		root2Comp[root] = setCCSize;//update root2Comp

		//need checked
		hasChecked[setCCSize] = checkedCC->tail;
		checkedCC->addItemAtLast(setCCSize);
	}
	else {//generateMaxTM case 2
		int ccPos = mapIter->second;
		auto generatedCC = (CComponentsFRTM*)setCC[ccPos];
		int edgeSize = (int)generatedCC->edges.size();
		generatedCC->addEdge(id);

		scope[id] = maxIntv[id];
		if (generatedCC->minEMaxIntvlEndT > scope[id].second)
			generatedCC->minEMaxIntvlEndT = scope[id].second;
		generatedCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		if (Setting::delta != 0 && Setting::c != 0) {
			if (!preMaxIntv[id].empty()) {
				auto intvf = preMaxIntv[id].front;
				while (intvf != preMaxIntv[id].rear) {
					auto& intv = preMaxIntv[id].q[intvf];
					if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
						generatedCC->preEMaxIntvlChangePos->emplace(intv.first, intv.second, id);
					}
					intvf++;
				}
			}

			now = vioT[numOfLabel * id + label].first;
			while (now != nullptr) {//O(delta log delta)
				int intvStartT, intvEndT;
				tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
				if (intvStartT > endTime) break;
				else if (intvStartT >= filterTime && intvStartT <= endTime) {
					generatedCC->noisePosQueue->emplace(
						min(intvEndT, endTime), // noise start position
						NoisePos(
							intvStartT, // noise end position
							edgeSize // edge position in newCC->edges
						)
					);
				}
				now = now->next;
			}
		}
		generatedCC->tabuTChangePos = -1;

		//need checked
		if (hasChecked.find(ccPos) == hasChecked.end()) {
			hasChecked[ccPos] = checkedCC->tail;
			checkedCC->addItemAtLast(ccPos);
		}
	}

}

void TGraph::maintainTempCCForUFWithTime(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
	DisjointSet* connectivity, i2iHMap& root2Comp, /*i2iHMap* subCCMap,*/ int*& subCCId, i2iHMap& root2Id, CComponentsFRTM* generatedCC, int currentTime) {
	//BEGIN_NSTIMER(p);

	DisjointSetWithTime* ds = (DisjointSetWithTime*)connectivity;

	//subCCMap->clear();
	root2Id.clear();
	//get the number of new connected components
	auto tempEdges = &generatedCC->edges;
	int size = (int)tempEdges->size();
	//int ufsetSize = min((int)size << 1, nNode);
	//ds->reset(ufsetSize);
	veciter(int) motifEdgesEnd = tempEdges->end();
	auto subCCIter = &subCCId[0];
	//int sId, tId;
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter) {//O(Em)
		if (*subCCIter != -1) {
			auto tempEdge = &edgeList[*iter];

			//int vertex = tempEdge->first;
			//auto vertexIt = subCCMap->find(vertex);//O(1)
			//sId = vertexIt->second;
			//vertex = tempEdge->second;
			//vertexIt = subCCMap->find(vertex);//O(1)
			//tId = vertexIt->second;

			/*union-find operation means that sId and tId are connected*/
			//ds->addE( sId, tId, currentTime);
			ds->addE(tempEdge->first, tempEdge->second, currentTime);
		}
	}
	//BEGIN_NSTIMER(q);
	//get subCCNum
	int subCCNum = 0;
	subCCIter = &subCCId[0];
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter) {//O(Em)
		if (*subCCIter != -1) {
			/*the root of the edge's vertex in the disjoint set*/
			//int root = ds->findRoot((*subCCMap)[edgeList[*iter].first]);
			int root = ds->findRoot(edgeList[*iter].first);
			auto ccIt = root2Id.find(root);
			if (ccIt == root2Id.end()) {//new subCC
				root2Id[root] = subCCNum;
				*subCCIter = subCCNum;
				++subCCNum;
			}
			else {
				*subCCIter = ccIt->second;
			}
		}
	}
	//Test::counter2 += END_NSTIMER(p);

	//get subCCs
	if (generatedCC->subCCs == nullptr) {
		generatedCC->subCCs = DBG_NEW vec(int)[subCCNum];
		generatedCC->subCCsaved = DBG_NEW bool[subCCNum];
	}
	else {
		generatedCC->subCCNum = (int)_msize(generatedCC->subCCsaved) / sizeof(bool);
		if (generatedCC->subCCNum < subCCNum) {
			delete[] generatedCC->subCCs;
			delete[] generatedCC->subCCsaved;
			generatedCC->subCCs = DBG_NEW vec(int)[subCCNum];
			generatedCC->subCCsaved = DBG_NEW bool[subCCNum];
		}
		else {
			for (int i = 0; i < subCCNum; i++) {
				generatedCC->subCCs[i].clear();
			}
		}
	}
	generatedCC->subCCNum = subCCNum;

	//BEGIN_NSTIMER(a);
	subCCIter = &subCCId[0];
	int edgePos = 0, ccid;
	for (auto iter = generatedCC->edges.begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
		ccid = *subCCIter;
		if (ccid != -1) {
			generatedCC->subCCs[ccid].emplace_back(*iter);
			generatedCC->subCCsaved[ccid] = edgePos >= generatedCC->newInsert;//right non-expandable
		}
	}
	//Test::counter3 += END_NSTIMER(a);
}


void TGraph::updateNewEdgeInfoFRTM(
	veciter(int)& infoBegin, veciter(int)& infoEnd,
	vec(CComponentsII*)& setCC,
	/*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivity, DisjointSet* tempDisjointSet,
	i2iHMap& root2Comp, /*i2iHMap* subCCMap,*/ int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue,
	vec(int)& saveCCPos, int startTime, int endTime, int k, ComponentsType type) {
	if (checkedCC->first == nullptr && infoBegin == infoEnd) return;

	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int filterTime = startTime + k - 1;

	//BEGIN_TIMER(a)
	for (auto infoIter = infoBegin;
		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets

		id = *infoIter;//edge's id

		/*the root of the edge's vertex in the disjoint set*/
		root = connectivity->findRoot(edgeList[id].first);
		//root = connectivity->findRoot(vertex2Pos[edgeList[id].first]);
		updateCCFRTM(setCC, root2Comp, checkedCC, hasChecked, infoIter,
			id, root, filterTime, startTime, endTime, k, type);
	}
	//Test::counter += END_TIMER(a);

	LinkedNode<int>* tempCCIter = checkedCC->first;
	if (tempCCIter == nullptr) return;
	int size;
	int ccPos;//cc's position in setCC
	/*the cc which has already generated/will generate*/
	veciter(int) motifEdgesEnd;
	veciter(int) subCCIter;

	i2iHMap_Iter vertexIt, ccIt;
	int removeEdge = 0;
	vec(int)* tempEdges;
	int currentTime = startTime * currNTimestamp + (endT - endTime);
	while (tempCCIter != nullptr) {
		ccPos = tempCCIter->item;
		auto generatedCC = (CComponentsFRTM*)setCC[ccPos];
		tempEdges = &generatedCC->edges;
		size = (int)tempEdges->size();
		if (generatedCC->newInsert == size) {//not updated edges
			if (generatedCC->tabuTChangePos != -1 && generatedCC->tabuTChangePos < endTime) {
				tempCCIter = tempCCIter->next;
				continue;
			}
		}
		
		if (Setting::delta != 0 && Setting::c != 0) {
			removeEdgesAndUpdateScopes(generatedCC, subCCId, tempRecordFromQueue, endTime, removeEdge);
		}

		//BEGIN_TIMER(c)
		if (removeEdge > 0) {
			if (removeEdge == size) {
				tempCCIter = tempCCIter->next;
				generatedCC->newInsert = size;
				continue;
			}
			//Test::counter3 += tempEdges->size();

			maintainTempCCForUFWithTime(tempCCIter, tempRecordFromQueue, tempDisjointSet, root2Comp, /*&vertex2Pos,*/ subCCId, root2Id, generatedCC, currentTime);

			saveCCPos.emplace_back(ccPos);
			tempCCIter = tempCCIter->next;
		}
		else {//removeEdge == 0: not edges removed
			generatedCC->subCCNum = 0;
			auto preNode = hasChecked[ccPos];
			if (preNode == nullptr) {
				checkedCC->deleteFirstNode();
				tempCCIter = checkedCC->first;
				if (tempCCIter != nullptr) hasChecked[tempCCIter->item] = nullptr;
			}
			else {
				checkedCC->deleteNextNode(preNode);
				auto nextNode = preNode->next;
				if (nextNode != nullptr) hasChecked[nextNode->item] = preNode;
				tempCCIter = nextNode;
			}
			hasChecked.erase(ccPos);
			saveCCPos.emplace_back(ccPos);
		}
		//Test::counter3 += END_TIMER(c);
	}
}

bool TGraph::checkLeftOrBothExpandable(vec(int)& edges, int motifStartT, int motifEndT,
	pair<int, int>& intersectIntv, bool*& expandMask) {
	
	if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
		int j = intersectIntv.second;
		for (; j >= motifEndT; --j) {
			CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
			if (bothFitDefAndSameLabelChangeStartTimePos(edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
				return false;
			}
		}
	}
	else {
		int j = intersectIntv.first;
		for (; j <= motifStartT - 1; ++j) {
			CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
			if (bothFitDefAndSameLabelChangeEndTimePos(edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
				return false;
			}
		}
	}
	return true;
}

bool TGraph::checkRightOrLeftOrBothExpandable(vec(int)& edges, int motifStartT, int motifEndT,
	pair<int, int>& intersectIntv, bool*& expandMask) {
	if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
		int j = intersectIntv.second;
		for (; j > motifEndT; --j) {
			CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
			if (bothFitDefAndSameLabelChangeStartTimePos(edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
				return false;
			}
		}
	}
	else {
		int j = intersectIntv.first;
		for (; j <= motifStartT; ++j) {
			CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
			if (bothFitDefAndSameLabelChangeEndTimePos(edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
				return false;
			}
		}
	}

	CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
	if (bothFitDefAndSameLabelChangeStartTimePos(edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
		return false;
	}
	return true;
}

bool TGraph::checkUpdateInMIntR(vec(int)& edges, veciter(int)& edgeEnd, int motifStartT) {
	for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
		if (newPosInEIntR[*edgeIter * numOfLabel + getEdgeLabel(*edgeIter, motifStartT)] == EMaxIntvlChange::INTVINIT) {
			return false;
		}
	}
	return true;
}

#pragma endregion

#pragma region maxCheck for FRTMOpt1 

void TGraph::updateCCOpt1(vec(CComponentsII*)& setCC, i2iHMap&  root2Comp, 
	/*unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,*/int*& edgesAdd, int*& haveNewEdgeCC,
	LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsType ccstype) {
	i2iHMap_Iter mapIter = root2Comp.find(root);//check which cc the edge belongs to 
	ForbidPairNode* now;
	int label = getEdgeLabel(id, startTime - startT);
	if (mapIter == root2Comp.end()) {//generateMaxTM case 1
			auto newCC = (CCTYPEOPT1*)CComponentsIIFactory::instance(root, ccstype);
		newCC->addEdge(id);

		scope[id] = maxIntv[id];
		newCC->minEMaxIntvlEndT = scope[id].second;
		newCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);

		if (Setting::delta != 0 && Setting::c != 0) {
			if (!preMaxIntv[id].empty()) {
				auto intvf = preMaxIntv[id].front;
				while (intvf != preMaxIntv[id].rear) {
					auto& intv = preMaxIntv[id].q[intvf];
					if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
						newCC->preEMaxIntvlChangePos->emplace(intv.first, intv.second, id);
					}
					intvf++;
				}
			}

			now = vioT[numOfLabel * id + label].first;
			while (now != nullptr) {//O(delta log delta)
				int intvStartT, intvEndT;
				tie(intvStartT, intvEndT) = now->item;
				if (intvStartT > endTime) break;
				else if (intvStartT >= filterTime && intvStartT <= endTime) {
					newCC->noisePosQueue->emplace(
						min(intvEndT, endTime), // noise start position
						NoisePos(
							intvStartT, // noise end position
							0 // edge position in newCC->edges
						)

					);
				}
				now = now->next;
			}
		}
		

		newCC->tabuTChangePos = -1;

		int setCCSize = (int)setCC.size();
		setCC.emplace_back(newCC);
		root2Comp[root] = setCCSize;//update root2Comp

		//if (edgesAdd.find(id) != edgesAddEnd) {
		if (edgesAdd[id] == startTime) {
			haveNewEdgeCC[setCCSize] = startTime;
			//need checked
			hasChecked[setCCSize] = checkedCC->tail;
			checkedCC->addItemAtLast(setCCSize);
		}
	}
	else {//generateMaxTM case 2
			int ccPos = mapIter->second;
		auto generatedCC = (CCTYPEOPT1*)setCC[ccPos];
		int edgeSize = (int)generatedCC->edges.size();
		generatedCC->addEdge(id);

		scope[id] = maxIntv[id];
		if (generatedCC->minEMaxIntvlEndT > scope[id].second)
			generatedCC->minEMaxIntvlEndT = scope[id].second;
		generatedCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		if (Setting::delta != 0 && Setting::c != 0) {
			if (!preMaxIntv[id].empty()) {
				auto intvf = preMaxIntv[id].front;
				while (intvf != preMaxIntv[id].rear) {
					auto& intv = preMaxIntv[id].q[intvf];
					if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
						generatedCC->preEMaxIntvlChangePos->emplace(intv.first, intv.second, id);
					}
					intvf++;
				}
			}

			now = vioT[numOfLabel * id + label].first;
			while (now != nullptr) {//O(delta log delta)
				int intvStartT, intvEndT;
				tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
				if (intvStartT > endTime) break;
				else if (intvStartT >= filterTime && intvStartT <= endTime) {
					generatedCC->noisePosQueue->emplace(
						min(intvEndT, endTime), // noise start position
						NoisePos(
							intvStartT, // noise end position
							edgeSize // edge position in newCC->edges
						)
					);
				}
				now = now->next;
			}
		}
		generatedCC->tabuTChangePos = -1;


		if (hasChecked.find(ccPos) == hasChecked.end()) {
			if (haveNewEdgeCC[ccPos] == startTime) {
				hasChecked[ccPos] = checkedCC->tail;
				checkedCC->addItemAtLast(ccPos);
			}
			else if (edgesAdd[id] == startTime) {
				haveNewEdgeCC[ccPos] = startTime;
				hasChecked[ccPos] = checkedCC->tail;
				checkedCC->addItemAtLast(ccPos);
			}
		}
	}
}

void TGraph::updateNewEdgeInfoOpt1(
	veciter(int)& infoBegin, veciter(int)& infoEnd, /* unordered_set<int>& edgesAdd,*/int*& edgesAdd,
	vec(CComponentsII*)& setCC,
	/*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivity, DisjointSet* tempDisjointSet,
	i2iHMap& root2Comp, int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& checkedCC, 
	unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue,
	vec(int) & saveCCPos, int startTime, int endTime, int k, int*& haveNewEdgeCC, int*& tempHaveNewEdgeCC, ComponentsType ccstype) {
	if (checkedCC->first == nullptr && infoBegin == infoEnd) return;

	//BEGIN_TIMER(a)
	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int filterTime = startTime + k - 1;
	//auto addEnd = edgesAdd.end();
	for (auto infoIter = infoBegin;
		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets

		id = *infoIter;//edge's id
		//cout << startTime<<" "<< endTime << " " << id << endl;
		/*the root of the edge's vertex in the disjoint set*/
		root = connectivity->findRoot(edgeList[id].first);
		updateCCOpt1(setCC, root2Comp, edgesAdd, haveNewEdgeCC, /*addEnd,*/ checkedCC, hasChecked,
			id, root, filterTime, startTime, endTime, k, ccstype);
	}
	//Test::counter += END_TIMER(a);


	LinkedNode<int>* tempCCIter = checkedCC->first;
	if (tempCCIter != nullptr) {

	

	int size;
	int ccPos;//cc's position in setCC
	/*the cc which has already generated/will generate*/
	veciter(int) motifEdgesEnd;
	veciter(int) subCCIter;

	i2iHMap_Iter vertexIt, ccIt;
	int removeEdge = 0;
	vec(int)* tempEdges;
	int currentTime = startTime * currNTimestamp + (endT - endTime);
	while (tempCCIter != nullptr) {

		//BEGIN_TIMER(b)

		ccPos = tempCCIter->item;
		auto generatedCC = (CCTYPEOPT1*)setCC[ccPos];
		tempEdges = &generatedCC->edges;

		size = (int)tempEdges->size();
		/*if (endTime == 153 && ccPos == 638) {
			exit(0);
		}*/
		if (generatedCC->newInsert == size) {//not updated edges
			if (generatedCC->tabuTChangePos != -1 && generatedCC->tabuTChangePos < endTime) {
				tempCCIter = tempCCIter->next;
				continue;
			}
		}

		//Test::counter++;
		//Test::counter2 += size;
		
		//BEGIN_TIMER(b)
		if (Setting::delta != 0 && Setting::c != 0) {
			removeEdgesAndUpdateScopes(generatedCC, subCCId, tempRecordFromQueue, endTime, removeEdge);
		}
		
		//cout << startTime << " " << endTime << " " << size << " " << removeEdge << endl;
		//BEGIN_TIMER(c)
		if (removeEdge > 0) {
			if (removeEdge == size) {
				tempCCIter = tempCCIter->next;
				generatedCC->newInsert = size;
				continue;
			}
			//Test::counter3 += tempEdges->size();

			if (maintainTempCCForUFOpt1WithTime(tempCCIter, tempRecordFromQueue, edgesAdd, /*addEnd,*/ tempDisjointSet,
				root2Comp, /*&vertex2Pos,*/ subCCId, root2Id, generatedCC, startTime, endTime, removeEdge,
				tempHaveNewEdgeCC, /*remainEdges,*/ currentTime)) {
				saveCCPos.emplace_back(ccPos);
			}
			else {
				generatedCC->newInsert = size;
			}
			tempCCIter = tempCCIter->next;
			
			/*if (generatedCC->subCCNum != 1) {
				Test::counter4++;
			}
			else {
				Test::counter5++;
			}
			Test::counter3++;*/
		}
		else {//removeEdge == 0: not edges removed
			generatedCC->subCCNum = 0;
			auto preNode = hasChecked[ccPos];
			if (preNode == nullptr) {
				checkedCC->deleteFirstNode();
				tempCCIter = checkedCC->first;
				if (tempCCIter != nullptr) hasChecked[tempCCIter->item] = nullptr;
			}
			else {
				checkedCC->deleteNextNode(preNode);
				auto nextNode = preNode->next;
				if (nextNode != nullptr) hasChecked[nextNode->item] = preNode;
				tempCCIter = nextNode;
			}
			hasChecked.erase(ccPos);
			saveCCPos.emplace_back(ccPos);
		}
		//Test::counter6 += END_TIMER(c);
	}
	}
}


#pragma endregion


#pragma region FRTMPlus
void TGraph::maintainConnectivityPlus(veciter(int)& infoBegin,
	veciter(int)& infoEnd, int*& addE, int& addENum, /*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivityShortIntv, /*int&vertexNum,*/
	bool*& hasE, bool*& saveR, DisjointSet*& connectivity,
	vec(int)& combineCCPos) {
	NodePair* temp;
	for (veciter(int) iter = infoBegin;
		iter != infoEnd; ++iter) {
		temp = &edgeList[*iter];//new edge from one R edge set
		if (hasE[*iter]) {
			//assert(vertex2Pos.find(vertex) != vertex2Pos.end());
			//int r = connectivity->findRoot(vertex2Pos[temp->first]);
			int r = connectivity->findRoot(temp->first);
			if (saveR[r]) {
				addE[addENum++] = (int)(iter - infoBegin);
				connected(temp, *iter, connectivityShortIntv, /*vertex2Pos,*/ combineCCPos/*, vertexNum*/);
			}
		}
		else {
			addE[addENum++] = (int)(iter - infoBegin);
			connected(temp, *iter, connectivityShortIntv, /*vertex2Pos,*/ combineCCPos/*, vertexNum*/);
		}
	}
}

void TGraph::combineComponentsPlus(vec(CComponentsShortIntv*)& tempComponents,
	DisjointSet*& connectivity,
	i2iHMap& root2Comp, vec(int)& combineCCPos) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in setCC
	CComponentsShortIntv* currentCComponents;//move currentCComponents's edges to newCComponent
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

			tempComponents[mapIter->second]->combineFrom(currentCComponents);

			delete currentCComponents;//currentCComponents need to be deleted
			tempComponents[combinePos] = nullptr;//this position will be deleted
		}
	}
#pragma endregion
	combineCCPos.clear();
}


void TGraph::updateNewEdgeInfoPlus(
	veciter(int)& infoBegin,
	int*& addE, int addENum,
	vec(CComponentsShortIntv*)& setCC,
	/*i2iHMap& vertex2PosShortIntv,*/ DisjointSet*& connectivity,
	i2iHMap& root2Comp, bool*& hasE,
	vec(int)& maxCC, int startTime, int endTime, ComponentsType type) {

#pragma region initialize
	i2bHMap hasSaved;//means whether the motif has been saved 
	i2iHMap_Iter mapIter;//root2Comp's iterator 
	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int eStartT;//the edge e's start time of maxIntv[e]
	int label;//the edge e's label
	int ccPos;//cc's position in setCC
	/*the cc which has already generated/will generate*/
	CComponentsShortIntv* generatedCC, *newCC;
	int tempComponentsSize;//size of setCC
	int tempScopeL;
#pragma endregion
	bool nonexpand;
	//auto vEnd = vertex2PosShortIntv.end();
	int count = 0;
	for (int addEIter = 0; addEIter < addENum; ++addEIter) {//new edges in R edge sets
		veciter(int) infoIter = infoBegin + addE[addEIter];
		id = *infoIter;//edge's id

		/*the root of the edge's vertex in the disjoint set*/
		//root = connectivity->findRoot(vertex2PosShortIntv[edgeList[id].first]);
		root = connectivity->findRoot(edgeList[id].first);

		mapIter = root2Comp.find(root);//check which cc the edge belongs to 
		eStartT = maxIntvShortIntv[id].first;//start time of the edge
		label = getEdgeLabel(id, eStartT); //label of the edge

		nonexpand = !hasE[id];

		if (mapIter == root2Comp.end()) {//generateMaxTM case 1
			newCC = CComponentsIIFactory::instanceShortIntv(eStartT, root, type);
			newCC->addEdge(id);
			newCC->haveNonExpand = nonexpand;
			if (maxIntv[id].second >= endTime) {
				newCC->scopeL = maxIntv[id].first;
				newCC->scopeR = maxIntv[id].second;
				if (!preMaxIntv[id].empty()) {
					auto intvf = preMaxIntv[id].front;
					while (intvf != preMaxIntv[id].rear) {
						auto& intv = preMaxIntv[id].q[intvf];
						if (intv.second >= endTime && intv.first < newCC->scopeL  &&getEdgeLabel(id, intv.second - startT) == label) {//need same label
							newCC->scopeL = intv.first;
						}
						intvf++;
					}
				}
			}
			else {
				newCC->scopeL = newCC->scopeR = -1;
			}


			tempComponentsSize = (int)setCC.size();
			root2Comp[root] = tempComponentsSize;//update root2Comp
			maxCC.emplace_back(tempComponentsSize);
			hasSaved[tempComponentsSize] = true;
			setCC.emplace_back(newCC);
		}
		else {//generateMaxTM case 2
			ccPos = mapIter->second;
			generatedCC = setCC[ccPos];
			generatedCC->haveNonExpand |= nonexpand;
			generatedCC->addEdge(id);

			if (generatedCC->scopeR != -1) {
				if (maxIntv[id].second >= endTime) {
					if (generatedCC->scopeR > maxIntv[id].second)
						generatedCC->scopeR = maxIntv[id].second;
					tempScopeL = maxIntv[id].first;
					if (!preMaxIntv[id].empty()) {
						auto intvf = preMaxIntv[id].front;
						while (intvf != preMaxIntv[id].rear) {
							auto& intv = preMaxIntv[id].q[intvf];
							if (intv.second >= endTime && intv.first < tempScopeL  &&getEdgeLabel(id, intv.second - startT) == label) {//need same label
								tempScopeL = intv.first;
							}
							intvf++;
						}
					}
					if (generatedCC->scopeL < tempScopeL)
						generatedCC->scopeL = tempScopeL;
				}
				else {
					generatedCC->scopeL = generatedCC->scopeR = -1;
				}
			}

			if (generatedCC->startT < eStartT) {
				generatedCC->startT = eStartT;
			}
			/* check whether this cc has already been inserted into save list */
			if (hasSaved.find(ccPos) == hasSaved.end()) {
				maxCC.emplace_back(ccPos);
				hasSaved[ccPos] = true;
			}
		}
	}
}

void TGraph::recordSavedRoot(vec(CComponentsII*)& setCC, unordered_map<int, LinkedNode<int>*>& hasChecked, int*& newE, int& newENum, bool*& saveR, DisjointSet*& connectivity/*, i2iHMap& vertex2Pos*/) {
	for (auto cc : hasChecked) {
		saveR[setCC[cc.first]->root] = true;
	}
	//auto vEnd = vertex2Pos.end();
	for (auto i = 0; i < newENum; i++) {
		int e = newE[i];
		//auto v = vertex2Pos.find(edgeList[e].first);
		//if (v != vEnd) {
			//int r = connectivity->findRoot(v->second);
			int r = connectivity->findRoot(edgeList[e].first);
			if (r != -1) {
				saveR[r] = true;
			}
		//}
		//v = vertex2Pos.find(edgeList[e].second);
		//if (v != vEnd) {
			//int r = connectivity->findRoot(v->second);
			r = connectivity->findRoot(edgeList[e].second);
			if (r != -1) {
				saveR[r] = true;
			}
		//}
	}
}

#pragma endregion

void TGraph::findRTMotifsDynamic(int k,
	vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed, ComponentsType cctype) {
#pragma region intialization
	//i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	CComponentsFRTM* emptyCC = nullptr;
	vec(CComponentsII*) setCC;//temporary component list CC
	vec(int)* edgeSetsR = DBG_NEW vec(int)[endT - oriEndT + 1];// R edge sets
	DisjointSet* disjointSet = DBG_NEW DisjointSet(nNode + 1);//disjoint set

	i2iHMap root2Id;
	DisjointSetWithTime* tempDisjointSet = DBG_NEW DisjointSetWithTime(nNode);
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	veciter(CComponentsII*) listEnd;//iterator for setCC
	int selectedNum;//the number of edges in edgeSetsR

	vec(int) saveMotifPos;//record the position of components which need to be saved
	/*record the position of components which need to combine with other components*/
	vec(int) combineMotifPos;
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	long long allSmk = 0, maxSmk = 0, smk;//the sum and the max number of S[m,m+k-1]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked; // map ccpos in setCC to previous node in checkedCC
	//int vertexNum;//the number of vertexs
	int lastTime = oriEndT - k + 1;//last start time 
	int edgeEndT;//record the currently considering end time of R edge sets
	int Te, oriEndTe;//the end time of interval
	int rightEndpoint;//max position of R set

	veciter(MidResult*) recordIter = EIntR->begin(), recordEnd = EIntR->end();
	veciter(bool) validRecordIter = validMidResult.begin();
	MidResult* resultPtr;
	vec(MotifPos)::iterator checkMotifIter = MIntR->begin(), checkMotifIterRecord, tempIter,
		checkMotifIterEnd = MIntR->end();
	int resultLastPos, resultFirstPos;
	int deltaT = endT - oriEndT;
	int transformPos = 0;
	bool* expandMask = DBG_NEW bool[currNTimestamp];
	int* subCCId = DBG_NEW int[nEdge];
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue;
	int id, Tsk = startT + k - 1;
	vec(int) edgeList;
	vec(int)* motifEdges = &edgeList;
	//traverse all possible intervals' start time
	TGraph::posUsedForEdgeFilter = startT * nEdge;
	for (int Ts = startT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge, transformPos += deltaT, Tsk++) {//O(T)
		while (recordIter != recordEnd) {//O(|L|*(R(1,T)¡È...R(T-k+1,T))) in all rows
			resultPtr = *recordIter;
			if (*validRecordIter) {
				if (resultPtr->row > Ts) break;

				if (!edgesFetched[resultPtr->edgeId]) {
					edgesFetched[resultPtr->edgeId] = true;
					edgesInEIntR->emplace_back(resultPtr->edgeId);
				}
				if (maxIntv[resultPtr->edgeId].first == -1) { // unused
					maxIntv[resultPtr->edgeId] = resultPtr->eMaxIntvl;
					resultPtr->preMaxIntv.copyTo(preMaxIntv[resultPtr->edgeId]);
				}
				int mainLabelPos = resultPtr->edgeId * numOfLabel + getEdgeLabel(resultPtr->edgeId, Ts - startT);
				vioT[mainLabelPos].copyTo(resultPtr->vioT);

				scanT[mainLabelPos] = resultPtr->scanT;
			}
			++recordIter;
			++validRecordIter;
		}
		BEGIN_TIMER(b)
		selectedNum = 0;
		rightEndpoint = 0;
		oriEndTe = oriEndT + 1;
		Te = Ts + k - 1;
		this->edgeFilterFRTMMidRForDYN(Ts, Te, oriEndTe,
			edgeSetsR, selectedNum, rightEndpoint,
			k, fixLabel, isEdgeTypeFixed);
			//initalize for every row
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;


		disjointSet->reset();
		//vertex2Pos.clear();
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		//vertexNum = 0;
		edgeEndT = endT;

		//Test S(m,T)
		smt = edgeSetsR[endT - oriEndTe].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		smk = edgeSetsR[0].size();
		maxSmk = max(maxSmk, smk);
		allSmk += smk;

		//R2L
		for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {

			if (tempT > rightEndpoint) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			BEGIN_TIMER(c)
#pragma region maxCheck

				if (iterStart != iterEnd) {//new edges fetched
					/*fetch new edges from one R edge set and insert into the disjoint set (maintain connectivity)*/
					maintainConnectivity(iterStart, iterEnd,
						/*vertex2Pos,*/ disjointSet, /*vertexNum,*/
						combineMotifPos);
					/*combine components which are connected after adding new edges from one R edge set,
						combine components before adding new edges for less computation
						(combine components in generateMaxTM case 3)*/
					combineComponentsFRTM(setCC,
						/*vertex2Pos,*/ disjointSet, root2Comp, checkedCC, hasChecked,
						/*setCCSize, realMotifNum,*/
						combineMotifPos);
				}

			/*add edge into generated motif or generate new motif and check left expandable
				(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)
			some ccs need to be saved though not new edges fetched*/
			updateNewEdgeInfoFRTM(iterStart, iterEnd,
				setCC, /*vertex2Pos,*/ disjointSet, tempDisjointSet, root2Comp, /*&subCCMap, */subCCId, root2Id, checkedCC, hasChecked, /*expandMask,
				 EdgeMaskLen,setCCSize,*/tempRecordFromQueue,
				/*realMotifNum,*/ saveMotifPos, Ts, edgeEndT, k, ComponentsType::CCFRTM);//O(¦¤Es)
			Test::gm += END_TIMER(c);

#pragma endregion

			BEGIN_TIMER(d)
				expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
					edgeEndT, expandMask, motifNumber, &TGraph::checkExpandableFRTMMidR<CComponentsFRTM>, emptyCC);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();

		}

		resultLastPos = resultPos(Ts, oriEndT, startT, oriEndT, k);
		resultFirstPos = resultPos(Ts, Tsk, startT, oriEndT, k);
		checkMotifIterRecord = checkMotifIter;
		if (checkMotifIter != checkMotifIterEnd) {
			int resultP = checkMotifIter->first, resultIndex = checkMotifIter->second;
			while (resultP <= resultLastPos) {
				int motifEndT = Tsk + resultP - resultFirstPos;
				int scanMotifEndT = max(oriEndT, motifEndT);
				resultP = checkMotifIter->first + transformPos;
				newResult[resultP][resultIndex]->getMotifEdge(newResult, motifEdges);
				auto motifEdgesEnd = motifEdges->end();
				auto motifEdgesIter = motifEdges->begin();
				bool flag = false;
				for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
					id = *motifEdgesIter;
					int inMidP = id * numOfLabel + getEdgeLabel(id, motifEndT);
					if (posInEIntR[inMidP] != EMaxIntvlChange::CHANGED) break;
				}

				if (motifEdgesIter == motifEdgesEnd) {//may change

					int maxLeft = 0, minRight = endT;
					motifEdgesIter = motifEdges->begin();
					for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
						id = *motifEdgesIter;
						int left = 0x7fffffff, right = 0;
						if (isSameLabel(id, maxIntv[id].second - startT, Ts - startT)) {
							left = maxIntv[id].first;
							right = maxIntv[id].second;
						}
						if (!preMaxIntv[id].empty()) {
							auto intvf = preMaxIntv[id].front;
							while (intvf != preMaxIntv[id].rear) {
								auto& intv = preMaxIntv[id].q[intvf];
								if (isSameLabel(id, intv.second - startT, Ts - startT)) {//need same label
									left = min(left, intv.first);
									right = max(right, intv.second);
								}
								intvf++;
							}
						}
						maxLeft = max(maxLeft, left);
						minRight = min(minRight, right);
					}

					bool isSave = true;
					if (Ts - maxLeft + 1 >= minRight - scanMotifEndT) {
						int j = minRight;
						for (; j > scanMotifEndT; --j) {
							CLEARALL(expandMask, true, Ts - maxLeft + 1, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(*motifEdges, maxLeft - startT, Ts - startT, j - startT, Ts - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}
					else {
						int j = maxLeft;
						for (; j <= Ts; ++j) {
							CLEARALL(expandMask, true, minRight - scanMotifEndT, bool);
							if (bothFitDefAndSameLabelChangeEndTimePos(*motifEdges, j - startT, scanMotifEndT + 1 - startT, minRight - startT, Ts - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}

					if (!isSave) {//expandable
						newResult[resultP][resultIndex] = nullptr;
					}
				}

				checkMotifIter++;
				if (checkMotifIter == checkMotifIterEnd) break;

				resultP = checkMotifIter->first;
				resultIndex = checkMotifIter->second;
			}

			if (checkMotifIter != checkMotifIterRecord) {//not empty
				tempIter = checkMotifIter - 1;
				while (tempIter != checkMotifIterRecord) {
					resultP = tempIter->first + transformPos;
					resultIndex = tempIter->second;
					if (newResult[resultP][resultIndex] == nullptr) {
						newResult[resultP][resultIndex] = newResult[resultP][newResult[resultP].size() - 1];
						newResult[resultP].pop_back();
						motifNumber--;
					}
					tempIter--;
				}
				//tempIter == checkMotifIterRecord
				resultP = tempIter->first + transformPos;
				resultIndex = tempIter->second;
				if (newResult[resultP][resultIndex] == nullptr) {
					newResult[resultP][resultIndex] = newResult[resultP][newResult[resultP].size() - 1];
					newResult[resultP].pop_back();
					motifNumber--;
				}
			}
		}

		//testing
		Test::updateMemoryUse();

		//release
		releaseCCs(setCC);

		if (END_TIMER(Test::startTime) > 5400000) {
			delete disjointSet;
			delete tempDisjointSet;
			delete[] expandMask;
			delete[] subCCId;
			delete[] edgeSetsR;
			delete checkedCC;
			return;
		}
	}
	delete disjointSet;
	delete tempDisjointSet;
	delete[] expandMask;
	delete[] subCCId;
	delete[] edgeSetsR;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime + 1) << endl;
	cout << "maxSmk: " << maxSmk << " averSmk: " << allSmk / (lastTime + 1) << endl;
}

void TGraph::findRTMotifsDynamicOpt1(int k,
	vec(TMotifII*)*& newResult, int oriEndT,
	long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed, ComponentsType cctype) {
#pragma region intialization
	//i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	CCTYPEOPT1* emptyCC = nullptr;
	vec(CComponentsII*) setCC;//temporary component list CC
	int realSize = endT - oriEndT + 1;
	vec(int)/*SAVESIMINFO_Vec*/* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	//unordered_set<int> edgesAdd(nEdge<<1);
	int* edgesAdd = DBG_NEW int[nEdge];
	CLEARALL(edgesAdd, 0, nEdge, int);
	DisjointSet* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	DisjointSetWithTime* tempDisjointSet = DBG_NEW DisjointSetWithTime(nNode);

	i2iHMap root2Id(nNode);
	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	int selectedNum;//the number of edges in edgeSetsR

	vec(int) saveMotifPos;//record the position of components which need to be saved
	/*record the position of components which need to combine with other components*/
	vec(int) combineMotifPos;
	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
	long long allSmk = 0, maxSmk = 0, smk;//the sum and the max number of S[m,m+k-1]
	LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map ccpos in setCC to previous node in checkedCC
	int lastTime = oriEndT - k + 1;//last start time 
	int edgeEndT;//record the currently considering end time of R edge sets
	int Te, oriEndTe;//the end time of interval
	int rightEndpoint;//max position of R set

	veciter(MidResult*) recordIter = EIntR->begin(), recordEnd = EIntR->end();
	veciter(bool) validRecordIter = validMidResult.begin();
	MidResult* resultPtr;
	vec(MotifPos)::iterator checkMotifIter = MIntR->begin(), checkMotifIterRecord, tempIter,
		checkMotifIterEnd = MIntR->end();
	int resultLastPos, resultFirstPos;
	int deltaT = endT - oriEndT;
	int transformPos = 0;
	bool* expandMask = DBG_NEW bool[currNTimestamp];
	int* subCCId = DBG_NEW int[nEdge]; 
	int* tempHaveNewEdgeCC = DBG_NEW int[nNode];
	int* haveNewEdgeCC = DBG_NEW int[nEdge];//have edge in R+
	CLEARALL(haveNewEdgeCC, -1, nEdge, int);
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue;
	vec(int) edgeList;
	vec(int)* motifEdges = &edgeList;

	int id, Tsk = startT + k - 1;
	//traverse all possible intervals' start time
	TGraph::posUsedForEdgeFilter = 0;
	for (int Ts = startT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge, transformPos += deltaT, Tsk++) {//O(T)
		//if (Ts > changeStartPosT && intvE < oriEndT) intvE++;
		while (recordIter != recordEnd) {//O(|L|*(R(1,T)¡È...R(T-k+1,T))) in all rows
			resultPtr = *recordIter;
			if (*validRecordIter) {
				if (resultPtr->row > Ts) break;

				if (!edgesFetched[resultPtr->edgeId]) {
					edgesFetched[resultPtr->edgeId] = true;
					edgesInEIntR->emplace_back(resultPtr->edgeId);
				}
				if (maxIntv[resultPtr->edgeId].first == -1) { // unused
					maxIntv[resultPtr->edgeId] = resultPtr->eMaxIntvl;
					resultPtr->preMaxIntv.copyTo(preMaxIntv[resultPtr->edgeId]);
				}
				int mainLabelPos = resultPtr->edgeId * numOfLabel + getEdgeLabel(resultPtr->edgeId, Ts - startT);
				vioT[mainLabelPos].copyTo(resultPtr->vioT);

				scanT[mainLabelPos] = resultPtr->scanT;
			}
			++recordIter;
			++validRecordIter;
		}
		BEGIN_TIMER(b)
			selectedNum = 0;
		rightEndpoint = 0;
		oriEndTe = oriEndT + 1;
		Te = Ts + k - 1;
		edgeFilterFRTMMidRForDYN(Ts, Te, oriEndTe,
			edgeSetsR, selectedNum, rightEndpoint,
			k, fixLabel, isEdgeTypeFixed, edgesAdd);
			//initalize for every row
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;


		disjointSet->reset();
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		edgeEndT = endT;

		//Test S(m,T)
		smt = edgeSetsR[endT - oriEndTe].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		smk = edgeSetsR[0].size();
		maxSmk = max(maxSmk, smk);
		allSmk += smk;

		//R2L
		//edgesAdd.clear();
		for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {

			if (tempT > rightEndpoint) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			BEGIN_TIMER(c)
#pragma region maxCheck

				if (iterStart != iterEnd) {//new edges fetched
					maintainConnectivity(iterStart, iterEnd,
						/*vertex2Pos,*/ disjointSet, /*vertexNum,*/
						combineMotifPos);
					combineComponentsFRTMOPT1(setCC,
						/*vertex2Pos,*/ disjointSet, root2Comp, checkedCC, hasChecked,
						/*setCCSize, realMotifNum,*/
						combineMotifPos, haveNewEdgeCC);
				}

			updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,
				setCC, /*vertex2Pos,*/ disjointSet, tempDisjointSet, root2Comp, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
				 saveMotifPos, Ts, edgeEndT, k, haveNewEdgeCC,tempHaveNewEdgeCC, ComponentsTypeOPT1);//O(¦¤Es)
			Test::gm += END_TIMER(c);

#pragma endregion

			BEGIN_TIMER(d)
				expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
					edgeEndT, expandMask, motifNumber, &TGraph::checkExpandableFRTMMidR<CCTYPEOPT1>,emptyCC);
			Test::gne += END_TIMER(d);
			edgeSetsR[tempT].clear();
		}

		resultLastPos = resultPos(Ts, oriEndT, startT, oriEndT, k);
		resultFirstPos = resultPos(Ts, Tsk, startT, oriEndT, k);
		checkMotifIterRecord = checkMotifIter;
		if (checkMotifIter != checkMotifIterEnd) {
			int resultP = checkMotifIter->first, resultIndex = checkMotifIter->second;
			while (resultP <= resultLastPos) {
				int motifEndT = Tsk + resultP - resultFirstPos;
				int scanMotifEndT = max(oriEndT, motifEndT);
				resultP = checkMotifIter->first + transformPos;
				newResult[resultP][resultIndex]->getMotifEdge(newResult, motifEdges);
				auto motifEdgesEnd = motifEdges->end();
				auto motifEdgesIter = motifEdges->begin();
				bool flag = false;
				for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
					id = *motifEdgesIter;
					int inMidP = id * numOfLabel + getEdgeLabel(id, motifEndT);
					
					if (posInEIntR[inMidP] != EMaxIntvlChange::CHANGED) break;
				}

				if (motifEdgesIter == motifEdgesEnd) {//may change

					int maxLeft = 0, minRight = endT;
					motifEdgesIter = motifEdges->begin();
					for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
						id = *motifEdgesIter;
						int left = 0x7fffffff, right = 0;
						if (isSameLabel(id, maxIntv[id].second - startT, Ts - startT)) {
							left = maxIntv[id].first;
							right = maxIntv[id].second;
						}
						if (!preMaxIntv[id].empty()) {
							auto intvf = preMaxIntv[id].front;
							while (intvf != preMaxIntv[id].rear) {
								auto& intv = preMaxIntv[id].q[intvf];
								if (isSameLabel(id, intv.second - startT, Ts - startT)) {//need same label
									left = min(left, intv.first);
									right = max(right, intv.second);
								}
								intvf++;
							}
						}
						maxLeft = max(maxLeft, left);
						minRight = min(minRight, right);
					}

					bool isSave = true;
					if (Ts - maxLeft + 1 >= minRight - scanMotifEndT) {
						int j = minRight;
						for (; j > scanMotifEndT; --j) {
							CLEARALL(expandMask, true, Ts - maxLeft + 1, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(*motifEdges, maxLeft - startT, Ts - startT, j - startT, Ts - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}
					else {
						int j = maxLeft;
						for (; j <= Ts; ++j) {
							CLEARALL(expandMask, true, minRight - scanMotifEndT, bool);
							if (bothFitDefAndSameLabelChangeEndTimePos(*motifEdges, j - startT, scanMotifEndT + 1 - startT, minRight - startT, Ts - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}

					if (!isSave) {//expandable
						newResult[resultP][resultIndex] = nullptr;
					}
				}

				checkMotifIter++;
				if (checkMotifIter == checkMotifIterEnd) break;

				resultP = checkMotifIter->first;
				resultIndex = checkMotifIter->second;
			}

			if (checkMotifIter != checkMotifIterRecord) {//not empty
				tempIter = checkMotifIter - 1;
				while (tempIter != checkMotifIterRecord) {
					resultP = tempIter->first + transformPos;
					resultIndex = tempIter->second;
					if (newResult[resultP][resultIndex] == nullptr) {
						newResult[resultP][resultIndex] = newResult[resultP][newResult[resultP].size() - 1];
						newResult[resultP].pop_back();
						motifNumber--;
					}
					tempIter--;
				}
				//tempIter == checkMotifIterRecord
				resultP = tempIter->first + transformPos;
				resultIndex = tempIter->second;
				if (newResult[resultP][resultIndex] == nullptr) {
					newResult[resultP][resultIndex] = newResult[resultP][newResult[resultP].size() - 1];
					newResult[resultP].pop_back();
					motifNumber--;
				}
			}
		}

		//testing
		Test::updateMemoryUse();

		//release
		releaseCCs(setCC);
	}
	delete disjointSet;
	delete[] subCCId;
	delete tempDisjointSet;
	delete[] expandMask;
	delete[] edgeSetsR;
	delete[] tempHaveNewEdgeCC;
	delete[] edgesAdd;
	delete[] haveNewEdgeCC;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime + 1) << endl;
	cout << "maxSmk: " << maxSmk << " averSmk: " << allSmk / (lastTime + 1) << endl;
}

void TGraph::findRTMotifsDynamicPlus(int k,
	vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed, ComponentsType cctype) {
	//for no noise field
	int newk = (int)(ceil(1 / Setting::delta) + 0.5);
	int maxLengthForNoNoise = newk - k;
	if (maxLengthForNoNoise >= 1) {
#pragma region intialization
		/*map the root in disjoint set to the position of its corresponding
		connected component in the component list*/
		i2iHMap root2Comp(nNode);

		CCTYPEOPT1* emptyCC = nullptr;
		vec(CComponentsII*) setCC;//temporary component list CC
		int realSize = endT - oriEndT + 1;
		vec(int)/*SAVESIMINFO_Vec*/* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
		//unordered_set<int> edgesAdd(nEdge<<1);
		int* edgesAdd = DBG_NEW int[nEdge];
		CLEARALL(edgesAdd, 0, nEdge, int);
		DisjointSet* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
		DisjointSetWithTime* tempDisjointSet = DBG_NEW DisjointSetWithTime(nNode);

		i2iHMap /*subCCMap(nNode),*/ root2Id(nNode);
		veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
		int selectedNum;//the number of edges in edgeSetsR

		vec(int) saveMotifPos;//record the position of components which need to be saved
		/*record the position of components which need to combine with other components*/
		vec(int) combineMotifPos;
		long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
		long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
		long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
		long long allSmk = 0, maxSmk = 0, smk;//the sum and the max number of S[m,m+k-1]
		LinkedList<int>* checkedCC = DBG_NEW LinkedList<int>();
		unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map ccpos in setCC to previous node in checkedCC
		//int vertexNum = 0;//the number of vertexs
		//int realMotifNum;//the number of generated connected component
		int lastTime = oriEndT - k + 1;//last start time 
		int edgeEndT;//record the currently considering end time of R edge sets
		int Te, TeNoNoise, oriEndTe;//the end time of interval
		int rightEndpoint;//max position of R set

		veciter(MidResult*) recordIter = EIntR->begin(), recordEnd = EIntR->end();
		veciter(bool) validRecordIter = validMidResult.begin();
		MidResult* resultPtr;
		vec(MotifPos)::iterator checkMotifIter = MIntR->begin(), checkMotifIterRecord, tempIter,
			checkMotifIterEnd = MIntR->end();
		int resultLastPos, resultFirstPos;
		int deltaT = endT - oriEndT;
		int transformPos = 0;
		bool* expandMask = DBG_NEW bool[currNTimestamp];
		int* subCCId = DBG_NEW int[nEdge];
		int* tempHaveNewEdgeCC = DBG_NEW int[nNode];
		int* haveNewEdgeCC = DBG_NEW int[nEdge];//have edge in R+
		CLEARALL(haveNewEdgeCC, -1, nEdge, int);
#pragma endregion
		vector<pair<int, int>> tempRecordFromQueue;

		vec(int)* edgeSetsRNoNoise = DBG_NEW vec(int)[maxLengthForNoNoise];// R edge sets
		i2iHMap root2CompNoNoise(nNode);
		vec(CComponentsShortIntv*) setCCShortIntv;//temporary component list CC
		DisjointSet* disjointSetNoNoise = DBG_NEW DisjointSet(nNode);//disjoint set
		veciter(int) iterEndNoNoise, iterStartNoNoise;//iterator for edgeSetsRNoNoise
		veciter(CComponentsShortIntv*) listEndNoNoise;//iterator for setCCShortIntv
		vec(int) saveMotifPosNoNoise;//record the position of components which need to be saved
		vec(int) combineMotifPosNoNoise;
		int selectedNumNoNoise;
		int lastTimeNoNoise = endT - newk + 1;
		bool* hasE = DBG_NEW bool[nEdge];
		//iSet hasE (nEdge);
		int* newE = DBG_NEW int[nEdge];
		int newENum;
		int* addE = DBG_NEW int[nEdge];
		int addENum;
		bool* saveR = DBG_NEW bool[nNode];
		int id, Tsk = startT + k - 1;
		vec(int) edgeList;
		vec(int)* motifEdges = &edgeList;

		//traverse all possible intervals' start time
		TGraph::posUsedForEdgeFilter = 0;
		TGraph::posUsedForEdgeFilterShortIntv = (startT + k - 1)*nEdge;
		for (int Ts = startT; Ts <= lastTime; Ts++, TGraph::posUsedForEdgeFilter += nEdge, TGraph::posUsedForEdgeFilterShortIntv += nEdge, transformPos += deltaT, Tsk++) {//O(T)
			while (recordIter != recordEnd) {//O(|L|*(R(1,T)¡È...R(T-k+1,T))) in all rows
				resultPtr = *recordIter;
				
				if (*validRecordIter) {
					if (resultPtr->row > Ts) break;
					if (!edgesFetched[resultPtr->edgeId]) {
						edgesFetched[resultPtr->edgeId] = true;
						edgesInEIntR->emplace_back(resultPtr->edgeId);
					}
					if (maxIntv[resultPtr->edgeId].first == -1) {
						maxIntv[resultPtr->edgeId] = resultPtr->eMaxIntvl;
						resultPtr->preMaxIntv.copyTo(preMaxIntv[resultPtr->edgeId]);
					}
					int mainLabelPos = resultPtr->edgeId * numOfLabel + getEdgeLabel(resultPtr->edgeId, Ts - startT);
					vioT[mainLabelPos].copyTo(resultPtr->vioT);
					scanT[mainLabelPos] = resultPtr->scanT;
				}
				++recordIter;
				++validRecordIter;
			}
			BEGIN_TIMER(b)
			selectedNum = 0;
			selectedNumNoNoise = 0;
			rightEndpoint = -1;
			oriEndTe = oriEndT + 1;
			TeNoNoise = Ts + k - 1;
			Te = Ts + newk - 1;
			
			if (Te > oriEndTe) {//deal with no distortion in special case
				if (Ts <= lastTimeNoNoise + 1) {
					CLEARALL(hasE, 0, nEdge, bool);
					CLEARALL(saveR, 0, nNode, bool);
					newENum = 0;
				}
				if (Ts <= lastTimeNoNoise) {
					if (Ts != startT) {
						root2Comp.clear();
						disjointSet->reset();
						hasChecked.clear();
						checkedCC->release();

						disjointSetNoNoise->reset();
						root2CompNoNoise.clear();
					}
					edgeFilterPlusMidRForDYN(Ts, Te, TeNoNoise, Te - oriEndTe - 1, oriEndTe, hasE, newE, newENum, edgesAdd,
						edgeSetsR, selectedNum, edgeSetsRNoNoise, selectedNumNoNoise, rightEndpoint, k, fixLabel, isEdgeTypeFixed);
				}
				else {//only no distortion
					if (Ts != startT) {
						disjointSetNoNoise->reset();
						root2CompNoNoise.clear();
					}
					edgeFilterShortIntvMidRForDYN(Ts, TeNoNoise, Te - oriEndTe - 1, hasE, newE, newENum, oriEndTe,
						lastTimeNoNoise, edgeSetsRNoNoise, selectedNumNoNoise, fixLabel, isEdgeTypeFixed);//not noise
				}
			}
			else {//same as FRTMOpt1
				if (Ts != startT) {
					root2Comp.clear();
					disjointSet->reset();
					hasChecked.clear();
					checkedCC->release();
				}
				edgeFilterFRTMMidRForDYN(Ts, Te, oriEndTe,
					edgeSetsR, selectedNum, rightEndpoint, k, fixLabel, isEdgeTypeFixed, edgesAdd);
			}
			//initalize for every row
			Test::compr += END_TIMER(b);
			maxSelect = max(maxSelect, selectedNum + selectedNumNoNoise);
			selectedSum += selectedNum + selectedNumNoNoise;

			
			
			edgeEndT = endT;

			//Test S(m,T)
			if (Te > oriEndTe) {//mix method
				if (Ts <= lastTimeNoNoise) {
					smt = edgeSetsR[endT - oriEndTe].size();
					maxSmt = max(maxSmt, smt);
					allSmt += smt;
				}
				else {
					smt = edgeSetsRNoNoise[endT - oriEndTe].size();
					maxSmt = max(maxSmt, smt);
					allSmt += smt;
				}
				smk = edgeSetsRNoNoise[0].size();
				maxSmk = max(maxSmk, smk);
				allSmk += smk;
			}
			else {
				smt = edgeSetsR[endT - oriEndTe].size();
				maxSmt = max(maxSmt, smt);
				allSmt += smt;
				smk = edgeSetsR[0].size();
				maxSmk = max(maxSmk, smk);
				allSmk += smk;
			}

			//R2L
			if (Te <= oriEndTe) {
				//edgesAdd.clear();
				for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {
					if (tempT > rightEndpoint) continue;
					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();
					BEGIN_TIMER(c)
#pragma region maxCheck
					if (iterStart != iterEnd) {//new edges fetched
						maintainConnectivity(iterStart, iterEnd,
							/*vertex2Pos,*/ disjointSet, /*vertexNum,*/
							combineMotifPos);
						combineComponentsFRTMOPT1(setCC,
							/*vertex2Pos,*/ disjointSet, root2Comp, checkedCC, hasChecked,
							combineMotifPos, haveNewEdgeCC);
					}
					updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,/**/
						setCC, /*vertex2Pos,*/ disjointSet, tempDisjointSet, root2Comp, /*&subCCMap,*/ subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
						saveMotifPos, Ts, edgeEndT, k, haveNewEdgeCC, tempHaveNewEdgeCC, ComponentsTypeOPT1);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion
					BEGIN_TIMER(d)
						expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
							edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableFRTMMidR<CCTYPEOPT1>,emptyCC);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();
				}
			}
			else {
			//edgesAdd.clear();
			for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {
				if (edgeEndT >= Te) {//large intervals
					if (tempT > rightEndpoint) continue;
					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();
					BEGIN_TIMER(c)
#pragma region maxCheck
						if (iterStart != iterEnd) {//new edges fetched
							maintainConnectivity(iterStart, iterEnd,
								/*vertex2Pos,*/ disjointSet,/* vertexNum,*/
								combineMotifPos);
							combineComponentsFRTMOPT1(setCC,
								/*vertex2Pos,*/ disjointSet, root2Comp, checkedCC, hasChecked,
								combineMotifPos, haveNewEdgeCC);
						}
					updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,/**/
						setCC, /*vertex2Pos,*/ disjointSet, tempDisjointSet, root2Comp, /*&subCCMap,*/ subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
						saveMotifPos, Ts, edgeEndT, k, haveNewEdgeCC, tempHaveNewEdgeCC, ComponentsTypeOPT1);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion
					BEGIN_TIMER(d)
						expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
							edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableFRTMMidR<CCTYPEOPT1>,emptyCC);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();
				}
				else {//short intervals

					if (Ts <= lastTimeNoNoise && edgeEndT == Te - 1) {//record root of disjointSet to be saved
						recordSavedRoot(setCC, hasChecked, newE, newENum, saveR, disjointSet/*, vertex2Pos*/);
					}
					iterStartNoNoise = edgeSetsRNoNoise[tempT].begin();
					iterEndNoNoise = edgeSetsRNoNoise[tempT].end();
					if (iterStartNoNoise == iterEndNoNoise) continue;
					BEGIN_TIMER(c)
#pragma region generateMaxTM
						addENum = 0;
					maintainConnectivityPlus(iterStartNoNoise, iterEndNoNoise, addE, addENum,
						/*vertex2Pos,*/ disjointSetNoNoise, /*vertexNum,*/ hasE, saveR, disjointSet, combineMotifPosNoNoise);
					combineComponentsPlus(setCCShortIntv,
						disjointSetNoNoise, root2CompNoNoise, combineMotifPosNoNoise);
					updateNewEdgeInfoPlus(iterStartNoNoise, addE, addENum,
						setCCShortIntv, /*vertex2Pos,*/ disjointSetNoNoise,
						root2CompNoNoise, hasE, saveMotifPosNoNoise, Ts, edgeEndT, ComponentsTypeOPT1);
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						expCheckFRTMPlusMidR(saveMotifPosNoNoise, setCCShortIntv,
							newResult, k, Ts, edgeEndT, expandMask, /*expCheck,*/ motifNumber, disjointSet,
							hasChecked, root2Comp, /*vertex2Pos,*/ setCC, emptyCC);
					Test::gne += END_TIMER(d);
					edgeSetsRNoNoise[tempT].clear();
				}
			}
			}

			resultLastPos = resultPos(Ts, oriEndT, startT, oriEndT, k);
			resultFirstPos = resultPos(Ts, Tsk, startT, oriEndT, k);
			checkMotifIterRecord = checkMotifIter;

			if (checkMotifIter != checkMotifIterEnd) {
				int resultP = checkMotifIter->first, resultIndex = checkMotifIter->second;
				while (resultP <= resultLastPos) {
					int motifEndT = Tsk + resultP - resultFirstPos;
					int scanMotifEndT = max(oriEndT, motifEndT);
					resultP = checkMotifIter->first + transformPos;
					newResult[resultP][resultIndex]->getMotifEdge(newResult, motifEdges);
					auto motifEdgesEnd = motifEdges->end();
					auto motifEdgesIter = motifEdges->begin();
					bool flag = false;
					for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
						id = *motifEdgesIter;
						int inMidP = id * numOfLabel + getEdgeLabel(id, motifEndT);
						if (posInEIntR[inMidP] != EMaxIntvlChange::CHANGED) break;
					}
					
					if (motifEdgesIter == motifEdgesEnd) {//may change

						int maxLeft = 0, minRight = endT;
						motifEdgesIter = motifEdges->begin();
						for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
							id = *motifEdgesIter;
							int left = 0x7fffffff, right = 0;
							if (maxIntvShortIntv[id].second >= motifEndT) {
								left = maxIntvShortIntv[id].first;
								right = maxIntvShortIntv[id].second;
							}
							if (maxIntv[id].second >= motifEndT && isSameLabel(id, maxIntv[id].second - startT, Ts - startT)) {//have noise
								left = min(left, maxIntv[id].first);
								right = max(right, maxIntv[id].second);
							}
							if (!preMaxIntv[id].empty()) {
								auto intvf = preMaxIntv[id].front;
								while (intvf != preMaxIntv[id].rear) {
									auto& intv = preMaxIntv[id].q[intvf];
									if (intv.second >= motifEndT && isSameLabel(id, intv.second - startT, Ts - startT)) {//need same label
										left = min(left, intv.first);
										right = max(right, intv.second);
									}
									intvf++;
								}
							}
							maxLeft = max(maxLeft, left);
							minRight = min(minRight, right);
						}

						bool isSave = true;
						if (Ts - maxLeft + 1 >= minRight - scanMotifEndT) {
							int j = minRight;
							for (; j > scanMotifEndT; --j) {
								CLEARALL(expandMask, true, Ts - maxLeft + 1, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(*motifEdges, maxLeft - startT, Ts - startT, j - startT, Ts - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						else {
							int j = maxLeft;
							for (; j <= Ts; ++j) {
								CLEARALL(expandMask, true, minRight - scanMotifEndT, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(*motifEdges, j - startT, scanMotifEndT + 1 - startT, minRight - startT, Ts - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}

						if (!isSave) {//expandable
							newResult[resultP][resultIndex] = nullptr;
						}
					}

					checkMotifIter++;
					if (checkMotifIter == checkMotifIterEnd) break;

					resultP = checkMotifIter->first;
					resultIndex = checkMotifIter->second;
				}

				if (checkMotifIter != checkMotifIterRecord) {//not empty
					tempIter = checkMotifIter - 1;
					while (tempIter != checkMotifIterRecord) {
						resultP = tempIter->first + transformPos;
						resultIndex = tempIter->second;
						if (newResult[resultP][resultIndex] == nullptr) {
							newResult[resultP][resultIndex] = newResult[resultP][newResult[resultP].size() - 1];
							newResult[resultP].pop_back();
							motifNumber--;
						}
						tempIter--;
					}
					resultP = tempIter->first + transformPos;
					resultIndex = tempIter->second;
					if (newResult[resultP][resultIndex] == nullptr) {
						newResult[resultP][resultIndex] = newResult[resultP][newResult[resultP].size() - 1];
						newResult[resultP].pop_back();
						motifNumber--;
					}
				}
			}

			//testing
			Test::updateMemoryUse();

			//release
			releaseCCs(setCC);
			releaseCCs(setCCShortIntv);
		}
		delete[] hasE;
		delete[] newE;
		delete[] addE;
		delete[] saveR;
		delete[] edgesAdd;
		delete[] haveNewEdgeCC;
		delete[] edgeSetsRNoNoise;
		delete disjointSetNoNoise;
		delete disjointSet;
		delete[] subCCId;
		delete tempDisjointSet;
		delete[] expandMask;
		delete[] edgeSetsR;
		delete[] tempHaveNewEdgeCC;
		delete checkedCC;
		cout << SELECT_EDGE << maxSelect << endl;
		cout << MEAN_SELECT_EDGE << selectedSum / (lastTime + 1) << endl;
		cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime + 1) << endl;
		cout << "maxSmk: " << maxSmk << " averSmk: " << allSmk / (lastTime + 1) << endl;
	}
	else {
		findRTMotifsDynamicOpt1(k, newResult, oriEndT, motifNumber, fixLabel, isEdgeTypeFixed, cctype);
	}
}


void TGraph::getIntersectionOfScope(vec(int)& edges, pair<int, int>& intersectIntv) {
	intersectIntv.first = 0;
	intersectIntv.second = 0x7fffffff;
	for (auto edgeId : edges) {
		intersectIntv.first = max(intersectIntv.first, scope[edgeId].first);
		intersectIntv.second = min(intersectIntv.second, scope[edgeId].second);
	}
}

