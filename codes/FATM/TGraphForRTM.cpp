#pragma once
#include"TGraph.h"

#pragma region maintain connectivity in maxCheck for FRTM and FRTMOPT1
/*fetch new edges from one R edge set and insert into the structure which maintains connectivity*/
void TGraph::maintainConnectivity(veciter(int)& infoBegin,
	veciter(int)& infoEnd, i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, int&vertexNum,
	vec(int)& combineCCPos) {
	NodePair* temp;
	int motifPos;
	i2iHMap_Iter vertexIt;
	int sId, tId, vertex;//vertexs'id in the disjoint set
	for (auto iter = infoBegin;
		iter != infoEnd; ++iter) {
		temp = &edgeList[*iter];//new edge from one R edge set
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

		motifPos = connectivity->addE(*iter, sId, tId);
		
		if (motifPos != -1)
			combineCCPos.emplace_back(motifPos);
	}
}
#pragma endregion

#pragma region maxCheck for FRTM

void TGraph::combineComponentsFRTM(vec(CComponentsII*)& setCC,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity,
	i2iHMap& root2Comp, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
	vec(int)& combineCCPos) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in setCC
	CComponentsFRTM* currentCComponents;//move currentCComponents's edges to newCComponent

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
			currentCComponents = (CComponentsFRTM*)setCC[combinePos];
			oldRoot = currentCComponents->root;
			newRoot = connectivity->findRoot(oldRoot);

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

			((CComponentsFRTM*)setCC[mapIter->second])->combineFrom(currentCComponents);

			delete currentCComponents;//currentCComponents need to be deleted
			setCC[combinePos] = nullptr;//this position will be deleted
		}
	}
#pragma endregion
	combineCCPos.clear();
}

void TGraph::updateCCFRTM(vec(CComponentsII*)& setCC, i2iHMap&  root2Comp,
	LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,/* int& realMotifNum,*/
	veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsD5Type type) {
	CComponentsFRTM* generatedCC, *newCC;
	i2iHMap_Iter mapIter;//root2Comp's iterator 
	mapIter = root2Comp.find(root);//check which cc the edge belongs to 
	int setCCSize;//size of setCC
	int edgeSize;
	int ccPos;
	ForbidPairNode* now;
	int label = getEdgeLabel(id, startTime - startT);
	scope[id] = make_pair(-1, -1);
	if (mapIter == root2Comp.end()) {//generateMaxTM case 1
		newCC = (CComponentsFRTM*)CComponentsIID5Factory::instance(root, type);
		newCC->edges.emplace_back(id);

		scope[id] = maxIntv[id];
		newCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		newCC->minEMaxIntvlEndT = scope[id].second;
		if (!preMaxIntv[id].empty()) {
			auto intvf = preMaxIntv[id].front;
			while (intvf != preMaxIntv[id].rear) {
				auto& intv = preMaxIntv[id].q[intvf];
				if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
					newCC->preEMaxIntvlChangePos.emplace(intv.first, intv.second, id);
				}
				intvf++;
			}
		}

		now = vioT[numOfLabel*id + label].first;
		while (now != nullptr) {//O(delta log delta)
			int intvStartT, intvEndT;
			tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
			if (intvStartT > endTime) break;
			else if (intvStartT >= filterTime && intvStartT <= endTime) {
				newCC->noisePosQueue.emplace(
					min(intvEndT, endTime), // noise start position
					NoisePos(
						intvStartT, // noise end position
						0 // edge position in newCC->edges
					)

				);
			}
			now = now->next;
		}

		newCC->tabuTChangePos = -1;

		setCCSize = (int)setCC.size();
		setCC.emplace_back(newCC);
		root2Comp[root] = setCCSize;//update root2Comp

		//need checked
		hasChecked[setCCSize] = checkedCC->tail;
		checkedCC->addItemAtLast(setCCSize);
	}
	else {//generateMaxTM case 2
		ccPos = mapIter->second;
		generatedCC = (CComponentsFRTM*)setCC[ccPos];
		edgeSize = (int)generatedCC->edges.size();
		generatedCC->edges.emplace_back(id/*, label*/);

		scope[id] = maxIntv[id];
		generatedCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		if (generatedCC->minEMaxIntvlEndT > scope[id].second)
			generatedCC->minEMaxIntvlEndT = scope[id].second;
		if (!preMaxIntv[id].empty()) {
			auto intvf = preMaxIntv[id].front;
			while (intvf != preMaxIntv[id].rear) {
				auto& intv = preMaxIntv[id].q[intvf];
				if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
					generatedCC->preEMaxIntvlChangePos.emplace(intv.first, intv.second, id);
				}
				intvf++;
			}
		}

		now = vioT[numOfLabel*id + label].first;
		while (now != nullptr) {//O(delta log delta)
			int intvStartT, intvEndT;
			tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
			if (intvStartT > endTime) break;
			else if (intvStartT >= filterTime && intvStartT <= endTime) {
				generatedCC->noisePosQueue.emplace(
					min(intvEndT, endTime), // noise start position
					NoisePos(
						intvStartT, // noise end position
						edgeSize // edge position in newCC->edges
					)
				);
			}
			now = now->next;
		}
		generatedCC->tabuTChangePos = -1;

		//need checked
		if (hasChecked.find(ccPos) == hasChecked.end()) {
			hasChecked[ccPos] = checkedCC->tail;
			checkedCC->addItemAtLast(ccPos);
		}
	}

}

void TGraph::maintainTempCCForUF(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
	i2iHMap& vertex2Pos, DynamicConnectivity* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTM* generatedCC) {
	//BEGIN_NSTIMER(p);

	DisjointSet* ds = (DisjointSet*)connectivity;

	subCCMap->clear();
	root2Id.clear();
	//get the number of new connected components
	auto tempEdges = &generatedCC->edges;
	int size = (int)tempEdges->size();
	int ufsetSize = min((int)size << 1, nNode);
	ds->reset(ufsetSize);
	veciter(int) motifEdgesEnd = tempEdges->end();
	auto subCCIter = &subCCId[0];
	int vertexNum = 0, sId, tId;
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter) {//O(Em)
		if (*subCCIter != -1) {
			auto tempEdge = &edgeList[*iter];

			int vertex = tempEdge->first;
			auto vertexIt = subCCMap->find(vertex);//O(1)
			if (vertexIt == subCCMap->end()) {
				sId = vertexNum;
				(*subCCMap)[vertex] = vertexNum++;
			}
			else sId = vertexIt->second;
			vertex = tempEdge->second;
			vertexIt = subCCMap->find(vertex);//O(1)
			if (vertexIt == subCCMap->end()) {
				tId = vertexNum;
				(*subCCMap)[vertex] = vertexNum++;
			}
			else tId = vertexIt->second;
			/*union-find operation means that sId and tId are connected*/
			ds->addE(*iter, sId, tId);
		}
	}
	//Test::counter+= END_NSTIMER(p);
	//BEGIN_NSTIMER(q);
	//get subCCNum
	int subCCNum = 0;
	subCCIter = &subCCId[0];
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter) {//O(Em)
		if (*subCCIter != -1) {
			/*the root of the edge's vertex in the disjoint set*/
			int root = ds->findRoot((*subCCMap)[edgeList[*iter].first]);
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
	//Test::counter2 += END_NSTIMER(q);

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

	subCCIter = &subCCId[0];
	int edgePos = 0, ccid;
	for (auto iter = generatedCC->edges.begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
		ccid = *subCCIter;
		if (ccid != -1) {
			generatedCC->subCCs[ccid].emplace_back(edgePos);
			generatedCC->subCCsaved[ccid] = edgePos >= generatedCC->newInsert ? true : false;//right non-expandable
		}
	}
}

void TGraph::updateNewEdgeInfoFRTM(
	veciter(int)& infoBegin, veciter(int)& infoEnd,
	vec(CComponentsII*)& setCC,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, DisjointSet* tempDisjointSet,
	i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue,
	vec(int)& saveCCPos, int startTime, int endTime, int k, MaintainCC mainStrategy, UpdateComponent updateStrategy, ComponentsD5Type type) {
	if (checkedCC->first == nullptr && infoBegin == infoEnd) return;

	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int filterTime = startTime + k - 1;

	//auto beginTest = std::chrono::steady_clock::now();
	for (auto infoIter = infoBegin;
		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets

		id = *infoIter;//edge's id

		/*the root of the edge's vertex in the disjoint set*/
		root = connectivity->findRoot(vertex2Pos[edgeList[id].first]);
		(this->*updateStrategy)(setCC, root2Comp, checkedCC, hasChecked, /*realMotifNum,*/ infoIter,
			id, root, filterTime, startTime, endTime, k, type);
	}

	LinkedNode<int>* tempCCIter = checkedCC->first;
	if (tempCCIter == nullptr) return;
	int size, queueSize;
	int ccPos;//cc's position in setCC
	/*the cc which has already generated/will generate*/
	CComponentsFRTM* generatedCC;
	veciter(int) motifEdgesEnd;
	veciter(int) subCCIter;

	i2iHMap_Iter vertexIt, ccIt;
	int removeEdge;
	vec(int)* tempEdges;
	while (tempCCIter != nullptr) {
		ccPos = tempCCIter->item;
		generatedCC = (CComponentsFRTM*)setCC[ccPos];
		tempEdges = &generatedCC->edges;

		size = (int)tempEdges->size();
		if (generatedCC->newInsert == size) {//not updated edges
			if (generatedCC->tabuTChangePos != -1 && generatedCC->tabuTChangePos < endTime) {
				tempCCIter = tempCCIter->next;
				continue;
			}
		}

		removeEdge = 0;
		queueSize = (int)generatedCC->noisePosQueue.size();
		if (queueSize != 0) {
			CLEARALL(subCCId, 0, size, int); //generatedCC->subCCId.assign(size, 0);
		}

		generatedCC->tabuTChangePos = -1;
		tempRecordFromQueue.clear();
		while (queueSize != 0) {//worst time O(Em log (delta Em))
			auto ccLastNoisePair = generatedCC->noisePosQueue.top();

			if (ccLastNoisePair.first < endTime)
				break;
			else {
				if (ccLastNoisePair.first == endTime) {//edge has noise in this time
					//cout << startTime << " " << endTime << " " << ccPos << " " << generatedCC->edges[ccLastNoisePair.second.pos] << endl;

					removeEdge++;
					subCCId[ccLastNoisePair.second.pos] = -1;//remove this edge temporarily
					tempRecordFromQueue.emplace_back(ccLastNoisePair.second.startTime, // noise end position
						ccLastNoisePair.second.pos);
					if (generatedCC->tabuTChangePos < ccLastNoisePair.second.startTime - 1) {
						generatedCC->tabuTChangePos = ccLastNoisePair.second.startTime - 1;
					}
					generatedCC->noisePosQueue.pop();
				}
				else {
					generatedCC->noisePosQueue.pop();
					if (endTime >= ccLastNoisePair.second.startTime) {
						removeEdge++;
						subCCId[ccLastNoisePair.second.pos] = -1;//remove this edge temporarily
						tempRecordFromQueue.emplace_back(ccLastNoisePair.second.startTime, // noise end position
							ccLastNoisePair.second.pos);
						if (generatedCC->tabuTChangePos < ccLastNoisePair.second.startTime - 1) {
							generatedCC->tabuTChangePos = ccLastNoisePair.second.startTime - 1;
						}
					}
				}
			}
			queueSize = (int)generatedCC->noisePosQueue.size();
		}

		if (generatedCC->tabuTChangePos != -1) {
			if (generatedCC->noisePosQueue.size() != 0) {
				auto ccLastNoisePair = generatedCC->noisePosQueue.top();
				//cout << ccLastNoisePair.first << " " << ccLastNoisePair.second.pos << " " << endTime << endl;
				generatedCC->tabuTChangePos = max(ccLastNoisePair.first, generatedCC->tabuTChangePos);
			}

			int nextTime = endTime - 1;
			for (auto item : tempRecordFromQueue) {
				//cout << item.first << " " << item.second << endl;
				if (nextTime >= item.first) {
					generatedCC->noisePosQueue.emplace(
						nextTime, // noise next position
						NoisePos(
							item.first, // noise end position
							item.second // edge position in newCC->edges
						)
					);
				}
			}
		}

		while (generatedCC->preEMaxIntvlChangePos.size() != 0) {//worst time O(log (delta Em))
			auto preEMaxIntvlPair = generatedCC->preEMaxIntvlChangePos.top();

			if (preEMaxIntvlPair.endTime < endTime)
				break;
			else {
				/*if (generatedCC->maxEMaxIntvlStartT == EMaxIntvlStartT[preEMaxIntvlPair.id]) {
					update = true;
				}*/
				id = preEMaxIntvlPair.id;
				scope[id].first = preEMaxIntvlPair.startTime;
				if (scope[id].second < preEMaxIntvlPair.endTime) {
					scope[id].second = preEMaxIntvlPair.endTime;
				}

				auto maxEMaxIntvlStartTPair = generatedCC->maxEMaxIntvlStartTQueue.top();
				if (maxEMaxIntvlStartTPair.second == id) {
					while (maxEMaxIntvlStartTPair.first != scope[id].first) {
						generatedCC->maxEMaxIntvlStartTQueue.pop();
						generatedCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
						maxEMaxIntvlStartTPair = generatedCC->maxEMaxIntvlStartTQueue.top();
						id = maxEMaxIntvlStartTPair.second;
					}
				}
				generatedCC->preEMaxIntvlChangePos.pop();
			}
		}
		if (removeEdge > 0) {
			if (removeEdge == size) {
				tempCCIter = tempCCIter->next;
				generatedCC->newInsert = size;
				continue;
			}

			if (tempDisjointSet)
				(this->*mainStrategy)(tempCCIter, tempRecordFromQueue, vertex2Pos, tempDisjointSet, root2Comp, subCCMap, subCCId, root2Id, generatedCC);
			else
				(this->*mainStrategy)(tempCCIter, tempRecordFromQueue, vertex2Pos, connectivity, root2Comp, subCCMap, subCCId, root2Id, generatedCC);

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
	}
}

void TGraph::checkExpandableFRTM(int savePos, CComponentsII* temp,
	int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber) {
	CComponentsFRTM*tempCC = (CComponentsFRTM*)temp;
	TMotifII* motif;
	//newly generated connected components
	bool isRightExpandable;
	pair<int, int> intersectIntv;
	//int edgeBefP = motifEndT * nEdge;

	if (tempCC->subCCNum == 0) {//not edges removed  (must be non right expandable)
		if (motifStartT - startT != 0) {
			intersectIntv.first = tempCC->maxEMaxIntvlStartTQueue.top().first;
			intersectIntv.second = tempCC->minEMaxIntvlEndT;
			if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
				Test::gne121 += tempCC->edges.size();
				Test::gnefield++;

				motif = DBG_NEW TMotifII(motifStartT, motifEndT);
				motif->copyEdges(tempCC->edges);
				result[savePos].emplace_back(motif);
				motifNumber++;
			}
			else {//left expandable or both expandable
				Test::gne122 += tempCC->edges.size();
				int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first + 1);
				Test::gnefield += tempField * tempCC->edges.size();
				Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/
				
				/*if (allNTimestamp > nEdge) {
					int maxP = 0x7fffffff, minP = -1;
					for (auto e : tempCC->edges) {
						int edgeP = edgeBefP + e;
						maxP = min(maxP, edgeBef[edgeP].first);
						minP = max(minP, edgeBef[edgeP].second);
						if (minP > maxP) break;
					}
					if (minP <= maxP && minP != -1) return;
				}*/
				bool isSave = true;

				if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
					//cout << Test::counter << " ! " << Test::counter2 << " " << Test::counter3 << " " << intersectIntv.first << " " << intersectIntv.second << endl;
					int j = intersectIntv.second;
					for (; j >= motifEndT; --j) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				else {
					int j = intersectIntv.first;
					for (; j <= motifStartT - 1; ++j) {
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
						if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				if (isSave) {
					motif = DBG_NEW TMotifII(motifStartT, motifEndT);
					motif->copyEdges(tempCC->edges);

					result[savePos].emplace_back(motif);
					motifNumber++;
				}
			}
		}
		else {
			Test::gne11 += tempCC->edges.size();
			Test::gnefield++;

			motif = DBG_NEW TMotifII(motifStartT, motifEndT);
			motif->copyEdges(tempCC->edges);
			result[savePos].emplace_back(motif);
			motifNumber++;
		}
		//Test::counter +=END_NSTIMER;

		/*if (allNTimestamp > nEdge) {
			for (auto e : tempCC->edges) {
				int edgeP = edgeBefP + e;
				if (edgeBef[edgeP].second + 1 == motifStartT) {
					if (motifStartT == 0) {
						edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
					}
					else edgeBef[edgeP].second = motifStartT;
				}
				else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
			}
		}*/
	}
	else {
		int subCCNum = tempCC->subCCNum;

		//check each subCC
		veciter(int) subCCEnd;
		for (int i = 0; i < subCCNum; ++i) {

			isRightExpandable = !tempCC->subCCsaved[i];//check whether exists e¡ÊR+[m,i]
			subCCEnd = tempCC->subCCs[i].end();
			if (!isRightExpandable) {//exists e¡ÊR+[m,i] (must be non right expandable)
				//BEGIN_NSTIMER; 
				if (motifStartT != 0) {

					getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
					if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
						Test::gne2121 += tempCC->subCCs[i].size();
						Test::gnefield++;

						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							motif->addEdge(tempCC->edges[*edgeIter]);
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}
					else {//left expandable or both expandable
						Test::gne2122 += tempCC->subCCs[i].size();
						int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first + 1);
						Test::gnefield += tempField * tempCC->subCCs[i].size();
						Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/
						//if (allNTimestamp > nEdge) {
						//	int maxP = 0x7fffffff, minP = -1;
						//	for (auto p : tempCC->subCCs[i]) {
						//		auto e = tempCC->edges[p];
						//		int edgeP = edgeBefP + e;
						//		maxP = min(maxP, edgeBef[edgeP].first);
						//		minP = max(minP, edgeBef[edgeP].second);
						//		if (minP > maxP) break;
						//	}
						//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
						//	if (minP <= maxP && minP != -1) return;
						//}

						bool isSave = true;
						if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
							int j = intersectIntv.second;
							for (; j >= motifEndT; --j) {
								CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						else {
							int j = intersectIntv.first; 
							for (; j <= motifStartT - 1; ++j) {
								CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}

						if (isSave) {
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							subCCEnd = tempCC->subCCs[i].end();
							for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
								motif->addEdge(tempCC->edges[*edgeIter]);
							}
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
					}
				}
				else {
					Test::gne211 += tempCC->subCCs[i].size();
					Test::gnefield++;

					motif = DBG_NEW TMotifII(motifStartT, motifEndT);
					subCCEnd = tempCC->subCCs[i].end();
					for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
						motif->addEdge(tempCC->edges[*edgeIter]);
					}
					result[savePos].emplace_back(motif);
					motifNumber++;
				}
				//Test::counter += END_NSTIMER;
			}
			else {
				//BEGIN_NSTIMER; 
				getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
				if (motifStartT == 0 || intersectIntv.first == motifStartT) {//only check right expandable

					Test::gne221 += tempCC->subCCs[i].size();
					int tempField;
					if (intersectIntv.second - motifEndT - 1 > 0) {
						tempField = (intersectIntv.second - motifEndT - 1) * (motifStartT - intersectIntv.first + 1);
						//cout << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << endl;
					}
					else tempField = 0;
					Test::gnefield += tempField * tempCC->subCCs[i].size();
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/

					if (motifStartT == 0 && intersectIntv.second == endT) continue;//right expandable
					CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
					if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
						//if(subCCNum == 1) cout << 1 << endl;

						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							motif->addEdge(tempCC->edges[*edgeIter]);
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}
					//Test::counter += END_NSTIMER;
				}
				else {//left expandable, right expandable or both expandable

					Test::gne222 += tempCC->subCCs[i].size();
					int tempField = (intersectIntv.second - motifEndT) * (motifStartT - intersectIntv.first);
					Test::gnefield += tempField * tempCC->subCCs[i].size();
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/
					//if (allNTimestamp > nEdge) {
					//	int maxP = 0x7fffffff, minP = -1;
					//	for (auto p : tempCC->subCCs[i]) {
					//		auto e = tempCC->edges[p];
					//		int edgeP = edgeBefP + e;
					//		maxP = min(maxP, edgeBef[edgeP].first);
					//		minP = max(minP, edgeBef[edgeP].second);
					//		if (minP > maxP) break;
					//	}
					//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
					//	if (minP <= maxP && minP != -1) return;
					//}

					bool isSave = true;
					if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
						int j = intersectIntv.second;
						for (; j > motifEndT; --j) {
							CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;

								break;
							}
						}
					}
					else {
						int j = intersectIntv.first;
						for (; j <= motifStartT; ++j) {
							CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
							if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}

					if (isSave) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;

						}

					}

					//if (subCCNum == 1)cout << isSave  << endl;
					if (isSave) {
						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							motif->addEdge(tempCC->edges[*edgeIter]);
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}

					//Test::counter2 += END_NSTIMER;
				}
			}

		/*	if (allNTimestamp > nEdge) {
				for (auto p : tempCC->subCCs[i]) {
					auto e = tempCC->edges[p];
					int edgeP = edgeBefP + e;
					if (edgeBef[edgeP].second + 1 == motifStartT) {
						if (motifStartT == 0) {
							edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
						}
						else edgeBef[edgeP].second = motifStartT;
					}
					else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
				}
			}*/
		}
	}
}

void TGraph::checkExpandableFRTMMidR(int savePos, CComponentsII* temp,
	int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/  vec(TMotifII*)*& result, long long& motifNumber) {
	CComponentsFRTM*tempCC = (CComponentsFRTM*)temp;
	TMotifII* motif;
	//newly generated connected components
	bool isRightExpandable;
	pair<int, int> intersectIntv;
	bool needCheck;
	int edgeId;
	auto edgeEnd = tempCC->edges.end();
	//int edgeBefP = motifEndT * nEdge;
	if (tempCC->subCCNum == 0) {//not edges removed  (must be non right expandable)
		if (motifStartT - startT != 0) {
			//getIntersectionOfScope(temp->edges, intersectIntv);
			intersectIntv.first = tempCC->maxEMaxIntvlStartTQueue.top().first;
			intersectIntv.second = tempCC->minEMaxIntvlEndT;

			//intersectIntv.second = tempCC->minEMaxIntvlEndTQueue.top().first;
			if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
				motif = DBG_NEW TMotifII(motifStartT, motifEndT);
				needCheck = true;
				for (auto edgeIter = tempCC->edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
					edgeId = *edgeIter;
					if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) { needCheck = false; break; }
				}
				if (needCheck) {
					newMIntR->emplace_back(savePos, (int)result[savePos].size());
				}
				motif->copyEdges(tempCC->edges);
				result[savePos].emplace_back(motif);
				motifNumber++;
			}
			else {//left expandable or both expandable
				bool isSave = true;

				/*if (allNTimestamp > nEdge) {
					int maxP = 0x7fffffff, minP = -1;
					for (auto e : tempCC->edges) {
						int edgeP = edgeBefP + e;
						maxP = min(maxP, edgeBef[edgeP].first);
						minP = max(minP, edgeBef[edgeP].second);
						if (minP > maxP) break;
					}
					if (minP <= maxP && minP != -1) return;
				}*/
				if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
					//cout << Test::counter << " ! " << Test::counter2 << " " << Test::counter3 << " " << intersectIntv.first << " " << intersectIntv.second << endl;
					int j = intersectIntv.second;
					for (; j >= motifEndT; --j) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				else {
					int j = intersectIntv.first;
					for (; j <= motifStartT - 1; ++j) {
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
						if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}

				if (isSave) {
					motif = DBG_NEW TMotifII(motifStartT, motifEndT);

					needCheck = true;
					for (auto edgeIter = tempCC->edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
						edgeId = *edgeIter;
						if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) { needCheck = false; break; }
					}
					if (needCheck) {
						newMIntR->emplace_back(savePos, (int)result[savePos].size());
					}
					motif->copyEdges(tempCC->edges);
					result[savePos].emplace_back(motif);
					motifNumber++;
				}
			}
		}
		else {
			motif = DBG_NEW TMotifII(motifStartT, motifEndT);

			needCheck = true;
			for (auto edgeIter = tempCC->edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
				edgeId = *edgeIter;
				if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) { needCheck = false; break; }
			}
			if (needCheck) {
				newMIntR->emplace_back(savePos, (int)result[savePos].size());
			}
			motif->copyEdges(tempCC->edges);
			result[savePos].emplace_back(motif);
			motifNumber++;
		}
		//Test::counter +=END_NSTIMER;

		/*if (allNTimestamp > nEdge) {
			for (auto e : tempCC->edges) {
				int edgeP = edgeBefP + e;
				if (edgeBef[edgeP].second + 1 == motifStartT) {
					if (motifStartT == 0) {
						edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
					}
					else edgeBef[edgeP].second = motifStartT;
				}
				else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
			}
		}*/
	}
	else {
		int subCCNum = tempCC->subCCNum;

		//check each subCC
		veciter(int) subCCEnd;
		for (int i = 0; i < subCCNum; ++i) {

			isRightExpandable = !tempCC->subCCsaved[i];//check whether exists e¡ÊR+[m,i]
			subCCEnd = tempCC->subCCs[i].end();

			if (!isRightExpandable) {//exists e¡ÊR+[m,i] (must be non right expandable)
				//BEGIN_NSTIMER; 
				if (motifStartT != 0) {
					getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
					if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();

						needCheck = true;
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							edgeId = tempCC->edges[*edgeIter];
							motif->addEdge(edgeId);
							if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
						}
						if (needCheck) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}
					else {//left expandable or both expandable
						//if (allNTimestamp > nEdge) {
						//	int maxP = 0x7fffffff, minP = -1;
						//	for (auto p : tempCC->subCCs[i]) {
						//		auto e = tempCC->edges[p];
						//		int edgeP = edgeBefP + e;
						//		maxP = min(maxP, edgeBef[edgeP].first);
						//		minP = max(minP, edgeBef[edgeP].second);
						//		if (minP > maxP) break;
						//	}
						//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
						//	if (minP <= maxP && minP != -1) return;
						//}
						bool isSave = true;
						if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
							int j = intersectIntv.second;
							for (; j >= motifEndT; --j) {
								CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						else {
							int j = intersectIntv.first;
							for (; j <= motifStartT - 1; ++j) {
								CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}

						if (isSave) {
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							subCCEnd = tempCC->subCCs[i].end();

							needCheck = true;
							for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
								edgeId = tempCC->edges[*edgeIter];
								motif->addEdge(edgeId);
								if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
							}
							if (needCheck) {
								newMIntR->emplace_back(savePos, (int)result[savePos].size());
							}
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
					}
				}
				else {
					motif = DBG_NEW TMotifII(motifStartT, motifEndT);
					subCCEnd = tempCC->subCCs[i].end();

					needCheck = true;
					for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
						edgeId = tempCC->edges[*edgeIter];
						motif->addEdge(edgeId);
						if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
					}
					if (needCheck) {
						newMIntR->emplace_back(savePos, (int)result[savePos].size());
					}
					result[savePos].emplace_back(motif);
					motifNumber++;
				}
				//Test::counter += END_NSTIMER;
			}
			else {
				//BEGIN_NSTIMER; 
				//if (tempCC->subCCs[i].size() > 1) {//only one edge must be expandable
				getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
				//cout << Test::counter << " # " << Test::counter2 << " " << Test::counter3 << " " << intersectIntv.first << " " << intersectIntv.second << endl;
				if (motifStartT == 0 || intersectIntv.first == motifStartT) {//only check right expandable
					CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
					if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();

						needCheck = true;
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							edgeId = tempCC->edges[*edgeIter];
							motif->addEdge(tempCC->edges[*edgeIter]);

							if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
						}
						if (needCheck) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}

						result[savePos].emplace_back(motif);
						motifNumber++;
					}
					//Test::counter += END_NSTIMER;
				}
				else {//left expandable, right expandable or both expandable
					//if (allNTimestamp > nEdge) {
					//	int maxP = 0x7fffffff, minP = -1;
					//	for (auto p : tempCC->subCCs[i]) {
					//		auto e = tempCC->edges[p];
					//		int edgeP = edgeBefP + e;
					//		maxP = min(maxP, edgeBef[edgeP].first);
					//		minP = max(minP, edgeBef[edgeP].second);
					//		if (minP > maxP) break;
					//	}
					//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
					//	if (minP <= maxP && minP != -1) return;
					//}
					bool isSave = true;
					if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
						int j = intersectIntv.second;
						for (; j > motifEndT; --j) {
							CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;

								break;
							}
						}
					}
					else {
						int j = intersectIntv.first;
						for (; j <= motifStartT; ++j) {
							CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
							if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}

					if (isSave) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;

						}
					}

					if (isSave) {
						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();

						needCheck = true;
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							edgeId = tempCC->edges[*edgeIter];
							motif->addEdge(tempCC->edges[*edgeIter]);

							if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
						}
						if (needCheck) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}
					//Test::counter2 += END_NSTIMER;
				}
			}

			/*if (allNTimestamp > nEdge) {
				for (auto p : tempCC->subCCs[i]) {
					auto e = tempCC->edges[p];
					int edgeP = edgeBefP + e;
					if (edgeBef[edgeP].second + 1 == motifStartT) {
						if (motifStartT == 0) {
							edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
						}
						else edgeBef[edgeP].second = motifStartT;
					}
					else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
				}
			}*/
		}
	}
}

#pragma endregion

#pragma region maxCheck for FRTMOpt1 

void TGraph::combineComponentsOpt1(vec(CComponentsII*)& setCC,
	i2iHMap& vertex2Pos, /*DisjointSet*& disjointSet,*/DynamicConnectivity*& connectivity,
	i2iHMap& root2Comp, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
	vec(int)& combineCCPos) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in setCC
	CComponentsFRTMOPT1* currentCComponents;//move currentCComponents's edges to newCComponent

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
			currentCComponents = (CComponentsFRTMOPT1*)setCC[combinePos];
			oldRoot = currentCComponents->root;
			newRoot = connectivity->findRoot(oldRoot);
			
			//root2Comp.erase(oldRoot);
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
			
			((CComponentsFRTMOPT1*)setCC[mapIter->second])->combineFrom(currentCComponents);

			delete currentCComponents;//currentCComponents need to be deleted
			setCC[combinePos] = nullptr;//this position will be deleted
			//realMotifNum--;
		}
	}
#pragma endregion
	combineCCPos.clear();
}


void TGraph::updateCCOpt1(vec(CComponentsII*)& setCC, i2iHMap&  root2Comp, unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd, /*int*& lastEdgeSetR,*/
	LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,/* int& realMotifNum,*/
	veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k) {
	CComponentsFRTMOPT1* generatedCC, *newCC;
	i2iHMap_Iter mapIter;//root2Comp's iterator 
	mapIter = root2Comp.find(root);//check which cc the edge belongs to 
	int setCCSize;//size of setCC
	int edgeSize;
	int ccPos;
	ForbidPairNode* now;
	int label = getEdgeLabel(id, startTime - startT);
	if (mapIter == root2Comp.end()) {//generateMaxTM case 1
		newCC = (CComponentsFRTMOPT1*)CComponentsIID5Factory::instance(root, ComponentsD5Type::CCPLUS);
		newCC->edges.emplace_back(id/*, label*/);

		scope[id] = maxIntv[id];
		newCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		newCC->minEMaxIntvlEndT = scope[id].second;
		if (!preMaxIntv[id].empty()) {
			auto intvf = preMaxIntv[id].front;
			while (intvf != preMaxIntv[id].rear) {
				auto& intv = preMaxIntv[id].q[intvf];
				if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
					newCC->preEMaxIntvlChangePos.emplace(intv.first, intv.second, id);
				}
				intvf++;
			}
		}

		now = vioT[numOfLabel*id + label].first;
		while (now != nullptr) {//O(delta log delta)
			int intvStartT, intvEndT;
			tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
			if (intvStartT > endTime) break;
			else if (intvStartT >= filterTime && intvStartT <= endTime) {
				newCC->noisePosQueue.emplace(
					min(intvEndT, endTime), // noise start position
					NoisePos(
						intvStartT, // noise end position
						0 // edge position in newCC->edges
					)

				);
			}
			now = now->next;
		}

		newCC->tabuTChangePos = -1;

		setCCSize = (int)setCC.size();
		setCC.emplace_back(newCC);
		root2Comp[root] = setCCSize;//update root2Comp

		if (edgesAdd.find(id) != edgesAddEnd) {
			newCC->haveNewEdges = 1;//lastEdgeSetR[id] + 1; //must > 0
			//need checked
			//cout << startTime<<" % "<<setCCSize << endl;
			hasChecked[setCCSize] = checkedCC->tail;
			checkedCC->addItemAtLast(setCCSize);
		}
		//realMotifNum++;
	}
	else {//generateMaxTM case 2
		ccPos = mapIter->second;
		generatedCC = (CComponentsFRTMOPT1*)setCC[ccPos];
		edgeSize = (int)generatedCC->edges.size();
		generatedCC->edges.emplace_back(id/*, label*/);

		scope[id] = maxIntv[id];

		generatedCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
		if (generatedCC->minEMaxIntvlEndT > scope[id].second)
			generatedCC->minEMaxIntvlEndT = scope[id].second;
		if (!preMaxIntv[id].empty()) {
			auto intvf = preMaxIntv[id].front;
			while (intvf != preMaxIntv[id].rear) {
				auto& intv = preMaxIntv[id].q[intvf];
				if (getEdgeLabel(id, intv.second - startT) == label) {//need same label
					generatedCC->preEMaxIntvlChangePos.emplace(intv.first, intv.second, id);
				}
				intvf++;
			}
		}

		now = vioT[numOfLabel*id + label].first;
		while (now != nullptr) {//O(delta log delta)
			int intvStartT, intvEndT;
			tie(intvStartT, intvEndT/*, std::ignore*/) = now->item;
			if (intvStartT > endTime) break;
			else if (intvStartT >= filterTime && intvStartT <= endTime) {
				generatedCC->noisePosQueue.emplace(
					min(intvEndT, endTime), // noise start position
					NoisePos(
						intvStartT, // noise end position
						edgeSize // edge position in newCC->edges
					)
				);
			}
			now = now->next;
		}
		generatedCC->tabuTChangePos = -1;

		if (!generatedCC->haveNewEdges&&edgesAdd.find(id) != edgesAddEnd) {
			generatedCC->haveNewEdges = true;
		}
		//need checked
		if (generatedCC->haveNewEdges && hasChecked.find(ccPos) == hasChecked.end()) {
			//cout << startTime << " % " << ccPos << endl;
			hasChecked[ccPos] = checkedCC->tail;
			checkedCC->addItemAtLast(ccPos);
		}
	}
}

bool TGraph::maintainTempCCForUFOpt1(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue, unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,  i2iHMap& vertex2Pos, DynamicConnectivity* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTMOPT1* generatedCC, int endTime, int removeEdge, i2iHMap& haveNewEdgeCC, int*& remainEdges) {
	//BEGIN_NSTIMER(p);
	DisjointSet* ds = (DisjointSet*)connectivity;
	auto tempEdges = &generatedCC->edges;
	int* subCCIter;
	int size = (int)tempEdges->size();
	veciter(int) motifEdgesEnd = tempEdges->end();
	int edgePos = 0, ccid;
	int subCCNum = 0, realSubCCNum = 0;
	int remainEdgesSize = size - removeEdge;
	
	if (remainEdgesSize == 1) {//one cc
		//BEGIN_NSTIMER(a1); Test::counter2++;
		if (generatedCC->subCCs == nullptr) {
			generatedCC->subCCs = DBG_NEW vec(int)[1];
			generatedCC->subCCsaved = DBG_NEW bool[1];
		}
		else {
			generatedCC->subCCs[0].clear();
		}
		generatedCC->subCCNum = 1;
		//generatedCC->subCCHaveNewEdges[0] = false;
		subCCIter = &subCCId[0];
		bool allCCHaveNoNewEdges = false;
		for (auto iter = generatedCC->edges.begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
			ccid = *subCCIter;
			if (ccid != -1) {
				if (edgesAdd.find(*iter) != edgesAddEnd) {
					generatedCC->subCCs[0].emplace_back(edgePos);
					generatedCC->subCCsaved[0] = edgePos >= generatedCC->newInsert ? true : false;//true: right non-expandable
					//Test::counter += END_NSTIMER(a1);
					return true;
				}
				else {
					//Test::counter += END_NSTIMER(a1);
					return false;
				}
			}
		}
		//Test::counter += END_NSTIMER(a1);
		return false;
	}
	else if (remainEdgesSize == 2) {
		//BEGIN_NSTIMER(a2); Test::counter4++;
		bool allCCHaveNoNewEdges = false, isNewEdges1, isNewEdges2;
		subCCIter = &subCCId[0];
		int edge1Pos = -1, edge2Pos = -1, edge1V, edge1U, edge2V, edge2U, edgeId1, edgeId2;
		for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
			if (*subCCIter != -1) {
				if (edge1Pos == -1) {
					edgeId1 = *iter;
					edge1Pos = edgePos;
					edge1V = edgeList[edgeId1].first;
					edge1U = edgeList[edgeId1].second;
					isNewEdges1 = edgesAdd.find(edgeId1) != edgesAddEnd;
				}
				else {
					edgeId2 = *iter;
					edge2Pos = edgePos;
					edge2V = edgeList[edgeId2].first;
					edge2U = edgeList[edgeId2].second;
					isNewEdges2 =  edgesAdd.find(edgeId2) != edgesAddEnd;
				}
			}
		}
		if (!isNewEdges1 && !isNewEdges2) return false;
		else {
			if (edge1V == edge2V || edge1V == edge2U || edge1U == edge2V || edge1U == edge2U) {
				generatedCC->subCCNum = 1;
				if (generatedCC->subCCs == nullptr) {
					generatedCC->subCCs = DBG_NEW vec(int)[1];
					generatedCC->subCCsaved = DBG_NEW bool[1];
				}
				else {
					generatedCC->subCCs[0].clear();
				}
				generatedCC->subCCs[0].emplace_back(edge1Pos);
				generatedCC->subCCs[0].emplace_back(edge2Pos);
				generatedCC->subCCsaved[0] = max(edge1Pos, edge2Pos) >= generatedCC->newInsert ? true : false;
				//Test::counter3 += END_NSTIMER(a2);
			}
			else {
				realSubCCNum = isNewEdges1 && isNewEdges2 ? 2 : 1;
				if (generatedCC->subCCs == nullptr) {
					generatedCC->subCCs = DBG_NEW vec(int)[realSubCCNum];
					generatedCC->subCCsaved = DBG_NEW bool[realSubCCNum];
				}
				else {
					generatedCC->subCCNum = (int)_msize(generatedCC->subCCsaved) / sizeof(bool);
					if (generatedCC->subCCNum < realSubCCNum) {
						delete[] generatedCC->subCCs;
						delete[] generatedCC->subCCsaved;
						generatedCC->subCCs = DBG_NEW vec(int)[realSubCCNum];
						generatedCC->subCCsaved = DBG_NEW bool[realSubCCNum];
					}
					if (realSubCCNum == 1) {
						generatedCC->subCCs[0].clear();
					}
					else {
						generatedCC->subCCs[0].clear();
						generatedCC->subCCs[1].clear();
					}
				}
				generatedCC->subCCNum = 0;
				if (isNewEdges1) {
					generatedCC->subCCs[generatedCC->subCCNum].emplace_back(edge1Pos);
					generatedCC->subCCsaved[generatedCC->subCCNum] = edge1Pos >= generatedCC->newInsert ? true : false;
					generatedCC->subCCNum++;
				}
				if (isNewEdges2) {
					generatedCC->subCCs[generatedCC->subCCNum].emplace_back(edge2Pos);
					generatedCC->subCCsaved[generatedCC->subCCNum] = edge2Pos >= generatedCC->newInsert ? true : false;
					generatedCC->subCCNum++;
				}
				//Test::counter3 += END_NSTIMER(a2);
			}
			return true;
		}
	}

	
	//BEGIN_NSTIMER(a3);
	subCCMap->clear();
	root2Id.clear();
	//get the number of new connected components
	int ufsetSize = min((int)size << 1, nNode);
	ds->reset(ufsetSize);
	int vertexNum = 0, sId, tId;
	subCCIter = &subCCId[0];
	int remainP = 0;
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
		if (*subCCIter != -1) {
			
			remainEdges[remainP++] = edgePos;
			auto tempEdge = &edgeList[*iter];

			int vertex = tempEdge->first;
			auto vertexIt = subCCMap->find(vertex);//O(1)
			if (vertexIt == subCCMap->end()) {
				sId = vertexNum;
				(*subCCMap)[vertex] = vertexNum++;
			}
			else sId = vertexIt->second;
			vertex = tempEdge->second;
			vertexIt = subCCMap->find(vertex);//O(1)
			if (vertexIt == subCCMap->end()) {
				tId = vertexNum;
				(*subCCMap)[vertex] = vertexNum++;
			}
			else tId = vertexIt->second;
			/*union-find operation means that sId and tId are connected*/
			ds->addE(*iter, sId, tId);
		}
	}
	//Test::counter+= END_NSTIMER(p);
	//BEGIN_NSTIMER(q);
	//get subCCNum
	subCCIter = &subCCId[0];
	haveNewEdgeCC.clear();
	bool allCCHaveNoNewEdges = false;
	for (int remainIter = 0; remainIter < remainP; remainIter++, subCCIter++) {
		int id = (*tempEdges)[remainEdges[remainIter]];
			/*the root of the edge's vertex in the disjoint set*/
			int root = ds->findRoot((*subCCMap)[edgeList[id].first]);
			auto ccIt = root2Id.find(root);
			if (ccIt == root2Id.end()) {//new subCC
				root2Id[root] = subCCNum;
				*subCCIter = subCCNum;
				++subCCNum;
			}
			else {
				*subCCIter = ccIt->second;
			}

			if (edgesAdd.find(id) != edgesAddEnd) {
				if (haveNewEdgeCC.find(*subCCIter) == haveNewEdgeCC.end()) {
					haveNewEdgeCC[*subCCIter] = realSubCCNum++;
				}
				allCCHaveNoNewEdges = true;
			}
	}

	if (realSubCCNum == 0) return false;
	//Test::counter2 += END_NSTIMER(q);

	//get subCCs
	if (generatedCC->subCCs == nullptr) {
		generatedCC->subCCs = DBG_NEW vec(int)[realSubCCNum];
		generatedCC->subCCsaved = DBG_NEW bool[realSubCCNum];
	}
	else {
		generatedCC->subCCNum = (int)_msize(generatedCC->subCCsaved) / sizeof(bool);
		if (generatedCC->subCCNum < realSubCCNum) {
			delete[] generatedCC->subCCs;
			delete[] generatedCC->subCCsaved;
			generatedCC->subCCs = DBG_NEW vec(int)[realSubCCNum];
			generatedCC->subCCsaved = DBG_NEW bool[realSubCCNum];
		}
		else {
			for (int i = 0; i < realSubCCNum; i++) {
				generatedCC->subCCs[i].clear();
			}
		}
	}
	generatedCC->subCCNum = realSubCCNum;

	subCCIter = &subCCId[0];
	auto newEdgesCCEnd = haveNewEdgeCC.end();
	//cout << Test::counter << " " << endTime << " " << Test::counter2 << endl;

	for (int remainIter = 0; remainIter < remainP; remainIter++, subCCIter++) {
		auto newEdgesCCIter = haveNewEdgeCC.find(*subCCIter);
		if (newEdgesCCIter != newEdgesCCEnd) {
			int realCCId = newEdgesCCIter->second;
			//cout << realCCId << " " << remainEdges[remainIter] << endl;
			generatedCC->subCCs[realCCId].emplace_back(remainEdges[remainIter]);
			generatedCC->subCCsaved[realCCId] = remainEdges[remainIter] >= generatedCC->newInsert ? true : false;//true: right non-expandable
		}
	}
	//Test::counter2++;
	//Test::counter5 += END_NSTIMER(a3);
	return allCCHaveNoNewEdges;
}

void TGraph::updateNewEdgeInfoOpt1(
	veciter(int)& infoBegin, veciter(int)& infoEnd, unordered_set<int>& edgesAdd,
	vec(CComponentsII*)& setCC,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, DisjointSet* tempDisjointSet,
	i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& checkedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue, vec(int) & saveCCPos, int startTime, int endTime, int k, i2iHMap& haveNewEdgeCC, int*& remainEdges, MaintainCCOPT1 mainStrategy) {
	if (checkedCC->first == nullptr && infoBegin == infoEnd) return;

	//BEGIN_TIMER(a)
	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int filterTime = startTime + k - 1;
	auto addEnd = edgesAdd.end();
	//auto beginTest = std::chrono::steady_clock::now();
	for (auto infoIter = infoBegin;
		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets

		id = *infoIter;//edge's id
		//cout << startTime<<" "<< endTime << " " << id << endl;
		/*the root of the edge's vertex in the disjoint set*/
		root = connectivity->findRoot(vertex2Pos[edgeList[id].first]);
		updateCCOpt1(setCC, root2Comp, edgesAdd, addEnd, /*lastEdgeSetR,*/ checkedCC, hasChecked, /*realMotifNum,*/ infoIter,
			id, root, filterTime, startTime, endTime, k);
	}
	//Test::counter += END_TIMER(a);


	LinkedNode<int>* tempCCIter = checkedCC->first;
	if (tempCCIter == nullptr) return;

	int size, queueSize;
	int ccPos;//cc's position in setCC
	/*the cc which has already generated/will generate*/
	CComponentsFRTMOPT1* generatedCC;
	veciter(int) motifEdgesEnd;
	veciter(int) subCCIter;

	i2iHMap_Iter vertexIt, ccIt;
	int removeEdge;
	vec(int)* tempEdges;
	//if (startTime == 3)exit(0);
	while (tempCCIter != nullptr) {

		//BEGIN_TIMER(b)

		ccPos = tempCCIter->item;
		generatedCC = (CComponentsFRTMOPT1*)setCC[ccPos];
		tempEdges = &generatedCC->edges;
		//cout << startTime << " " << endTime << " " << ccPos << endl;

		size = (int)tempEdges->size();
		if (generatedCC->newInsert == size) {//not updated edges
			if (generatedCC->tabuTChangePos != -1 && generatedCC->tabuTChangePos < endTime) {
				tempCCIter = tempCCIter->next;
				continue;
			}
		}

		removeEdge = 0;
		queueSize = (int)generatedCC->noisePosQueue.size();
		if (queueSize != 0) {
			CLEARALL(subCCId, 0, size, int);
		}

		generatedCC->tabuTChangePos = -1;
		tempRecordFromQueue.clear();
		while (queueSize != 0) {//worst time O(Em log (delta Em))
			auto ccLastNoisePair = generatedCC->noisePosQueue.top();

			if (ccLastNoisePair.first < endTime)
				break;
			else {
				if (ccLastNoisePair.first == endTime) {//edge has noise in this time
					//cout << startTime << " " << endTime << " " << ccPos << " " << generatedCC->edges[ccLastNoisePair.second.pos] << endl;

					removeEdge++;
					subCCId[ccLastNoisePair.second.pos] = -1;//remove this edge temporarily
					tempRecordFromQueue.emplace_back(ccLastNoisePair.second.startTime, // noise end position
						ccLastNoisePair.second.pos);
					if (generatedCC->tabuTChangePos < ccLastNoisePair.second.startTime - 1) {
						generatedCC->tabuTChangePos = ccLastNoisePair.second.startTime - 1;
					}
					/*if (ccPos == 13 && startTime == 42) {
						cout << endTime << " " << ccLastNoisePair.second.startTime << endl;
					}*/
					generatedCC->noisePosQueue.pop();
				}
				else {
					generatedCC->noisePosQueue.pop();
					if (endTime >= ccLastNoisePair.second.startTime) {
						//	cout << startTime << " " << endTime << " " << ccPos << " " << generatedCC->edges[ccLastNoisePair.second.pos] << endl;
						removeEdge++;
						subCCId[ccLastNoisePair.second.pos] = -1;//remove this edge temporarily
						tempRecordFromQueue.emplace_back(ccLastNoisePair.second.startTime, // noise end position
							ccLastNoisePair.second.pos);
						if (generatedCC->tabuTChangePos < ccLastNoisePair.second.startTime - 1) {
							generatedCC->tabuTChangePos = ccLastNoisePair.second.startTime - 1;
						}
					}
				}

				/*if (ccPos == 13 && startTime == 2) {
					cout << generatedCC->edges[ccLastNoisePair.second.pos] << "&"<< endl;
				}*/

			}
			queueSize = (int)generatedCC->noisePosQueue.size();
		}

		if (generatedCC->tabuTChangePos != -1) {
			if (generatedCC->noisePosQueue.size() != 0) {
				auto ccLastNoisePair = generatedCC->noisePosQueue.top();
				
				//cout << ccLastNoisePair.first << " " << ccLastNoisePair.second.pos << " " << endTime << endl;
				generatedCC->tabuTChangePos = max(ccLastNoisePair.first, generatedCC->tabuTChangePos);
			}

			int nextTime = endTime - 1;
			for (auto item : tempRecordFromQueue) {
				//cout << item.first << " " << item.second << endl;
				if (nextTime >= item.first) {
					generatedCC->noisePosQueue.emplace(
						nextTime, // noise next position
						NoisePos(
							item.first, // noise end position
							item.second // edge position in newCC->edges
						)
					);
				}
			}
		}

		while (generatedCC->preEMaxIntvlChangePos.size() != 0) {//worst time O(log (delta Em))
			auto preEMaxIntvlPair = generatedCC->preEMaxIntvlChangePos.top();

			if (preEMaxIntvlPair.endTime < endTime)
				break;
			else {
				id = preEMaxIntvlPair.id;
				scope[id].first = preEMaxIntvlPair.startTime;
				if (scope[id].second < preEMaxIntvlPair.endTime) {
					scope[id].second = preEMaxIntvlPair.endTime;
				}

				auto maxEMaxIntvlStartTPair = generatedCC->maxEMaxIntvlStartTQueue.top();
				if (maxEMaxIntvlStartTPair.second == id) {
					while (maxEMaxIntvlStartTPair.first != scope[id].first) {
						generatedCC->maxEMaxIntvlStartTQueue.pop();
						generatedCC->maxEMaxIntvlStartTQueue.emplace(scope[id].first, id);
						maxEMaxIntvlStartTPair = generatedCC->maxEMaxIntvlStartTQueue.top();
						id = maxEMaxIntvlStartTPair.second;
					}
				}

				generatedCC->preEMaxIntvlChangePos.pop();
			}
		}
		
		//Test::counter2 += END_TIMER(b);
		//BEGIN_TIMER(c)
		if (removeEdge > 0) {
			if (removeEdge == size) {
				tempCCIter = tempCCIter->next;
				generatedCC->newInsert = size;
				continue;
			}
			if ((this->*mainStrategy)(tempCCIter, tempRecordFromQueue, edgesAdd, addEnd,/* lastEdgeSetR,*/ vertex2Pos, tempDisjointSet, root2Comp, subCCMap, subCCId, root2Id, generatedCC, endTime, removeEdge,
				haveNewEdgeCC, remainEdges)) {
				saveCCPos.emplace_back(ccPos);
			}
			else {
				generatedCC->newInsert = size;
			}
			
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

#pragma endregion

#pragma region expandCheck
void TGraph::checkExpandableOpt1(int savePos, CComponentsII* temp,
	int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber) {
	CComponentsFRTMOPT1*tempCC = (CComponentsFRTMOPT1*)temp;
	TMotifII* motif;
	//newly generated connected components
	bool isRightExpandable;
	//if (Test::counter2 == 1) Test::counter3++;
	//else Test::counter = 2;
	//int edgeBefP = motifEndT * nEdge;
	pair<int, int> intersectIntv;
	//if (motifStartT == 41 && motifEndT == 172)exit(0);

	
	if (tempCC->subCCNum == 0) {//not edges removed  (must be non right expandable)
		if (motifStartT - startT != 0) {
			intersectIntv.first = tempCC->maxEMaxIntvlStartTQueue.top().first;
			intersectIntv.second = tempCC->minEMaxIntvlEndT;
			
			if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
				Test::gne121+= tempCC->edges.size();
				Test::gnefield++;
				//Test::gne121++;

				motif = DBG_NEW TMotifII(motifStartT, motifEndT);
				motif->copyEdges(tempCC->edges);
				result[savePos].emplace_back(motif);
				motifNumber++;

			}
			else {//left expandable or both expandable

				/*if (allNTimestamp > nEdge) {
					int maxP = 0x7fffffff, minP = -1;
					for (auto e : tempCC->edges) {
						int edgeP = edgeBefP + e;
						maxP = min(maxP, edgeBef[edgeP].first);
						minP = max(minP, edgeBef[edgeP].second);
						if (minP > maxP) break;
					}
					if (minP <= maxP && minP != -1 ) return;
				}*/
				Test::gne122 += tempCC->edges.size();
				//Test::gne122++;
				int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first);
				Test::gnefield += tempField * tempCC->edges.size();
				//Test::gnefield += tempField;
				Test::gnemaxfield = max(Test::gnemaxfield,tempField);
				
				bool isSave = true;

				if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
					//cout << Test::counter << " ! " << Test::counter2 << " " << Test::counter3 << " " << intersectIntv.first << " " << intersectIntv.second << endl;
					int j = intersectIntv.second;
					for (; j >= motifEndT; --j) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				else {
					int j = intersectIntv.first;
					for (; j <= motifStartT - 1; ++j) {
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
						if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				//cout << "$" << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << " "<< isSave<< endl;
				if (isSave) {
					//Test::counter++;
					//assert(intersectIntv.first != startT);
					motif = DBG_NEW TMotifII(motifStartT, motifEndT);
					motif->copyEdges(tempCC->edges);

					result[savePos].emplace_back(motif);
					motifNumber++;

				}
			}
		}
		else {
			Test::gne11+= tempCC->edges.size();
			//Test::gne11++;
			Test::gnefield++;

			motif = DBG_NEW TMotifII(motifStartT, motifEndT);
			motif->copyEdges(tempCC->edges);
			result[savePos].emplace_back(motif);
			motifNumber++;

		}

		/*if (allNTimestamp > nEdge) {
			for (auto e : tempCC->edges) {
				int edgeP = edgeBefP + e;
				if (edgeBef[edgeP].second + 1 == motifStartT) {
					if (motifStartT == 0) {
						edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
					}
					else edgeBef[edgeP].second = motifStartT;
				}
				else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
			}
		}*/
		//Test::counter +=END_NSTIMER;
	}
	else {
	
		int subCCNum = tempCC->subCCNum;
		
		//check each subCC
		veciter(int) subCCEnd;
		for (int i = 0; i < subCCNum; ++i) {
			//if (tempCC->subCCHaveNewEdges[i]) {//might expandable
			isRightExpandable = !tempCC->subCCsaved[i];//check whether exists e¡ÊR+[m,i]
			subCCEnd = tempCC->subCCs[i].end();
			if (!isRightExpandable) {//exists e¡ÊR+[m,i] (must be non right expandable)
				//BEGIN_NSTIMER; 
				if (motifStartT != 0) {

					getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);

					if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
						Test::gne2121 += tempCC->subCCs[i].size();
						//Test::gne2121++;
						Test::gnefield++;

						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							motif->addEdge(tempCC->edges[*edgeIter]);
						}
						result[savePos].emplace_back(motif);
						motifNumber++;

					}
					else {//left expandable or both expandable
						//if (allNTimestamp > nEdge) {
						//	int maxP = 0x7fffffff , minP = -1;
						//	for (auto p : tempCC->subCCs[i]) {
						//		auto e = tempCC->edges[p];
						//		int edgeP = edgeBefP + e;
						//		maxP = min(maxP, edgeBef[edgeP].first);
						//		minP = max(minP, edgeBef[edgeP].second);
						//		if (minP > maxP) break;
						//	}
						//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
						//	if (minP <= maxP && minP != -1) return;
						//}

						Test::gne2122 += tempCC->subCCs[i].size();
						//Test::gne2122++;
						int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first);
						Test::gnefield += tempField * tempCC->subCCs[i].size();
						//Test::gnefield += tempField;
						Test::gnemaxfield = max(Test::gnemaxfield, tempField);

						bool isSave = true;
						if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
							int j = intersectIntv.second;
							for (; j >= motifEndT; --j) {
								CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						else {
							int j = intersectIntv.first;
							for (; j <= motifStartT - 1; ++j) {
								CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
							/*if (motifStartT == 2 && ccPos->first == 13) {
								for (int i = motifEndT; i <= intersectIntv.second; i++) {
									cout << expandMask[i - motifEndT] << "^" << endl;
								}
							}*/
						}

						//cout << "$1 " << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << " " << isSave << endl;
						if (isSave) {
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							subCCEnd = tempCC->subCCs[i].end();
							for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
								motif->addEdge(tempCC->edges[*edgeIter]);
							}
							result[savePos].emplace_back(motif);
							motifNumber++;

						}
					}
				}
				else {
					Test::gne211 += tempCC->subCCs[i].size();
					//Test::gne211++;
					Test::gnefield++;

					motif = DBG_NEW TMotifII(motifStartT, motifEndT);
					subCCEnd = tempCC->subCCs[i].end();
					for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
						motif->addEdge(tempCC->edges[*edgeIter]);
					}
					result[savePos].emplace_back(motif);
					motifNumber++;

				}
			}
			else {
				//BEGIN_NSTIMER; 
				getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
				
				//cout << Test::counter << " # " << Test::counter2 << " " << Test::counter3 << " " << intersectIntv.first << " " << intersectIntv.second << endl;
				if (/*motifStartT == 0 || */intersectIntv.first == motifStartT) {//only check right expandable
					
					Test::gne221 += tempCC->subCCs[i].size();
					//Test::gne221++;
					int tempField;
					tempField = intersectIntv.second - motifEndT;
					Test::gnefield += tempField * tempCC->subCCs[i].size();
					//Test::gnefield += tempField;
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);
					CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
					if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
						//cout << "$2 " << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << " " << 1 << endl;
							
						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							motif->addEdge(tempCC->edges[*edgeIter]);
						}
						result[savePos].emplace_back(motif);
						motifNumber++;

					}
						//cout << "$3 " << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << " " << 0 << endl;



					//Test::counter += END_NSTIMER;
				}
				else {//left expandable, right expandable or both expandable
					//if (allNTimestamp > nEdge) {
					//	int maxP = 0x7fffffff, minP = -1;
					//	for (auto p : tempCC->subCCs[i]) {
					//		auto e = tempCC->edges[p];
					//		int edgeP = edgeBefP + e;
					//		maxP = min(maxP, edgeBef[edgeP].first);
					//		minP = max(minP, edgeBef[edgeP].second);
					//		if (minP > maxP) break;
					//	}
					//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
					//	if (minP <= maxP && minP != -1) return;
					//}

					Test::gne222 += tempCC->subCCs[i].size();
					//Test::gne222++;
					int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first + 1) - 1;
					Test::gnefield += tempField * tempCC->subCCs[i].size();
					//Test::gnefield += tempField;
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/

					bool isSave = true;
					if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
						int j = intersectIntv.second;
						for (; j > motifEndT; --j) {
							CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}
					}
					else {
						int j = intersectIntv.first;
						for (; j <= motifStartT; ++j) {
							CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
							if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
								break;
							}
						}

						/*if (!isSave && motifStartT == 2 && ccPos->first == 13) {
							for (int s = motifEndT + 1; s <= intersectIntv.second; s++) {
								cout << j << " " << s << " " << expandMask[s - motifEndT - 1] << "^" << endl;
							}
						}*/
					}

					if (isSave) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
						}
						/*if (!isSave && motifStartT == 2 && ccPos->first == 13) {
							for (int i = intersectIntv.first; i <= motifStartT - 1; i++) {
								cout << i<< " "<<expandMask[i - intersectIntv.first] << "^" << endl;
							}
						}*/
					}

					//cout << "$4 " << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << " " << isSave << endl;
						
					if (isSave) {
						//Test::counter3++;
						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							motif->addEdge(tempCC->edges[*edgeIter]);
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}

					//Test::counter2 += END_NSTIMER;
				}
				//}
			}
			//}

			/*if (allNTimestamp > nEdge) {
				for (auto p : tempCC->subCCs[i]) {
					auto e = tempCC->edges[p];
					int edgeP = edgeBefP + e;
					if (edgeBef[edgeP].second + 1 == motifStartT) {
						if (motifStartT == 0) {
							edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
						}
						else edgeBef[edgeP].second = motifStartT;
					}
					else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
				}
			}*/
		}
		//}
	}
}

void TGraph::checkExpandableOpt1MidR(int savePos, CComponentsII* temp,
	int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber) {
	CComponentsFRTMOPT1*tempCC = (CComponentsFRTMOPT1*)temp;
	TMotifII* motif;
	bool needCheck;
	int edgeId;
	bool isRightExpandable;
	//int edgeBefP = motifEndT * nEdge;
	pair<int, int> intersectIntv;
	auto edgeEnd = tempCC->edges.end(); 
	if (tempCC->subCCNum == 0) {//not edges removed  (must be non right expandable)
		if (motifStartT - startT != 0) {
			intersectIntv.first = tempCC->maxEMaxIntvlStartTQueue.top().first;
			intersectIntv.second = tempCC->minEMaxIntvlEndT;

			if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
				//Test::gne121++;

				motif = DBG_NEW TMotifII(motifStartT, motifEndT);
				needCheck = true;
				for (auto edgeIter = tempCC->edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
					edgeId = *edgeIter;
					if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) { needCheck = false; break; }
				}
				if (needCheck) {
					newMIntR->emplace_back(savePos, (int)result[savePos].size());
				}
				motif->copyEdges(tempCC->edges);
				result[savePos].emplace_back(motif);
				motifNumber++;
			}
			else {//left expandable or both expandable
				//Test::gne122++;
				//int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first + 1);
				//Test::gnefield += tempField;
				//Test::gnemaxfield = max(Test::gnemaxfield,tempField);
				//cout << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << endl;
				/*if (allNTimestamp > nEdge) {
					int maxP = 0x7fffffff, minP = -1;
					for (auto e : tempCC->edges) {
						int edgeP = edgeBefP + e;
						maxP = min(maxP, edgeBef[edgeP].first);
						minP = max(minP, edgeBef[edgeP].second);
						if (minP > maxP) break;
					}
					if (minP <= maxP && minP != -1) return;
				}*/
				bool isSave = true;

				if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
					//cout << Test::counter << " ! " << Test::counter2 << " " << Test::counter3 << " " << intersectIntv.first << " " << intersectIntv.second << endl;
					int j = intersectIntv.second;
					for (; j >= motifEndT; --j) {
						CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				else {
					int j = intersectIntv.first;
					for (; j <= motifStartT - 1; ++j) {
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
						if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							isSave = false;
							break;
						}
					}
				}
				if (isSave) {
					motif = DBG_NEW TMotifII(motifStartT, motifEndT);
					needCheck = true;
					for (auto edgeIter = tempCC->edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
						edgeId = *edgeIter;
						if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) { needCheck = false; break; }
					}
					if (needCheck) {
						newMIntR->emplace_back(savePos, (int)result[savePos].size());
					}
					motif->copyEdges(tempCC->edges);

					result[savePos].emplace_back(motif);
					motifNumber++;
				}
			}
		}
		else {
			//Test::gne11++; 

			motif = DBG_NEW TMotifII(motifStartT, motifEndT);
			needCheck = true;
			for (auto edgeIter = tempCC->edges.begin(); edgeIter != edgeEnd; ++edgeIter) {//O(motif number)
				edgeId = *edgeIter;
				if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) { needCheck = false; break; }
			}
			if (needCheck) {
				newMIntR->emplace_back(savePos, (int)result[savePos].size());
			}
			motif->copyEdges(tempCC->edges);
			result[savePos].emplace_back(motif);
			motifNumber++;
		}

		/*if (allNTimestamp > nEdge) {
			for (auto e : tempCC->edges) {
				int edgeP = edgeBefP + e;
				if (edgeBef[edgeP].second + 1 == motifStartT) {
					if (motifStartT == 0) {
						edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
					}
					else edgeBef[edgeP].second = motifStartT;
				}
				else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
			}
		}*/
		//Test::counter +=END_NSTIMER;
	}
	else {
		int subCCNum = tempCC->subCCNum;
		//check each subCC
		veciter(int) subCCEnd;
		for (int i = 0; i < subCCNum; ++i) {
			//if (tempCC->subCCHaveNewEdges[i]) {//might expandable
				isRightExpandable = !tempCC->subCCsaved[i];//check whether exists e¡ÊR+[m,i]
				subCCEnd = tempCC->subCCs[i].end();
				
				if (!isRightExpandable) {//exists e¡ÊR+[m,i] (must be non right expandable)
					//BEGIN_NSTIMER; 
					if (motifStartT != 0) {

						getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
						if (intersectIntv.first == motifStartT) {//non left expandable => non both expandable
							//Test::gne2121++;

							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							subCCEnd = tempCC->subCCs[i].end();

							needCheck = true;
							for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
								edgeId = tempCC->edges[*edgeIter]; 
								motif->addEdge(edgeId);
								if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
							}
							if (needCheck) {
								newMIntR->emplace_back(savePos, (int)result[savePos].size());
							}
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
						else {//left expandable or both expandable
							//Test::gne2122++;
							//int tempField = (intersectIntv.second - motifEndT+1) * (motifStartT - intersectIntv.first + 1);
							//Test::gnefield += tempField;
							//Test::gnemaxfield = max(Test::gnemaxfield, tempField);
							//cout << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << endl;
							//if (allNTimestamp > nEdge) {
							//	int maxP = 0x7fffffff, minP = -1;
							//	for (auto p : tempCC->subCCs[i]) {
							//		auto e = tempCC->edges[p];
							//		int edgeP = edgeBefP + e;
							//		maxP = min(maxP, edgeBef[edgeP].first);
							//		minP = max(minP, edgeBef[edgeP].second);
							//		if (minP > maxP) break;
							//	}
							//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
							//	if (minP <= maxP && minP != -1) return;
							//}

							bool isSave = true;
							if (motifStartT - intersectIntv.first >= intersectIntv.second - motifEndT + 1) {
								int j = intersectIntv.second;
								for (; j >= motifEndT; --j) {
									CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
									if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
										isSave = false;
										break;
									}
								}
							}
							else {
								int j = intersectIntv.first;
								for (; j <= motifStartT - 1; ++j) {
									CLEARALL(expandMask, true, intersectIntv.second - motifEndT + 1, bool);
									if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
										isSave = false;
										break;
									}
								}
							}

							if (isSave) {
								motif = DBG_NEW TMotifII(motifStartT, motifEndT);
								subCCEnd = tempCC->subCCs[i].end();
								needCheck = true;
								for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
									edgeId = tempCC->edges[*edgeIter];
									motif->addEdge(edgeId);
									if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
								}
								if (needCheck) {
									newMIntR->emplace_back(savePos, (int)result[savePos].size());
								}
								result[savePos].emplace_back(motif);
								motifNumber++;
							}
						}
					}
					else {
						//Test::gne211++;

						motif = DBG_NEW TMotifII(motifStartT, motifEndT);
						subCCEnd = tempCC->subCCs[i].end();
						needCheck = true;
						for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
							edgeId = tempCC->edges[*edgeIter];
							motif->addEdge(edgeId);
							if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
						}
						if (needCheck) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}
						result[savePos].emplace_back(motif);
						motifNumber++;
					}
					//Test::counter += END_NSTIMER;
				}
				else {
					//BEGIN_NSTIMER; 
					getIntersectionOfScope(tempCC->subCCs[i], tempCC->edges, intersectIntv);
					if (motifStartT == 0 || intersectIntv.first == motifStartT) {//only check right expandable

						//Test::gne221++;
						//int tempField;
						//if (intersectIntv.second - motifEndT - 1 > 0) {
						//	tempField = (intersectIntv.second - motifEndT - 1) * (motifStartT - intersectIntv.first + 1);
							//cout << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << endl;
						//}else tempField = 0;
						//Test::gnefield += tempField;
						//Test::gnemaxfield = max(Test::gnemaxfield, tempField);

						//if (motifStartT == 0 && intersectIntv.second == endT) continue;//right expandable
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
						if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)

							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							subCCEnd = tempCC->subCCs[i].end();
							needCheck = true;
							for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
								edgeId = tempCC->edges[*edgeIter];
								motif->addEdge(edgeId);
								if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
							}
							if (needCheck) {
								newMIntR->emplace_back(savePos, (int)result[savePos].size());
							}
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
						//Test::counter += END_NSTIMER;
					}
					else {//left expandable, right expandable or both expandable

						//Test::gne222++;
						//int tempField = (intersectIntv.second - motifEndT) * (motifStartT - intersectIntv.first);
						//Test::gnefield += tempField;
						//Test::gnemaxfield = max(Test::gnemaxfield, tempField);
						//cout << tempField << " " << intersectIntv.first << " " << intersectIntv.second << " " << motifStartT << " " << motifEndT << endl;
						//if (allNTimestamp > nEdge) {
						//	int maxP = 0x7fffffff, minP = -1;
						//	for (auto p : tempCC->subCCs[i]) {
						//		auto e = tempCC->edges[p];
						//		int edgeP = edgeBefP + e;
						//		maxP = min(maxP, edgeBef[edgeP].first);
						//		minP = max(minP, edgeBef[edgeP].second);
						//		if (minP > maxP) break;
						//	}
						//	//if (motifStartT == 41 && motifEndT == 173) cout << minP << " @ " << maxP << endl;
						//	if (minP <= maxP && minP != -1) return;
						//}

						bool isSave = true;
						if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
							int j = intersectIntv.second;
							for (; j > motifEndT; --j) {
								CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;

									break;
								}
							}
						}
						else {
							int j = intersectIntv.first;
							for (; j <= motifStartT; ++j) {
								CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}

						if (isSave) {
							CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->subCCs[i], tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
							}

						}

						if (isSave) {
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							subCCEnd = tempCC->subCCs[i].end();
							needCheck = true;
							for (auto edgeIter = tempCC->subCCs[i].begin(); edgeIter != subCCEnd/*edgeEnd*/; ++edgeIter/*, ++subCCIter*/) {//O(motif number)
								edgeId = tempCC->edges[*edgeIter];
								motif->addEdge(tempCC->edges[*edgeIter]);
								if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == -3) needCheck = false;
							}
							if (needCheck) {
								newMIntR->emplace_back(savePos, (int)result[savePos].size());
							}
							result[savePos].emplace_back(motif);
							motifNumber++;
						}

						//Test::counter2 += END_NSTIMER;
					}
				}

				/*if (allNTimestamp > nEdge) {
					for (auto p : tempCC->subCCs[i]) {
						auto e = tempCC->edges[p];
						int edgeP = edgeBefP + e;
						if (edgeBef[edgeP].second + 1 == motifStartT) {
							if (motifStartT == 0) {
								edgeBef[edgeP].first = edgeBef[edgeP].second = 0;
							}
							else edgeBef[edgeP].second = motifStartT;
						}
						else edgeBef[edgeP].first = edgeBef[edgeP].second = motifStartT;
					}
				}*/
		}
	}
}

void TGraph::expCheckFRTM(vec(int)&saveCCPos,
	vec(CComponentsII*)& setCC,
	vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
	long long& motifNumber, CheckExpandable checkExpandable) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	auto saveMotifEnd = saveCCPos.end();
	CComponentsII* tempCC;
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {//need check expandable
		//Util::output(motifStartT, motifEndT, *saveMotifIter);
		auto a = std::chrono::steady_clock::now();
		tempCC = setCC[*saveMotifIter];
		
		(this->*checkExpandable)(savePos, tempCC, motifStartT, motifEndT, expandMask, /*expCheck, saveMotifIter,*/ result, motifNumber);

		int edgeNum = (int)tempCC->edges.size();
		//cout << saveMotifIter->first <<  " @ "<<std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count()<< endl;
		
		if (edgeNum > tempCC->newInsert) {
			tempCC->newInsert = edgeNum;
			//expCheck.erase(saveMotifIter->first);
		}
	}
	saveCCPos.clear();
}

void TGraph::expCheckFRTMPlus(vec(int)&saveCCPos,
	vec(CComponentsShortIntv*)& setCCShortIntv,
	vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
	long long& motifNumber, DynamicConnectivity*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked, i2iHMap& root2Comp, i2iHMap& vertex2Pos, vec(CComponentsII*)& setCC) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	veciter(int) saveMotifEnd = saveCCPos.end();
	TMotifII* motif;
	CComponentsShortIntv* tempCC;
	pair<int, int> intersectIntv;
	auto mapEnd = root2Comp.end();
	auto checkedEnd = hasChecked.end();
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		tempCC = setCCShortIntv[*saveMotifIter];
		if (tempCC->startT == motifStartT) {//left unexpandable
			intersectIntv.first = tempCC->scopeL;
			intersectIntv.second = tempCC->scopeR;
			
			if (intersectIntv.second == -1) {//only in [motifStartT, motifStartT+limited)
				Test::gnenonoise += tempCC->edges.size();
				//Test::gnenonoise++;
				//Test::gnefield++;

				motif = DBG_NEW TMotifII(motifStartT, motifEndT);
				motif->copyEdges(tempCC->edges);
				result[savePos].emplace_back(motif);
				motifNumber++;
			}
			else {
				bool unsavedBefore = false;
				if (!tempCC->haveNonExpand) {
					int firstENode = edgeList[*tempCC->edges.begin()].first;
					int cc = root2Comp[connectivity->findRoot(vertex2Pos[firstENode])];
					CComponentsFRTMOPT1* checkedCC = (CComponentsFRTMOPT1*)setCC[cc];
					auto checkIter = hasChecked.find(cc);
					unsavedBefore = checkIter == checkedEnd || (checkedCC->tabuTChangePos != -1 && checkedCC->tabuTChangePos < motifEndT);
				}
				
				if (tempCC->haveNonExpand || !unsavedBefore) {//check expandable
					Test::gnenonoise+= tempCC->edges.size();
					//Test::gnenonoise++;
					if (intersectIntv.first == motifStartT) {//only check right expandable
						int tempField;
						if (intersectIntv.second - motifEndT - 1 > 0) {
							tempField = intersectIntv.second - motifEndT - 1;
						}
						else tempField = 0;
						Test::gnefield += tempField * tempCC->edges.size();
						//Test::gnefield += tempField;
						Test::gnemaxfield = max(Test::gnemaxfield, tempField);
						//Test::counter4 += tempField;
						//Test::counter5++;
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
						if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							motif->copyEdges(tempCC->edges);
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
					}
					else {//left expandable, right expandable or both expandable
						//Test::gnenonoise+=tempCC->edges.size();
						int tempField = (intersectIntv.second - motifEndT) * (motifStartT - intersectIntv.first);
						Test::gnefield += tempField * tempCC->edges.size();
						//Test::gnefield += tempField;
						Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/
						
						bool isSave = true;
						if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
							int j = intersectIntv.second;
							for (; j > motifEndT; --j) {
								CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						else {
							int j = intersectIntv.first;
							for (; j <= motifStartT; ++j) {
								CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						if (isSave) {
							CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
							}
							if (isSave) {
								motif = DBG_NEW TMotifII(motifStartT, motifEndT);
								motif->copyEdges(tempCC->edges);
								result[savePos].emplace_back(motif);
								motifNumber++;
							}
						}
					}
				}
			}
		}
	}
	saveCCPos.clear();
}

void TGraph::expCheckFRTMPlusMidR(vec(int)&saveCCPos,
	vec(CComponentsShortIntv*)& setCCShortIntv,
	vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
	long long& motifNumber, DynamicConnectivity*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked, i2iHMap& root2Comp, i2iHMap& vertex2Pos, vec(CComponentsII*)& setCC) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	veciter(int) saveMotifEnd = saveCCPos.end();
	TMotifII* motif;
	CComponentsShortIntv* tempCC;
	pair<int, int> intersectIntv;
	bool needCheck;
	auto mapEnd = root2Comp.end();
	auto checkedEnd = hasChecked.end();
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		tempCC = setCCShortIntv[*saveMotifIter];
		
		if (tempCC->startT == motifStartT) {//left unexpandable
			intersectIntv.first = tempCC->scopeL;
			intersectIntv.second = tempCC->scopeR;
			if (intersectIntv.second == -1) {//only in [motifStartT, motifStartT+limited)
				motif = DBG_NEW TMotifII(motifStartT, motifEndT);

				needCheck = true;
				for (auto edgeId : tempCC->edges) {//O(motif number)
					if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == EMaxIntvlChange::INTVINIT) { needCheck = false; break; }
				}
				if (needCheck) {
					newMIntR->emplace_back(savePos, (int)result[savePos].size());
				}

				motif->copyEdges(tempCC->edges);
				result[savePos].emplace_back(motif);
				motifNumber++;
			}
			else {
				bool unsavedBefore = false;
				if (!tempCC->haveNonExpand) {
					int firstENode = edgeList[*tempCC->edges.begin()].first;
					int cc = root2Comp[connectivity->findRoot(vertex2Pos[firstENode])];
					CComponentsFRTMOPT1* checkedCC = (CComponentsFRTMOPT1*)setCC[cc];
					auto checkIter = hasChecked.find(cc);
					unsavedBefore = checkIter == checkedEnd || (checkedCC->tabuTChangePos != -1 && checkedCC->tabuTChangePos < motifEndT);
				}
				
				if (tempCC->haveNonExpand || !unsavedBefore) {//check expandable
					
					if (intersectIntv.first == motifStartT) {//only check right expandable
						CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
						if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);
							
							needCheck = true;
							for (auto edgeId : tempCC->edges) {//O(motif number)
								if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == EMaxIntvlChange::INTVINIT) { needCheck = false; break; }
							}
							if (needCheck) {
								newMIntR->emplace_back(savePos, (int)result[savePos].size());
							}

							motif->copyEdges(tempCC->edges);
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
					}
					else {//left expandable, right expandable or both expandable
						bool isSave = true;
						if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
							int j = intersectIntv.second;
							for (; j > motifEndT; --j) {
								CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
								if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						else {
							int j = intersectIntv.first;
							for (; j <= motifStartT; ++j) {
								CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
								if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
									isSave = false;
									break;
								}
							}
						}
						if (isSave) {
							CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
							if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
								isSave = false;
							}
							if (isSave) {
								motif = DBG_NEW TMotifII(motifStartT, motifEndT);

								needCheck = true;
								for (auto edgeId : tempCC->edges) {//O(motif number)
									if (needCheck&& newPosInEIntR[edgeId*numOfLabel + getEdgeLabel(edgeId, motifStartT)] == EMaxIntvlChange::INTVINIT) { needCheck = false; break; }
								}
								if (needCheck) {
									newMIntR->emplace_back(savePos, (int)result[savePos].size());
								}

								motif->copyEdges(tempCC->edges);
								result[savePos].emplace_back(motif);
								motifNumber++;
							}
						}
					}
				}
			}
		}
	}
	saveCCPos.clear();
}


#pragma region FRTMPlus
void TGraph::maintainConnectivityPlus(veciter(int)& infoBegin,
	veciter(int)& infoEnd, int*& addE, int& addENum, i2iHMap& vertex2Pos, DynamicConnectivity*& connectivityShortIntv, int&vertexNum,
	bool*& hasE, bool*& saveR, DynamicConnectivity*& connectivity,
	vec(int)& combineCCPos) {
	NodePair* temp;
	int motifPos;
	i2iHMap_Iter vertexIt;
	int sId, tId, vertex;//vertexs'id in the disjoint set

	for (veciter(int) iter = infoBegin;
		iter != infoEnd; ++iter) {
		temp = &edgeList[*iter];//new edge from one R edge set
		vertex = temp->first;
		if (hasE[*iter]) {
			assert(vertex2Pos.find(vertex) != vertex2Pos.end());
			int r = connectivity->findRoot(vertex2Pos[vertex]);
			if (saveR[r]) {
				addE[addENum++] = iter - infoBegin;
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

				motifPos = connectivityShortIntv->addE(*iter, sId, tId);

				if (motifPos != -1)
					combineCCPos.emplace_back(motifPos);
			}
		}
		else {
			addE[addENum++] = iter - infoBegin;
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

			motifPos = connectivityShortIntv->addE(*iter, sId, tId);

			if (motifPos != -1)
				combineCCPos.emplace_back(motifPos);
		}
	}
}

void TGraph::combineComponentsPlus(vec(CComponentsShortIntv*)& tempComponents,
	i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity,
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
	i2iHMap& vertex2PosShortIntv, DynamicConnectivity*& connectivity,
	i2iHMap& root2Comp, bool*& hasE,
	vec(int)& maxCC, int startTime, int endTime) {

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
	auto vEnd = vertex2PosShortIntv.end();
	int count = 0;
	for (int addEIter = 0; addEIter < addENum; ++addEIter) {//new edges in R edge sets
		veciter(int) infoIter = infoBegin + addE[addEIter];
		id = *infoIter;//edge's id
		//cout << startTime << " " << endTime << " " << id << endl;

		/*the root of the edge's vertex in the disjoint set*/
		root = connectivity->findRoot(vertex2PosShortIntv[edgeList[id].first]);

		mapIter = root2Comp.find(root);//check which cc the edge belongs to 
		eStartT = maxIntvShortIntv[id].first;//start time of the edge
		label = getEdgeLabel(id, eStartT); //label of the edge

		nonexpand = !hasE[id];

		if (mapIter == root2Comp.end()) {//generateMaxTM case 1
			newCC = DBG_NEW CComponentsShortIntv(eStartT, root);
			newCC->edges.emplace_back(id);
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
			generatedCC->edges.emplace_back(id);

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
				//if (generatedCC->startT == startTime) {//check left unexpandable
				maxCC.emplace_back(ccPos);
				hasSaved[ccPos] = true;
				//}
			}
		}
	}
}

void TGraph::recordSavedRoot(vec(CComponentsII*)& setCC, unordered_map<int, LinkedNode<int>*>& hasChecked, int*& newE, int& newENum, bool*& saveR, DynamicConnectivity*& connectivity, i2iHMap& vertex2Pos) {
	for (auto cc : hasChecked) {
		saveR[setCC[cc.first]->root] = true;
	}
	auto vEnd = vertex2Pos.end();
	for (auto i = 0; i < newENum; i++) {
		int e = newE[i];
		auto v = vertex2Pos.find(edgeList[e].first);
		if (v != vEnd) {
			int r = connectivity->findRoot(v->second);
			if (r != -1) {
				saveR[r] = true;
			}
		}
		v = vertex2Pos.find(edgeList[e].second);
		if (v != vEnd) {
			int r = connectivity->findRoot(v->second);
			if (r != -1) {
				saveR[r] = true;
			}
		}
	}
}

#pragma endregion

void TGraph::findRTMotifsDynamic(int k,
	vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed) {
#pragma region intialization
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(CComponentsII*) setCC;//temporary component list CC
	vec(int)/*SAVESIMINFO_Vec*/* edgeSetsR = DBG_NEW vec(int)[endT - oriEndT + 1];// R edge sets
	DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode + 1);//disjoint set

	i2iHMap subCCMap, root2Id;
	DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);
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
	int vertexNum;//the number of vertexs
	//int realMotifNum;//the number of generated connected component
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
	//i2tupHMap expCheck;
	int id, Tsk = startT + k - 1;
	//traverse all possible intervals' start time
	//TGraph::tempUsedForComputeRESDYN = (oriEndT+1)*nEdge;
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
		//cout << Ts << " currentComputeRES" << endl;
		this->edgeFilterFRTMMidRForDYN(Ts, Te, oriEndTe,
			edgeSetsR, selectedNum, rightEndpoint/*, EdgeMaskLen*/,
			k, fixLabel, isEdgeTypeFixed);
			//initalize for every row
		Test::compr += END_TIMER(b);
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;


		disjointSet->reset();
		vertex2Pos.clear();
		root2Comp.clear();
		hasChecked.clear();
		checkedCC->release();
		vertexNum = 0;
		//realMotifNum = 0;
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
					//cout << Ts << " " << tempT << " maintainConnectivity" << endl;
					/*fetch new edges from one R edge set and insert into the disjoint set (maintain connectivity)*/
					maintainConnectivity(iterStart, iterEnd,
						vertex2Pos, disjointSet, vertexNum,
						combineMotifPos);
					/*combine components which are connected after adding new edges from one R edge set,
						combine components before adding new edges for less computation
						(combine components in generateMaxTM case 3)*/
					combineComponentsFRTM(setCC,
						vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
						/*setCCSize, realMotifNum,*/
						combineMotifPos);
				}
			//cout << Ts << " " << tempT << " updateNewEdgeInfoD4" << endl;

			/*add edge into generated motif or generate new motif and check left expandable
				(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)
			some ccs need to be saved though not new edges fetched*/
			updateNewEdgeInfoFRTM(iterStart, iterEnd,
				setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, /*expandMask,
				 EdgeMaskLen,setCCSize,*/tempRecordFromQueue,
				 ///*realMotifNum,*/ saveMotifPos, Ts, edgeEndT, k, &TGraph::maintainTempCCForUF, &TGraph::updateCCD5, ComponentsD5Type::CCFRTM);//O(¦¤Es)
				/*realMotifNum,*/ saveMotifPos, Ts, edgeEndT, k, &TGraph::maintainTempCCForUF, &TGraph::updateCCFRTM, ComponentsD5Type::CCFRTM);//O(¦¤Es)
			Test::gm += END_TIMER(c);
			//cout << Ts << " " << tempT << " " << Test::gm << endl;

#pragma endregion

			BEGIN_TIMER(d)
				expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
					edgeEndT, expandMask,/* expCheck,*/ motifNumber, &TGraph::checkExpandableFRTMMidR);
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
				auto motifEdges = newResult[resultP][resultIndex]->getMotifEdge();
				auto motifEdgesEnd = motifEdges->end();
				auto motifEdgesIter = motifEdges->begin();
				bool flag = false;
				for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
					id = *motifEdgesIter;
					int inMidP = id * numOfLabel + getEdgeLabel(id, motifEndT);
					if (posInEIntR[inMidP] != EMaxIntvlChange::CHANGED) break;
				}
				//cout << Ts << " " << motifEndT << " " << resultIndex << endl;

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

					//if (Ts == 176 && motifEndT == 198) cout << maxLeft << " @ " << minRight << endl;
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
	long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed) {
#pragma region intialization
	i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp(nNode);

	vec(CComponentsII*) setCC;//temporary component list CC
	int realSize = endT - oriEndT + 1;
	vec(int)/*SAVESIMINFO_Vec*/* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
	vec(int) *edgeSetsRAdd = DBG_NEW vec(int)[realSize]/*, *edgeSetsRDel = DBG_NEW vec(int)[realSize]*/;
	unordered_set<int> edgesAdd(nEdge<<1);
	DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
	DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);

	i2iHMap subCCMap(nNode), root2Id(nNode);
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
	unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map ccpos in setCC to previous node in checkedCC
	int vertexNum = 0;//the number of vertexs
	//int realMotifNum;//the number of generated connected component
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
	i2iHMap haveNewEdgeCC;
	int* remainEdges = DBG_NEW int[nEdge];
#pragma endregion
	vector<pair<int, int>> tempRecordFromQueue;
	//i2tupHMap expCheck;

	int id, Tsk = startT + k - 1;
	//traverse all possible intervals' start time
	//TGraph::tempUsedForComputeRESDYN = (oriEndT+1)*nEdge;
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
		//cout << Ts << " currentComputeRES" << endl;
		edgeFilterOpt1MidRForDYN(Ts, Te, oriEndTe, edgeSetsRAdd,/**/
			edgeSetsR, selectedNum, rightEndpoint/*, EdgeMaskLen*/,
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
		//realMotifNum = 0;
		edgeEndT = endT;

		//Test S(m,T)
		smt = edgeSetsR[endT - oriEndTe].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;
		smk = edgeSetsR[0].size();
		maxSmk = max(maxSmk, smk);
		allSmk += smk;

		//R2L
		//if(Ts > changeStartPosT) startPos--;
		//int stopT = oriEndTe - Te;
		edgesAdd.clear();
		for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {

			if (tempT > rightEndpoint) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			
			//maintain sets which need to add and delete O(delta |R|)
			for (auto add : edgeSetsRAdd[tempT]) {
				edgesAdd.emplace(add);
			}
			edgeSetsRAdd[tempT].clear();
			BEGIN_TIMER(c)
#pragma region maxCheck

				if (iterStart != iterEnd) {//new edges fetched
					//cout << Ts << " " << tempT << " maintainConnectivity" << endl;
					maintainConnectivity(iterStart, iterEnd,
						vertex2Pos, disjointSet, vertexNum,
						combineMotifPos);
					combineComponentsOpt1(setCC,
						vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
						/*setCCSize, realMotifNum,*/
						combineMotifPos);
				}
			//cout << Ts << " " << tempT << " updateNewEdgeInfoD4" << endl;

			updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,/**/
				setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, /*expandMask,
				 EdgeMaskLen,setCCSize,*/tempRecordFromQueue,
				 ///*realMotifNum,*/ saveMotifPos, Ts, edgeEndT, k, &TGraph::maintainTempCCForUF, &TGraph::updateCCD5, ComponentsD5Type::CCFRTM);//O(¦¤Es)
				/*realMotifNum,*/ saveMotifPos, Ts, edgeEndT, k, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
			Test::gm += END_TIMER(c);
			//cout << Ts << " " << tempT << " " << Test::gm << endl;

#pragma endregion

			BEGIN_TIMER(d)
				expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
					edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1MidR);
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
				auto motifEdges = newResult[resultP][resultIndex]->getMotifEdge();
				auto motifEdgesEnd = motifEdges->end();
				auto motifEdgesIter = motifEdges->begin();
				bool flag = false;
				for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
					id = *motifEdgesIter;
					int inMidP = id * numOfLabel + getEdgeLabel(id, motifEndT);
					
					if (posInEIntR[inMidP] != EMaxIntvlChange::CHANGED) break;
				}
				//cout << Ts << " " << motifEndT << " " << resultIndex << endl;

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
	delete disjointSet;
	delete[] subCCId;
	delete tempDisjointSet;
	delete[] expandMask;
	delete[] edgeSetsR;
	delete[] edgeSetsRAdd;
	delete[] remainEdges;
	delete checkedCC;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime + 1) << endl;
	cout << "maxSmk: " << maxSmk << " averSmk: " << allSmk / (lastTime + 1) << endl;
}

void TGraph::findRTMotifsDynamicPlus(int k,
	vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed) {
	//for no noise field
	int newk = (int)(ceil(1 / Setting::delta) + 0.5);
	int maxLengthForNoNoise = newk - k;
	if (maxLengthForNoNoise >= 1) {
#pragma region intialization
		i2iHMap vertex2Pos(nNode << 1);//map the vertex's id to the position in disjoint set
		/*map the root in disjoint set to the position of its corresponding
		connected component in the component list*/
		i2iHMap root2Comp(nNode);

		vec(CComponentsII*) setCC;//temporary component list CC
		int realSize = endT - oriEndT + 1;
		vec(int)/*SAVESIMINFO_Vec*/* edgeSetsR = DBG_NEW vec(int)[realSize];// R edge sets
		vec(int) *edgeSetsRAdd = DBG_NEW vec(int)[realSize]/*, *edgeSetsRDel = DBG_NEW vec(int)[realSize]*/;
		unordered_set<int> edgesAdd(nEdge << 1);
		DynamicConnectivity* disjointSet = DBG_NEW DisjointSet(nNode);//disjoint set
		DisjointSet* tempDisjointSet = DBG_NEW DisjointSet(nNode);

		i2iHMap subCCMap(nNode), root2Id(nNode);
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
		unordered_map<int, LinkedNode<int>*> hasChecked(nEdge); // map ccpos in setCC to previous node in checkedCC
		int vertexNum = 0;//the number of vertexs
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
		i2iHMap haveNewEdgeCC;
		int* remainEdges = DBG_NEW int[nEdge];
#pragma endregion
		vector<pair<int, int>> tempRecordFromQueue;

		vec(int)* edgeSetsRNoNoise = DBG_NEW vec(int)[maxLengthForNoNoise];// R edge sets
		i2iHMap root2CompNoNoise(nNode);
		vec(CComponentsShortIntv*) setCCShortIntv;//temporary component list CC
		DynamicConnectivity* disjointSetNoNoise = DBG_NEW DisjointSet(nNode);//disjoint set
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
		//iSet saveR(nNode);
		bool* saveR = DBG_NEW bool[nNode];
		int id, Tsk = startT + k - 1;
		//i2tupHMap expCheck;

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
						disjointSet->reset(vertexNum);
						hasChecked.clear();
						checkedCC->release();

						disjointSetNoNoise->reset(vertexNum);
						root2CompNoNoise.clear();
					}
					edgeFilterPlusMidRForDYN(Ts, Te, TeNoNoise, Te - oriEndTe - 1, oriEndTe, hasE, newE, newENum, edgeSetsRAdd,
						edgeSetsR, selectedNum, edgeSetsRNoNoise, selectedNumNoNoise, rightEndpoint, k, fixLabel, isEdgeTypeFixed);
				}
				else {//only no distortion
					if (Ts != startT) {
						disjointSetNoNoise->reset(vertexNum);
						root2CompNoNoise.clear();
					}
					edgeFilterShortIntvMidRForDYN(Ts, TeNoNoise, Te - oriEndTe - 1, hasE, newE, newENum, oriEndTe, lastTimeNoNoise, edgeSetsRNoNoise, selectedNumNoNoise, fixLabel, isEdgeTypeFixed);//not noise
				}
			}
			else {//same as FRTMOpt1
				if (Ts != startT) {
					root2Comp.clear();
					disjointSet->reset(vertexNum);
					hasChecked.clear();
					checkedCC->release();
				}
				edgeFilterOpt1MidRForDYN(Ts, Te, oriEndTe, edgeSetsRAdd,
					edgeSetsR, selectedNum, rightEndpoint, k, fixLabel, isEdgeTypeFixed);
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
				edgesAdd.clear();
				for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {
					if (tempT > rightEndpoint) continue;
					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();

					//maintain sets which need to add and delete O(delta |R|)
					for (auto add : edgeSetsRAdd[tempT]) {
						edgesAdd.emplace(add);
					}
					edgeSetsRAdd[tempT].clear();
					BEGIN_TIMER(c)
#pragma region maxCheck
					if (iterStart != iterEnd) {//new edges fetched
						maintainConnectivity(iterStart, iterEnd,
							vertex2Pos, disjointSet, vertexNum,
							combineMotifPos);
						combineComponentsOpt1(setCC,
							vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
							combineMotifPos);
					}
					updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,/**/
						setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
						saveMotifPos, Ts, edgeEndT, k, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion
					BEGIN_TIMER(d)
						expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
							edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1MidR);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();
				}
			}
			else {
			edgesAdd.clear();
			for (int tempT = endT - oriEndTe; tempT >= 0; tempT--, edgeEndT--) {
				if (edgeEndT >= Te) {//large intervals
					if (tempT > rightEndpoint) continue;
					iterStart = edgeSetsR[tempT].begin();
					iterEnd = edgeSetsR[tempT].end();

					//maintain sets which need to add and delete O(delta |R|)
					for (auto add : edgeSetsRAdd[tempT]) {
						edgesAdd.emplace(add);
					}
					edgeSetsRAdd[tempT].clear();
					BEGIN_TIMER(c)
#pragma region maxCheck
						if (iterStart != iterEnd) {//new edges fetched
							maintainConnectivity(iterStart, iterEnd,
								vertex2Pos, disjointSet, vertexNum,
								combineMotifPos);
							combineComponentsOpt1(setCC,
								vertex2Pos, disjointSet, root2Comp, checkedCC, hasChecked,
								/*setCCSize, realMotifNum,*/
								combineMotifPos);
						}
					updateNewEdgeInfoOpt1(iterStart, iterEnd, edgesAdd,/**/
						setCC, vertex2Pos, disjointSet, tempDisjointSet, root2Comp, &subCCMap, subCCId, root2Id, checkedCC, hasChecked, tempRecordFromQueue,
						saveMotifPos, Ts, edgeEndT, k, haveNewEdgeCC, remainEdges, &TGraph::maintainTempCCForUFOpt1);//O(¦¤Es)
					Test::gm += END_TIMER(c);
#pragma endregion
					BEGIN_TIMER(d)
						expCheckFRTM(saveMotifPos, setCC, newResult, k, Ts,
							edgeEndT, expandMask, /*expCheck,*/ motifNumber, &TGraph::checkExpandableOpt1MidR);
					Test::gne += END_TIMER(d);
					edgeSetsR[tempT].clear();
				}
				else {//short intervals

					if (Ts <= lastTimeNoNoise && edgeEndT == Te - 1) {//record root of disjointSet to be saved
						recordSavedRoot(setCC, hasChecked, newE, newENum, saveR, disjointSet, vertex2Pos);
					}
					iterStartNoNoise = edgeSetsRNoNoise[tempT].begin();
					iterEndNoNoise = edgeSetsRNoNoise[tempT].end();
					if (iterStartNoNoise == iterEndNoNoise) continue;
					BEGIN_TIMER(c)
#pragma region generateMaxTM
						addENum = 0;
					maintainConnectivityPlus(iterStartNoNoise, iterEndNoNoise, addE, addENum,
						vertex2Pos, disjointSetNoNoise, vertexNum, hasE, saveR, disjointSet, combineMotifPosNoNoise);
					combineComponentsPlus(setCCShortIntv,
						vertex2Pos, disjointSetNoNoise, root2CompNoNoise, combineMotifPosNoNoise);
					updateNewEdgeInfoPlus(iterStartNoNoise, addE, addENum,
						setCCShortIntv, vertex2Pos, disjointSetNoNoise,
						root2CompNoNoise, hasE, saveMotifPosNoNoise, Ts, edgeEndT);
					Test::gm += END_TIMER(c);
#pragma endregion

					BEGIN_TIMER(d)
						expCheckFRTMPlusMidR(saveMotifPosNoNoise, setCCShortIntv,
							newResult, k, Ts, edgeEndT, expandMask, /*expCheck,*/ motifNumber, disjointSet, hasChecked, root2Comp, vertex2Pos, setCC);
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
					auto motifEdges = newResult[resultP][resultIndex]->getMotifEdge();
					auto motifEdgesEnd = motifEdges->end();
					auto motifEdgesIter = motifEdges->begin();
					bool flag = false;
					for (; motifEdgesIter != motifEdgesEnd; ++motifEdgesIter) {
						id = *motifEdgesIter;
						int inMidP = id * numOfLabel + getEdgeLabel(id, motifEndT);
						//if ( (id == 1334 || id == 174)){
						//cout << Ts << " " << posInEIntR[inMidP] << " "<< maxIntv[id].first << " "<< maxIntv[id].second<<endl; }
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
			listEnd = setCC.end();
			for (veciter(CComponentsII*) listIter = setCC.begin();
				listIter != listEnd; ++listIter) {
				if (*listIter != nullptr) {
					CComponentsFRTM* listIt = (CComponentsFRTM*)*listIter;
					delete *listIter;
				}
			}
			setCC.clear();

			listEndNoNoise = setCCShortIntv.end();
			for (veciter(CComponentsShortIntv*) listIterNoNoise = setCCShortIntv.begin();
				listIterNoNoise != listEndNoNoise; ++listIterNoNoise) {
				if (*listIterNoNoise != nullptr)
					delete *listIterNoNoise;
			}
			setCCShortIntv.clear();
		}
		delete[] hasE;
		delete[] newE;
		delete[] addE;
		delete[] saveR;
		delete[] edgeSetsRNoNoise;
		delete disjointSetNoNoise;
		delete disjointSet;
		delete[] subCCId;
		delete tempDisjointSet;
		delete[] expandMask;
		delete[] edgeSetsR;
		//delete[] edgeSetsRAdd;
		delete[] remainEdges;
		delete checkedCC;
		cout << SELECT_EDGE << maxSelect << endl;
		cout << MEAN_SELECT_EDGE << selectedSum / (lastTime + 1) << endl;
		cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime + 1) << endl;
		cout << "maxSmk: " << maxSmk << " averSmk: " << allSmk / (lastTime + 1) << endl;
	}
	else {
		findRTMotifsDynamicOpt1(k, newResult, oriEndT, motifNumber, fixLabel, isEdgeTypeFixed);
	}
}


void TGraph::getIntersectionOfScope(vec(int)& edges, pair<int, int>& intersectIntv) {
	veciter(int) edgeEnd = edges.end();
	int edgeId;
	auto edgeIter = edges.begin();
	edgeId = *edgeIter;
	intersectIntv.first = scope[edgeId].first;
	intersectIntv.second = scope[edgeId].second;
	++edgeIter;
	for (; edgeIter != edgeEnd; ++edgeIter) {
		edgeId = *edgeIter;

		intersectIntv.first = max(intersectIntv.first, scope[edgeId].first);
		intersectIntv.second = min(intersectIntv.second, scope[edgeId].second);
	}
}

void TGraph::getIntersectionOfScope(vec(int)& edges, int*&subCCId, pair<int, int>& intersectIntv) {
	veciter(int) edgeEnd = edges.end();
	int edgeId;
	intersectIntv.first = 0;
	intersectIntv.second = 0x7fffffff;
	auto subCCIter = &subCCId[0];
	for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter,++subCCIter) {
		if (*subCCIter != -1) {
			edgeId = *edgeIter;
			intersectIntv.first = max(intersectIntv.first, scope[edgeId].first);
			intersectIntv.second = min(intersectIntv.second, scope[edgeId].second);
		}
	}
}

void TGraph::getIntersectionOfScope(vec(int)& subCCs, vec(int)& edges, pair<int, int>& intersectIntv) {
	auto edgeEnd = edges.end();
	auto subCCEnd = subCCs.end(); 
	int edgeId;
	auto subCCIter = subCCs.begin();
	edgeId = edges[*subCCIter];
	intersectIntv.first = scope[edgeId].first;
	intersectIntv.second = scope[edgeId].second;
	++subCCIter;
	for (; subCCIter != subCCEnd; ++subCCIter) {
		edgeId = edges[*subCCIter];
		intersectIntv.first = max(intersectIntv.first, scope[edgeId].first);
		intersectIntv.second = min(intersectIntv.second, scope[edgeId].second);
	}
}
