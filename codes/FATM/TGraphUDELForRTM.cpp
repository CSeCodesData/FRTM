#include "TGraphUDEL.h"


#pragma region FRTMPlus For intvB > lastTimeShortIntv
void TGraphUDEL::edgeFilterShortIntv(int intvB, int intvE, int limited, 
	vec(int)*& edgeSetsR, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed) {
	int timePos = intvE - startT;
	int check;
	int edgeType;
	int intvLen = intvE - intvB + 1;
	int endPos = endT - startT;
	int intvStart;
	for (int i = 0; i < nEdge; i++) {

		edgeType = lab[TGraph::posUsedForEdgeFilterShortIntv + i];//label at intvE
		if (isEdgeTypeFixed && !fixLabel[lab[TGraph::posUsedForEdgeFilter + i]])  continue;//label at intvB
		/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
			put here for less time*/
		if (max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) < intvLen) {
			continue;
		}
		intvStart = intvE - max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) + 1;
		int j = timePos + 1;
		if (maxIntvShortIntv[i].second >= timePos && maxIntvShortIntv[i].first <= intvB - startT) {//case 4
			check = min(maxIntvShortIntv[i].second - timePos, limited);
			selectedNum++;
			edgeSetsR[check].emplace_back(i);
		}
		else if (maxIntvShortIntv[i].second >= intvB - startT && maxIntvShortIntv[i].first <= timePos) {//case 3
			continue;
		}
		else {//case 1 and 2
			maxIntvShortIntv[i].first = intvStart;
			selectedNum++;

			check = max(aft[TGraph::posUsedForEdgeFilterShortIntv + i], 1) - 1;
			maxIntvShortIntv[i].second = timePos + check + startT;

			check = min(check, limited);
			edgeSetsR[check].emplace_back(i);

		}
	}
}
void TGraphUDEL::edgeFilterShortIntvMidR(int intvB, int intvE, int limited, int lastTimeShortIntv,
	vec(int)*& edgeSetsR, int& selectedNum, int choiceEndT, bool*& fixLabel, bool isEdgeTypeFixed) {
	int timePos = intvE - startT;
	int check;
	int edgeType;
	int intvLen = intvE - intvB + 1;
	int endPos = endT - startT;
	int intvStart;
	int labelPosForEdge = 0;
	int mainLabelPos;
	for (int i = 0; i < nEdge; labelPosForEdge += numOfLabel, i++) {
		edgeType = lab[TGraph::posUsedForEdgeFilterShortIntv + i]; //label at intvE
		if (isEdgeTypeFixed && !fixLabel[lab[TGraph::posUsedForEdgeFilter + i]])  continue;//label at intvB
	
		//dynamic
		int pos = intvB * nEdge + i;
		int mainLabel = lab[pos];
		mainLabelPos = labelPosForEdge + mainLabel;
		int checkE = scanT[mainLabelPos];//scan at the previous row
		if (checkE == choiceEndT) {
			int localNoise = 0;
			pos = checkE * nEdge + i;
			if (lab[pos] != mainLabel) {
				int localNoise = 0;
				int gap;
				while (lab[pos]!=mainLabel) {
					gap = max(bef[pos], 1);
					localNoise += gap;
					pos -= gap * nEdge;
				}
				
				if (localNoise <= Setting::c) {
					if (newPosInEIntR[mainLabelPos] < 0) {
						auto& now = vioT[mainLabelPos];

						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i, maxIntvShortIntv[i], preMaxIntv[i], checkE, now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] < 0) {
					auto& now = vioT[mainLabelPos];

					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i, maxIntvShortIntv[i], preMaxIntv[i], checkE, now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
		}
		else if (checkE < intvB) {//from scratch
			auto& now = vioT[mainLabelPos]; 
			now.removeNodeLessThan(intvB);

			pos = intvB * nEdge + i;
			int currentPos = intvB - startT;
			lazyUpdate(currentPos, pos, i);//update aft
			int nextLabPos = min(currentPos + max(aft[pos], 1) - 1, endT - startT) + 1;
			int localNoise;
			int labelsSum = 0, noiseNum = 0, forbidTimeStartT;
			while (currentPos < currNTimestamp) {
				int intvL = nextLabPos - currentPos;
				edgeType = lab[pos];
				if (edgeType == mainLabel) {
					localNoise = 0;
					double minNoise = Setting::delta * (labelsSum + 1);
					if (LESSEQ(noiseNum, minNoise)) {
						if (forbidTimeStartT != -1) {
							now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
							forbidTimeStartT = -1;
						}
					}
					else if (MORE(noiseNum, Setting::delta * (labelsSum + intvL))) {
						if (forbidTimeStartT == -1) {
							forbidTimeStartT = max(endPos, currentPos) + startT;
						}
					}
					else {
						double checkNum = (noiseNum - minNoise) / Setting::delta;
						int noiseT;
						if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
						else noiseT = (int)checkNum + currentPos + startT;
						if (forbidTimeStartT == -1) {
							now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
						}
						else {
							now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
							forbidTimeStartT = -1;
						}
					}
				}
				else {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
					localNoise += intvL;
					noiseNum += intvL;
					if (localNoise > Setting::c) break;
					//if (MORE(noiseNum, Setting::delta * allLen)) break;
				}
				labelsSum += intvL;

				currentPos = nextLabPos;
				pos = currentPos * nEdge + i;
				if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
					if (currentPos < currNTimestamp){
						lazyUpdate(currentPos - 1, pos - nEdge, i);//update aft
						nextLabPos = min(nextLabPos - 2 - aft[pos - nEdge], endT - startT) + 1;
					}
				}
				else {
					if (currentPos < currNTimestamp) {
						lazyUpdate(currentPos, pos, i);//update aft
						nextLabPos = min(currentPos + max(aft[pos], 1) - 1, endT - startT) + 1;
					}
				}
			}
			if (forbidTimeStartT != -1) {
				now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
			}
			scanT[mainLabelPos] = nextLabPos - 1;

			if (localNoise <= Setting::c) {
				if (newPosInEIntR[mainLabelPos] < 0) {
					auto& now = vioT[mainLabelPos];

					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i, maxIntvShortIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else if (newPosInEIntR[mainLabelPos] >= 0) {
				newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
				newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
			}
		}

		/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
			put here for less time*/
		if (max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) < intvLen) {
			continue;
		}
		intvStart = intvE - max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) + 1;
		int j = timePos + 1;
		if (maxIntvShortIntv[i].second >= timePos && maxIntvShortIntv[i].first <= intvB - startT) {//case 4
			check = min(maxIntvShortIntv[i].second - timePos, limited);
			selectedNum++;
			edgeSetsR[check].emplace_back(i);
		}
		else if (maxIntvShortIntv[i].second >= intvB - startT && maxIntvShortIntv[i].first <= timePos) {//case 3
			continue;
		}
		else {//case 1 and 2
			maxIntvShortIntv[i].first = intvStart;
			selectedNum++;
			lazyUpdate(intvE, TGraph::posUsedForEdgeFilterShortIntv + i, i);//update aft
			check = max(aft[TGraph::posUsedForEdgeFilterShortIntv + i], 1) - 1;

			maxIntvShortIntv[i].second = timePos + check + startT;
			check = min(check, limited);
			edgeSetsR[check].emplace_back(i); 
		}
	}
}

void TGraphUDEL::edgeFilterShortIntvMidRForDYN(int intvB, int intvE, int limited, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, int oriEndTe, int lastTimeNoNoise,
	vec(int)*& edgeSetsR, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed) {
	int timePos = intvE - startT;
	int check;
	int edgeType;
	int intvLen = intvE - intvB + 1;
	int endPos = endT - startT;
	int intvStart;
	int labelPosForEdge = 0;
	int mainLabelPos;
	auto fromMidREnd = edgesInEIntR->end(), fromMidRIter = edgesInEIntR->begin();
	for (; fromMidRIter != fromMidREnd; ++fromMidRIter) {
		int i = *fromMidRIter;
		labelPosForEdge = i * numOfLabel;

		int pos = TGraph::posUsedForEdgeFilter + i;
		int mainLabel = lab[pos];
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;

		mainLabelPos = labelPosForEdge + mainLabel;
		if (posInEIntR[mainLabelPos] == EMaxIntvlChange::INTVINIT) continue;//initial 

		//dynamic
		int checkE = scanT[mainLabelPos];//scan at the previous row
		if (checkE == endT) {
			int localNoise = 0;
			pos = checkE * nEdge + i;
			if (lab[pos] != mainLabel) {
				int localNoise = 0;
				int gap;
				while (lab[pos] != mainLabel) {
					gap = max(bef[pos], 1);
					localNoise += gap;
					pos -= gap * nEdge;
				}
				if (localNoise <= Setting::c) {
					if (newPosInEIntR[mainLabelPos] < 0) {
						auto& now = vioT[mainLabelPos];

						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i, maxIntvShortIntv[i], preMaxIntv[i], checkE, now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] < 0) {
					auto& now = vioT[mainLabelPos];

					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i, maxIntvShortIntv[i], preMaxIntv[i], checkE, now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
		}

		edgeType = lab[TGraph::posUsedForEdgeFilterShortIntv + i];
		/*not exists the maximum interval containing [intvB,intvE]*/
		if (max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) < intvLen) {
			if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
			continue;
		}
		intvStart = intvE - max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) + 1;
		int j = timePos + 1;
		if (maxIntvShortIntv[i].second >= timePos && maxIntvShortIntv[i].first <= intvB - startT) {//case 4
			if (maxIntvShortIntv[i].second >= oriEndTe) {
				check = min(maxIntvShortIntv[i].second - oriEndTe, limited);
				selectedNum++;
				edgeSetsR[check].emplace_back(i);

				posInEIntR[mainLabelPos] = -1;//R for edge i is changed
			}
		}
		else if (maxIntvShortIntv[i].second >= intvB - startT && maxIntvShortIntv[i].first <= timePos) {//case 3
			if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
			continue;
		}
		else {//case 1 and 2
			maxIntvShortIntv[i].first = intvStart;
			selectedNum++;
			lazyUpdate(intvE, TGraph::posUsedForEdgeFilterShortIntv + i, i);//update aft
			check = max(aft[TGraph::posUsedForEdgeFilterShortIntv + i], 1) - 1;
			maxIntvShortIntv[i].second = timePos + check + startT;

			if (maxIntvShortIntv[i].second >= oriEndTe) {
				check = min(maxIntvShortIntv[i].second - oriEndTe, limited);
				edgeSetsR[check].emplace_back(i);

				posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R for edge i is changed
			}
			else if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R for edge i is unchanged 
		}
	}
}
#pragma endregion

#pragma region FRTM

void TGraphUDEL::edgeFilterFRTM(int intvB, int intvE,
	vec(int)*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;

	int timestampPosForEMaxIntvlEndT;
	for (int i = 0; i < nEdge; i++, labelPosForEdge += numOfLabel) {//O(|E|)

		mainLabel = lab[posUsedForEdgeFilter + i];
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		mainLabelPos = labelPosForEdge + mainLabel;

		auto& now = vioT[mainLabelPos];
		currentPos = scanT[mainLabelPos];
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		
		if (currentPos != 0) {
			if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {
				if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid intervals
					continue;
				}
				else {
					auto intvItem = now.first;

					bool shrink = false;
					while (intvItem != nullptr) {

						if (intvItem->item.first > maxIntv[i].second) break;
						else if (intvItem->item.second >= maxIntv[i].second) {//shrink
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}
						
						int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
						if (lab[tempEWPos] == mainLabel) {
							intvItem->item.second++;
						}
						else {
							int labelSum = intvItem->item.second - intvB + 2;
							if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
								intvItem->item.second++;
							}
						}

						if (intvItem->next != nullptr) {
							if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
								intvItem->item.second = intvItem->next->item.second;
								now.deleteNextNode(intvItem);
							}
							else intvItem = intvItem->next;
						}
						else if (intvItem->item.second >= maxIntv[i].second) {//shrink
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp); 
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}
						else intvItem = intvItem->next;
					}

					if (!shrink) {
						savePos = maxIntv[i].second - intvE;
						if (savePos >= 0) {
							for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
								temp++) {
								auto& intv = preEMaxIntvlPtr.q[temp];
								if (intv.second < intvE) {//no overlap
									preEMaxIntvlPtr.swapToTop(temp);
									preEMaxIntvlPtr.pop();
								}
								//else break;
							}
							selectedNum++;
							if (rightEndpoint < savePos) rightEndpoint = savePos;
							edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);

						}
					}
				}

				//Test::comprp1n++;
				//Test::comprp1t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
				now.removeNodeLessThan(intvE);
				continue;
			}
			else {
				//Test::comprp2n++;
				now.removeNodeLessThan(intvE);

				if (maxIntv[i].second == -1) {
					EMaxIntvlEndT = -1;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					EMaxIntvlEndT = maxIntv[i].second;
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
				}
				lastMainLabelT = max(EMaxIntvlEndT, intvB);
				if (currentPos < beginPos) { //does not scan
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					noiseNum = 0;
				}
				else {
					auto forbidIntv = now.first;
					if (forbidIntv == nullptr) {//does not need update
						if (currentPos != currNTimestamp) {
							int timestampPosForEL = currentPos * nEdge + i;
							edgeType = lab[timestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[timestampPosForEL], 1);
								int tempTPos = timestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = currentPos - eLvalue;
								localNoise = eLvalue;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								localNoise = 0;
								noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
							}
						}
					}
					else {
						//update vioT

						//get noiseNum, localNoise
						int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
						int tempPos = tempT - startT, intvLen = tempT - intvB;
						labelsSum = tempT - intvB;
						int tempTimestampPosForEL = tempPos * nEdge + i;
						edgeType = lab[tempTimestampPosForEL];
						if (edgeType != mainLabel) {
							int eLvalue = max(bef[tempTimestampPosForEL], 1);
							int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
							int tempPosForLabels = tempPos - eLvalue;
							localNoise = eLvalue - 1;
							while (lab[tempTPos] != mainLabel) {
								eLvalue = max(bef[tempTPos], 1);
								tempPosForLabels -= eLvalue;
								localNoise += eLvalue;
								tempTPos -= eLvalue * nEdge;
							}
							noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
						}
						else {
							noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
							localNoise = 0;
						}

						auto tempIntv = forbidIntv;
						while (forbidIntv != nullptr) {

							forbidTimeStartT = -1;
							int nextLabT = min(tempT + max(aft[tempTimestampPosForEL],1)-1, endT) + 1;
							while (tempT <= stopT) {

								intvLen = nextLabT - tempT;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType == mainLabel) {
									localNoise = 0;
									minNoise = Setting::delta * (labelsSum + 1);
									if (LESSEQ(noiseNum, minNoise)) {
										//lastMainLabelPos = nextLabPos - 1;
										//remove = 0;
										if (forbidTimeStartT != -1) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											forbidTimeStartT = -1;
											tempIntv = tempIntv->next;
										}
									}
									else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
										//remove += intvLen;
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
									}
									else {
										//remove = 0;
										//lastMainLabelPos = nextLabPos - 1;
										checkNum = (noiseNum - minNoise) / Setting::delta;
										if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
										else noiseT = (int)checkNum + tempT;
										if (forbidTimeStartT == -1) {
											now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT,stopT)));
											tempIntv = tempIntv->next;
										}
										else {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
											forbidTimeStartT = -1;
										}
									}
								}
								else {
									if (forbidTimeStartT == -1) {
										forbidTimeStartT = tempT + startT;
									}
									localNoise += intvLen;
									noiseNum += intvLen;
								}
								labelsSum += intvLen;

								tempT = nextLabT;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
									if (tempT <= stopT)
										nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
								}
								else {
									if (tempT <= stopT)
										nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
								}
							}
							if (forbidTimeStartT != -1) {
								//now.vioT.addNodeAft(tempIntv, make_pair(forbidTimeStartT,
								if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
									now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
									tempIntv = tempIntv->next;

									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;
								}
								else {
									tempIntv = forbidIntv = forbidIntv->next;
								}
							}
							else {
								forbidIntv = forbidIntv->pre;
								tempIntv = tempIntv->next;
								if (forbidIntv == nullptr) {
									now.deleteFirstNode();
								}
								else now.deleteNextNode(forbidIntv);
								forbidIntv = tempIntv;
							}

							if (forbidIntv != nullptr) {
								tempT = forbidIntv->item.first;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								labelsSum = tempT - intvB;
								localNoise = 0;
								stopT = forbidIntv->item.second;
							}
						}

						//get the last time where the edge has main label
						auto lastIntv = now.tail;
						int updateEndT = min(currentPos + startT, endT);
						if (lastIntv != nullptr) {
							if (lastIntv->item.second != updateEndT) {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
							else {
								lastMainLabelT = lastIntv->item.first - 1;
							}
						}
						else {
							lastMainLabelT = updateEndT;
							localNoise = 0;
						}
					}
					if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
						if (lastMainLabelT >= intvE) {
							if (lastMainLabelT == maxIntv[i].second) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							else if (lastMainLabelT < maxIntv[i].second) {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}

								if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
									if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							else {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}
								if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
									preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
								}

								if (EMaxIntvlEndT < lastMainLabelT) {
									maxIntv[i].first = intvB;
									maxIntv[i].second = lastMainLabelT;
								}
								else {
									maxIntv[i].first = EMaxIntvlStartT;
									maxIntv[i].second = EMaxIntvlEndT;
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
						}

						continue;
					}
					else {
						currentPos = max(currentPos + 1, beginPos);
						lastMainLabelPos = lastMainLabelT - startT;
					}
				}


				//Test::comprp2t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
			}

			forbidTimeStartT = -1;
			//auto tail = now.vioT.tail;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second+1 == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
		}
		else {
			localNoise = 0;
			currentPos = lastMainLabelPos = beginPos;
			forbidTimeStartT = -1;
			noiseNum = 0;
		}
		labelsSum = currentPos - beginPos;
		
		int tempTimestampPosForIT = currentPos * nEdge + i;
		
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT],1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {
				if (currentPos < currNTimestamp)
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
			}
			else {
				if (currentPos < currNTimestamp)
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;
		
		
		if (lastMainLabelPos < endPos) {
			continue;
		}

		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			selectedNum++;
			savePos = lastMainLabelPos - endPos;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);

		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
		}
		else {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
		}
	}
}

void TGraphUDEL::edgeFilterFRTMMidR(int intvB, int intvE,
	vec(int)*& edgeSetsR, int& selectedNum, int& rightEndpoint,
	int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;

	int timestampPosForEMaxIntvlEndT;
	for (int i = 0; i < nEdge; i++, labelPosForEdge += numOfLabel) {//O(|E|)
		mainLabel = lab[posUsedForEdgeFilter + i];
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		mainLabelPos = labelPosForEdge + mainLabel;

		auto& now = vioT[mainLabelPos];
		currentPos = scanT[mainLabelPos];
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];

		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		
		if (currentPos != 0) {
			
			if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {
				if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid intervals
					continue;
				}
				else {
					auto intvItem = now.first;

					bool shrink = false;
					while (intvItem != nullptr) {

						if (intvItem->item.first > maxIntv[i].second) break;
						else if (intvItem->item.second >= maxIntv[i].second) {//shrink
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preMaxIntv[i].swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}

						int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
						if (lab[tempEWPos] == mainLabel) {
							intvItem->item.second++;
						}
						else {
							int labelSum = intvItem->item.second - intvB + 2;
							if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
								intvItem->item.second++;
							}
						}

						if (intvItem->next != nullptr) {
							if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
								intvItem->item.second = intvItem->next->item.second;
								now.deleteNextNode(intvItem);
							}
							else intvItem = intvItem->next;
						}
						else if (intvItem->item.second >= maxIntv[i].second) {//shrink
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}
						else intvItem = intvItem->next;
					}

					if (!shrink) {
						savePos = maxIntv[i].second - intvE;
						if (savePos >= 0) {
							for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
								temp++) {
								auto& intv = preEMaxIntvlPtr.q[temp];
								if (intv.second < intvE) {//no overlap
									preMaxIntv[i].swapToTop(temp);
									preEMaxIntvlPtr.pop();
								}
								//else break;
							}
							selectedNum++;
							if (rightEndpoint < savePos) rightEndpoint = savePos;
							edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);

						}
					}
				}

				now.removeNodeLessThan(intvE);
				continue;
			}
			else {
				//Test::comprp2n++;
				//now.vioT.removeNodeLessThan(intvE);
				now.removeNodeLessThan(intvE);

				if (maxIntv[i].second == -1) {
					EMaxIntvlEndT = -1;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					EMaxIntvlEndT = maxIntv[i].second;
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
				}
				lastMainLabelT = max(EMaxIntvlEndT, intvB);
				if (currentPos < beginPos) { //does not scan
					//now.localNoise = 0;
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					noiseNum = 0;
				}
				else {
					auto forbidIntv = now.first;
					if (forbidIntv == nullptr) {//does not need update
						if (currentPos != currNTimestamp) {
							int timestampPosForEL = currentPos * nEdge + i;
							edgeType = lab[timestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[timestampPosForEL], 1);
								int tempTPos = timestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = currentPos - eLvalue;
								localNoise = eLvalue;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
								localNoise = 0;
							}
						}
					}
					else {
						//update vioT

						//get noiseNum, localNoise
						int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
						int tempPos = tempT - startT, intvLen = tempT - intvB;
						labelsSum = tempT - intvB;
						int tempTimestampPosForEL = tempPos * nEdge + i;
						edgeType = lab[tempTimestampPosForEL];
						if (edgeType != mainLabel) {
							int eLvalue = max(bef[tempTimestampPosForEL], 1);
							int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
							int tempPosForLabels = tempPos - eLvalue;
							localNoise = eLvalue - 1;
							while (lab[tempTPos] != mainLabel) {
								eLvalue = max(bef[tempTPos], 1);
								tempPosForLabels -= eLvalue;
								localNoise += eLvalue;
								tempTPos -= eLvalue * nEdge;
							}
							noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
						}
						else {
							noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
							localNoise = 0;
						}

						auto tempIntv = forbidIntv;
						while (forbidIntv != nullptr) {

							forbidTimeStartT = -1;
							lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
							int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
							while (tempT <= stopT) {

								intvLen = nextLabT - tempT;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType == mainLabel) {
									localNoise = 0;
									minNoise = Setting::delta * (labelsSum + 1);
									if (LESSEQ(noiseNum, minNoise)) {
										//lastMainLabelPos = nextLabPos - 1;
										//remove = 0;
										if (forbidTimeStartT != -1) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											forbidTimeStartT = -1;
											tempIntv = tempIntv->next;
										}
									}
									else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
										//remove += intvLen;
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
									}
									else {
										//remove = 0;
										//lastMainLabelPos = nextLabPos - 1;
										checkNum = (noiseNum - minNoise) / Setting::delta;
										if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
										else noiseT = (int)checkNum + tempT;
										if (forbidTimeStartT == -1) {
											now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
										}
										else {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
											forbidTimeStartT = -1;
										}
									}
								}
								else {
									if (forbidTimeStartT == -1) {
										forbidTimeStartT = tempT + startT;
									}
									localNoise += intvLen;
									noiseNum += intvLen;
								}
								labelsSum += intvLen;

								tempT = nextLabT;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
									if (tempT <= stopT) {
										lazyUpdate(tempT - 1, tempTimestampPosForEL - nEdge, i);//update aft
										nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
									}
								}
								else {
									if (tempT <= stopT) {
										lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
										nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
									}
								}
							}
							if (forbidTimeStartT != -1) {
								if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
									now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
									tempIntv = tempIntv->next;

									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;

								}
								else {
									tempIntv = forbidIntv = forbidIntv->next;
								}
							}
							else {
								forbidIntv = forbidIntv->pre;
								tempIntv = tempIntv->next;
								if (forbidIntv == nullptr) {
									now.deleteFirstNode();
								}
								else now.deleteNextNode(forbidIntv);
								forbidIntv = tempIntv;
							}

							if (forbidIntv != nullptr) {
								tempT = forbidIntv->item.first;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								labelsSum = tempT - intvB;
								localNoise = 0;
								stopT = forbidIntv->item.second;
							}
						}

						//get the last time where the edge has main label
						auto lastIntv = now.tail;
						int updateEndT = min(currentPos + startT, endT);
						if (lastIntv != nullptr) {
							//tie(intvStartT, intvEndT/*, std::ignore*/) = lastIntv->item;
							if (lastIntv->item.second != updateEndT) {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
							else {
								lastMainLabelT = lastIntv->item.first - 1;
							}
						}
						else {
							lastMainLabelT = updateEndT;
							localNoise = 0;
						}
					}
					if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
						int checkPos = (lastMainLabelT - startT)*nEdge + i;
						lazyUpdate(lastMainLabelT, checkPos, i);//update aft
						int labelEnd = lastMainLabelT - startT + max(aft[checkPos], 1) - 1;
						//dynamic
						if (localNoise <= Setting::c) {
							checkPos = labelEnd * nEdge + i;
							lazyUpdate(labelEnd, checkPos, i);//update aft
							if (aft[checkPos] != -MYINFINITE) {
								if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
									newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
									MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
									newEIntR->emplace_back(rd4);
									newValidMidResult.emplace_back(true);
								}
							}
							else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
								newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
								MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
								newEIntR->emplace_back(rd4);
								newValidMidResult.emplace_back(true);
							}
						}
						else {
							if (newPosInEIntR[mainLabelPos] >= 0) {
								newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
								newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
							}
						}

						if (lastMainLabelT >= intvE) {
							if (lastMainLabelT == maxIntv[i].second) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							else if (lastMainLabelT < maxIntv[i].second) {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									//if (intv.endT < intvE) {//no overlap
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}

								if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
									if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							else {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									//if (intv.endT < intvE) {//no overlap
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;

										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}
								if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
									preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
								}

								if (EMaxIntvlEndT < lastMainLabelT) {
									maxIntv[i].first = intvB;
									maxIntv[i].second = lastMainLabelT;
								}
								else {
									maxIntv[i].first = EMaxIntvlStartT;
									maxIntv[i].second = EMaxIntvlEndT;
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
						}

						continue;
					}
					else {
						currentPos = max(currentPos + 1, beginPos);
						lastMainLabelPos = lastMainLabelT - startT;
					}
				}
				//Test::comprp2t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
			}

			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
		}
		else {
			localNoise = 0;
			currentPos = lastMainLabelPos = beginPos;
			forbidTimeStartT = -1;
			noiseNum = 0;
		}
		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;

		lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos - 1, tempTimestampPosForIT - nEdge, i);//update aft
					nextLabPos = min(nextLabPos - aft[tempTimestampPosForIT - nEdge] - 2, endT - startT) + 1;
				}
			}
			else {
				if (currentPos < currNTimestamp){
					lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
				}
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;


		if (lastMainLabelPos < endPos) {
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			
			continue;
		}
		//vioT[i].removeNodeLessThan(intvE);

		//if ((isEdgeTypeFixed && fixLabel.find(idToLabel[mainLabel[i]]) == fixLabelEnd)) continue;
		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			selectedNum++;
			savePos = lastMainLabelPos - endPos;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);

		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
		}
		else {

			//dynamic
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			

			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}


			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
		}
	}
}

void TGraphUDEL::edgeFilterFRTMMidRForDYN(int intvB, int intvE, int oriEndTE,
	vec(int)*& edgeSetsR, int& selectedNum, int& rightEndpoint, 
	int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT, checkEndPos = oriEndTE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;

	int timestampPosForEMaxIntvlEndT;
	auto fromMidREnd = edgesInEIntR->end(), fromMidRIter = edgesInEIntR->begin();
	for (; fromMidRIter != fromMidREnd; ++fromMidRIter) {

		int i = *fromMidRIter;
		labelPosForEdge = i * numOfLabel;
		//auto beginTest = std::chrono::steady_clock::now();

		mainLabel = lab[posUsedForEdgeFilter + i];
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		mainLabelPos = labelPosForEdge + mainLabel;

		auto& now = vioT[mainLabelPos];
		currentPos = scanT[mainLabelPos];
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];
		
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		if (posInEIntR[mainLabelPos] == EMaxIntvlChange::INTVINIT) {//initial 
			continue;
		}
		else if (posInEIntR[mainLabelPos] >= 0) {//continue to scan
			if (currentPos == 0) continue; //initial

			//recover localNoise , noiseNum and lastMainLabelPos
			int timestampPosForEL = currentPos * nEdge + i;
			edgeType = lab[timestampPosForEL];
			if (edgeType != mainLabel) {
				int eLvalue = max(bef[timestampPosForEL], 1);
				int tempTPos = timestampPosForEL - eLvalue * nEdge;
				int tempPosForLabels = currentPos - eLvalue;
				localNoise = eLvalue;
				while (lab[tempTPos] != mainLabel) {
					eLvalue = max(bef[tempTPos], 1);
					tempPosForLabels -= eLvalue;
					localNoise += eLvalue;
					tempTPos -= eLvalue * nEdge;
				}
				noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
			}
			else {
				localNoise = 0;
				noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
			}

			if (maxIntv[i].second == -1) {
				lastMainLabelPos = beginPos;
			}
			else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
				lastMainLabelPos = max(maxIntv[i].second, beginPos);
			}
			else {
				EMaxIntvlEndT = -1;
				int last = preEMaxIntvlPtr.rear;
				if (last != preEMaxIntvlPtr.front) {
					int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
					for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
						EMaxIntvlEndT = preEMaxIntvlPtr.q[temp].second;
						if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
					}
					if (temp == preEMaxIntvlPtr.front) {
						EMaxIntvlEndT = preEMaxIntvlPtr.q[temp].second;
						if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
					}
				}
				lastMainLabelPos = max(EMaxIntvlEndT, beginPos);
			}

			//recover vioT
			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				lastMainLabelPos = forbidTimeStartT - 1;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
			currentPos = max(currentPos + 1, beginPos);
		}
		else {
			if (currentPos != 0) {
				
				if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {
					if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid intervals
						continue;
					}
					else {
						auto intvItem = now.first;

						bool shrink = false;
						while (intvItem != nullptr) {

							if (intvItem->item.first > maxIntv[i].second) break;
							else if (intvItem->item.second >= maxIntv[i].second) {//shrink
								shrink = true;
								int tempEndT = intvItem->item.first - 1;
								if (tempEndT >= intvE) {
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
										temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preMaxIntv[i].swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										//else break;
									}
									savePos = tempEndT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
								break;
							}

							int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
							if (lab[tempEWPos] == mainLabel) {
								intvItem->item.second++;
							}
							else {
								int labelSum = intvItem->item.second - intvB + 2;
								if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
									intvItem->item.second++;
								}
							}

							if (intvItem->next != nullptr) {
								if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
									intvItem->item.second = intvItem->next->item.second;
									now.deleteNextNode(intvItem);
								}
								else intvItem = intvItem->next;
							}
							else if (intvItem->item.second >= maxIntv[i].second) {//shrink
								shrink = true;
								int tempEndT = intvItem->item.first - 1;
								if (tempEndT >= intvE) {
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
										temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										//else break;
									}
									savePos = tempEndT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
								break;
							}
							else intvItem = intvItem->next;
						}

						if (!shrink) {
							if (maxIntv[i].second >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}

								savePos = maxIntv[i].second - oriEndTE;
								if (savePos >= 0) {
									selectedNum++;
									if (rightEndpoint < savePos) rightEndpoint = savePos;
									edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								}
							}
						}
					}

					//Test::comprp1n++;
					//Test::comprp1t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
					now.removeNodeLessThan(intvE);
					continue;
				}
				else {
					//Test::comprp2n++;
					now.removeNodeLessThan(intvE);

					if (maxIntv[i].second == -1) {
						EMaxIntvlEndT = -1;
					}
					else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
						EMaxIntvlEndT = maxIntv[i].second;
					}
					else {
						EMaxIntvlEndT = -1;
						int last = preEMaxIntvlPtr.rear;
						if (last != preEMaxIntvlPtr.front) {
							int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
							for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
								tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
								if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
							}
							if (temp == preEMaxIntvlPtr.front) {
								tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
								if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
							}
						}
					}
					lastMainLabelT = max(EMaxIntvlEndT, intvB);
					if (currentPos < beginPos) { //does not scan
						localNoise = 0;
						currentPos = lastMainLabelPos = beginPos;
						noiseNum = 0;
					}
					else {
						auto forbidIntv = now.first;
						if (forbidIntv == nullptr) {//does not need update
							if (currentPos != currNTimestamp) {
								int timestampPosForEL = currentPos * nEdge + i;
								edgeType = lab[timestampPosForEL];
								if (edgeType != mainLabel) {
									int eLvalue = max(bef[timestampPosForEL], 1);
									int tempTPos = timestampPosForEL - eLvalue * nEdge;
									int tempPosForLabels = currentPos - eLvalue;
									localNoise = eLvalue;
									while (lab[tempTPos] != mainLabel) {
										eLvalue = max(bef[tempTPos], 1);
										tempPosForLabels -= eLvalue;
										localNoise += eLvalue;
										tempTPos -= eLvalue * nEdge;
									}
									noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
								}
								else {
									localNoise = 0;
									noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
								}
							}
						}
						else {
							//update vioT

							//get noiseNum, localNoise
							int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
							int tempPos = tempT - startT, intvLen = tempT - intvB;
							labelsSum = tempT - intvB;
							int tempTimestampPosForEL = tempPos * nEdge + i;
							edgeType = lab[tempTimestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[tempTimestampPosForEL], 1);
								int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = tempPos - eLvalue;
								localNoise = eLvalue - 1;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
								localNoise = 0;
							}

							auto tempIntv = forbidIntv;
							while (forbidIntv != nullptr) {

								forbidTimeStartT = -1;
								lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
								int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
								while (tempT <= stopT) {

									intvLen = nextLabT - tempT;
									edgeType = lab[tempTimestampPosForEL];
									if (edgeType == mainLabel) {
										localNoise = 0;
										minNoise = Setting::delta * (labelsSum + 1);
										if (LESSEQ(noiseNum, minNoise)) {
											//lastMainLabelPos = nextLabPos - 1;
											//remove = 0;
											if (forbidTimeStartT != -1) {
												now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
												forbidTimeStartT = -1;
												tempIntv = tempIntv->next;
											}
										}
										else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
											//remove += intvLen;
											if (forbidTimeStartT == -1) {
												forbidTimeStartT = tempT + startT;
											}
										}
										else {
											//remove = 0;
											//lastMainLabelPos = nextLabPos - 1;
											checkNum = (noiseNum - minNoise) / Setting::delta;
											if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
											else noiseT = (int)checkNum + tempT;
											if (forbidTimeStartT == -1) {
												now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
												tempIntv = tempIntv->next;
											}
											else {
												now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
												tempIntv = tempIntv->next;
												forbidTimeStartT = -1;
											}
										}
									}
									else {
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
										localNoise += intvLen;
										noiseNum += intvLen;
									}
									labelsSum += intvLen;

									tempT = nextLabT;
									tempTimestampPosForEL = (tempT - startT) * nEdge + i;
									if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
										if (tempT <= stopT) {
											lazyUpdate(tempT - 1, tempTimestampPosForEL - nEdge, i);//update aft
											nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
										}
									}
									else {
										if (tempT <= stopT){
											lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
											nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
										}
									}
								}
								if (forbidTimeStartT != -1) {
									if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
										now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
										tempIntv = tempIntv->next;

										forbidIntv = forbidIntv->pre;
										tempIntv = tempIntv->next;
										if (forbidIntv == nullptr) {
											now.deleteFirstNode();
										}
										else now.deleteNextNode(forbidIntv);
										forbidIntv = tempIntv;

									}
									else {
										tempIntv = forbidIntv = forbidIntv->next;
									}
								}
								else {
									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;
								}

								if (forbidIntv != nullptr) {
									tempT = forbidIntv->item.first;
									tempTimestampPosForEL = (tempT - startT) * nEdge + i;
									labelsSum = tempT - intvB;
									localNoise = 0;
									stopT = forbidIntv->item.second;
								}
							}

							//get the last time where the edge has main label
							auto lastIntv = now.tail;
							int updateEndT = min(currentPos + startT, endT);
							if (lastIntv != nullptr) {
								//tie(intvStartT, intvEndT/*, std::ignore*/) = lastIntv->item;
								if (lastIntv->item.second != updateEndT) {
									lastMainLabelT = updateEndT;
									localNoise = 0;
								}
								else {
									lastMainLabelT = lastIntv->item.first - 1;
								}
							}
							else {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
						}
						if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
							int checkPos = (lastMainLabelT - startT) *nEdge + i;
							lazyUpdate(lastMainLabelT, checkPos, i);//update aft
							int labelEnd = lastMainLabelT - startT + max(aft[checkPos], 1) - 1;
							//dynamic
							if (localNoise <= Setting::c) {
								checkPos = labelEnd * nEdge + i;
								lazyUpdate(labelEnd, checkPos, i);//update aft
								if (aft[checkPos] != -MYINFINITE) {
									if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
										newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
										MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
										newEIntR->emplace_back(rd4);
										newValidMidResult.emplace_back(true);
									}
								}
								else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
									newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
									MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
									newEIntR->emplace_back(rd4);
									newValidMidResult.emplace_back(true);
								}
							}
							else {
								if (newPosInEIntR[mainLabelPos] >= 0) {
									newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
									newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
								}
							}
							
							if (lastMainLabelT >= intvE) {
								if (lastMainLabelT >= oriEndTE)
									posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R# for edge i is changed
								else if (posInEIntR[mainLabelPos] != EMaxIntvlChange::CHANGED) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R# for edge i is unchanged 

								if (lastMainLabelT == maxIntv[i].second) {
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
									}
									savePos = lastMainLabelT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
								else if (lastMainLabelT < maxIntv[i].second) {
									int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										else {
											int tempTimestampPos = (intv.second - startT)*nEdge + i;
											if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
												preMaxEMaxIntvlEndT = intv.second;
												preMaxEMaxIntvlStartT = intv.first;
											}
										}
									}

									if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
										if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
										if (EMaxIntvlEndT < lastMainLabelT) {
											maxIntv[i].first = intvB;
											maxIntv[i].second = lastMainLabelT;
										}
										else {
											maxIntv[i].first = EMaxIntvlStartT;
											maxIntv[i].second = EMaxIntvlEndT;
										}
									}

									savePos = lastMainLabelT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
								else {
									int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										else {
											int tempTimestampPos = (intv.second - startT)*nEdge + i;

											if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
												preMaxEMaxIntvlEndT = intv.second;
												preMaxEMaxIntvlStartT = intv.first;
											}
										}
									}
									if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
										preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									}

									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}

									savePos = lastMainLabelT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
							}

							continue;
						}
						else {
							currentPos = max(currentPos + 1, beginPos);
							lastMainLabelPos = lastMainLabelT - startT;
						}

					}

					//Test::comprp2t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
				}

				forbidTimeStartT = -1;
				auto tail = now.tail;
				if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
					forbidTimeStartT = tail->item.first;
					if (tail->pre != nullptr)
						now.deleteNextNode(tail->pre);
					else now.deleteFirstNode();
				}
			}
			else {
				localNoise = 0;
				currentPos = lastMainLabelPos = beginPos;
				forbidTimeStartT = -1;
				noiseNum = 0;
			}
		}

		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;

		lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp){
					lazyUpdate(currentPos - 1, tempTimestampPosForIT - nEdge, i);//update aft
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
				}
			}
			else {
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
				}
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;

		if (lastMainLabelPos >= checkEndPos)
			posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R# for edge i is changed
		else if (posInEIntR[mainLabelPos] != EMaxIntvlChange::CHANGED) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R# for edge i is unchanged 
		if (lastMainLabelPos < endPos) {

			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			
			continue;
		}
		//if ((isEdgeTypeFixed && fixLabel.find(idToLabel[mainLabel[i]]) == fixLabelEnd)) continue;
		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			}

		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}
			//else now.addItemAtLast(make_pair(currentT + 1, maxIntv[i].second));

			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			}
		}
		else {
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
		

			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}

			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			}
		}
	}
}
#pragma endregion

#pragma region FRTMOpt1

void TGraphUDEL::edgeFilterOpt1(int intvB, int intvE,
	vec(int)*& edgeSetsRAdd, vec(int)*& edgeSetsR, int& selectedNum,  int& rightEndpoint, 
	int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;

	int timestampPosForEMaxIntvlEndT;
	for (int i = 0; i < nEdge; i++, labelPosForEdge += numOfLabel) {//O(|E|)
		//auto beginTest = std::chrono::steady_clock::now();
		
		mainLabel = lab[posUsedForEdgeFilter + i];//mainL = L^intvB(e)
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		mainLabelPos = labelPosForEdge + mainLabel;//position of <e,mainL>

		auto& now = vioT[mainLabelPos];//tabuT[e,mainL] before updated
		currentPos = scanT[mainLabelPos];//scanT[e,mainL] before updated
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];//maxIntv[e,mainL] 
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		
		if (currentPos != 0) {

			if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {// L^intvB(e) = L^(intvB-1)(e)
				if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid R sets for e
					continue;
				}
				else {
					auto intvItem = now.first;

					bool shrink = false;
					while (intvItem != nullptr) {//check tabuT[e,mainL]

						if (intvItem->item.first > maxIntv[i].second) break;//no noise
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							savePos = tempEndT - intvE;
							if (savePos >= 0) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}

						int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
						if (lab[tempEWPos] == mainLabel) {
							intvItem->item.second++;
						}
						else {
							int labelSum = intvItem->item.second - intvB + 2;
							if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
								intvItem->item.second++;
							}
						}

						if (intvItem->next != nullptr) {
							if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
								intvItem->item.second = intvItem->next->item.second;
								now.deleteNextNode(intvItem);
							}
							else intvItem = intvItem->next;
						}
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {//maintain maxIntv[e,mainL]
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}
						else intvItem = intvItem->next;
					}

					if (!shrink) {
						savePos = maxIntv[i].second - intvE;
						if (savePos >= 0) {//R set does not change
							for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
								temp++) {
								auto& intv = preEMaxIntvlPtr.q[temp];
								if (intv.second < intvE) {//no overlap
									preEMaxIntvlPtr.swapToTop(temp);
									preEMaxIntvlPtr.pop();
								}
								//else break;
							}
							selectedNum++;
							if (rightEndpoint < savePos) rightEndpoint = savePos;
							edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
						}
					}
				}

				//Test::comprp1n++;
				//Test::comprp1t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
				
				now.removeNodeLessThan(intvE);
				continue;
			}
			else {
				//Test::comprp2n++;
				//now.vioT.removeNodeLessThan(intvE);
				now.removeNodeLessThan(intvE);

				if (maxIntv[i].second == -1) {
					EMaxIntvlEndT = -1;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					EMaxIntvlEndT = maxIntv[i].second;
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
				}
				lastMainLabelT = max(EMaxIntvlEndT, intvB);
				if (currentPos < beginPos) { //does not scan
					//now.localNoise = 0;
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					noiseNum = 0;
				}
				else {
					auto forbidIntv = now.first;
					if (forbidIntv == nullptr) {//does not need update
						if (currentPos != currNTimestamp) {
							int timestampPosForEL = currentPos * nEdge + i;
							edgeType = lab[timestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[timestampPosForEL], 1);
								int tempTPos = timestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = currentPos - eLvalue;
								localNoise = eLvalue;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								localNoise = 0;
								noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
							}
						}
					}
					else {
						//update vioT

						//get noiseNum, localNoise
						int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
						int tempPos = tempT - startT, intvLen = tempT - intvB;
						labelsSum = tempT - intvB;
						int tempTimestampPosForEL = tempPos * nEdge + i;
						edgeType = lab[tempTimestampPosForEL];
						if (edgeType != mainLabel) {
							int eLvalue = max(bef[tempTimestampPosForEL], 1);
							int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
							int tempPosForLabels = tempPos - eLvalue;
							localNoise = eLvalue - 1;
							while (lab[tempTPos] != mainLabel) {
								eLvalue = max(bef[tempTPos], 1);
								tempPosForLabels -= eLvalue;
								localNoise += eLvalue;
								tempTPos -= eLvalue * nEdge;
							}
							noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
						}
						else {
							noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
							localNoise = 0;
						}

						auto tempIntv = forbidIntv;
						while (forbidIntv != nullptr) {

							forbidTimeStartT = -1;
							int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
							while (tempT <= stopT) {

								intvLen = nextLabT - tempT;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType == mainLabel) {
									localNoise = 0;
									minNoise = Setting::delta * (labelsSum + 1);
									if (LESSEQ(noiseNum, minNoise)) {
										//lastMainLabelPos = nextLabPos - 1;
										//remove = 0;
										if (forbidTimeStartT != -1) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											forbidTimeStartT = -1;
											tempIntv = tempIntv->next;
										}
									}
									else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
										//remove += intvLen;
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
									}
									else {
										//remove = 0;
										//lastMainLabelPos = nextLabPos - 1;
										checkNum = (noiseNum - minNoise) / Setting::delta;
										if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
										else noiseT = (int)checkNum + tempT;
										if (forbidTimeStartT == -1) {
											now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
										}
										else {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
											forbidTimeStartT = -1;
										}
									}
								}
								else {
									if (forbidTimeStartT == -1) {
										forbidTimeStartT = tempT + startT;
									}
									localNoise += intvLen;
									noiseNum += intvLen;
								}
								labelsSum += intvLen;

								tempT = nextLabT;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
									if (tempT <= stopT)
										nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
								}
								else {
									if (tempT <= stopT)
										nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
								}
							}
							if (forbidTimeStartT != -1) {
								//now.vioT.addNodeAft(tempIntv, make_pair(forbidTimeStartT,
								if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
									now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
									tempIntv = tempIntv->next;

									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;

								}
								else {
									tempIntv = forbidIntv = forbidIntv->next;
								}
							}
							else {
								forbidIntv = forbidIntv->pre;
								tempIntv = tempIntv->next;
								if (forbidIntv == nullptr) {
									now.deleteFirstNode();
								}
								else now.deleteNextNode(forbidIntv);
								forbidIntv = tempIntv;
							}

							if (forbidIntv != nullptr) {
								tempT = forbidIntv->item.first;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								labelsSum = tempT - intvB;
								localNoise = 0;
								stopT = forbidIntv->item.second;
							}
						}

						//get the last time where the edge has main label
						auto lastIntv = now.tail;
						int updateEndT = min(currentPos + startT, endT);
						if (lastIntv != nullptr) {
							//tie(intvStartT, intvEndT/*, std::ignore*/) = lastIntv->item;
							if (lastIntv->item.second != updateEndT) {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
							else {
								lastMainLabelT = lastIntv->item.first - 1;
							}
						}
						else {
							lastMainLabelT = updateEndT;
							localNoise = 0;
						}
					}
					if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
						if (lastMainLabelT >= intvE) {
							if (lastMainLabelT == maxIntv[i].second) {
								
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
							}
							else if (lastMainLabelT < maxIntv[i].second) {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}

								if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
									if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									/*if (Test::testingMode == 3) {
										if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
											cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << lastMainLabelT << "," << idToLabel[mainLabel] << "]" << endl;
										}
										else {
											cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << lastMainLabelT << "," << idToLabel[mainLabel] << "]" << endl;
										}
									}*/
									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
								
							}
							else {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;

										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}
								if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
									preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
								}
								/*if (Test::testingMode ==  3 && maxIntv[i].second != -1 && lab[timestampPosForEMaxIntvlEndT] != mainLabel && intvB <= maxIntv[i].second) {
									cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << lastMainLabelT << "," << idToLabel[mainLabel] << "]" << endl;
								}*/
								if (EMaxIntvlEndT < lastMainLabelT) {
									maxIntv[i].first = intvB;
									maxIntv[i].second = lastMainLabelT;
								}
								else {
									maxIntv[i].first = EMaxIntvlStartT;
									maxIntv[i].second = EMaxIntvlEndT;
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
								
							}
						}

						continue;
					}
					else {
						currentPos = max(currentPos + 1, beginPos);
						lastMainLabelPos = lastMainLabelT - startT;
					}
				}

				//Test::comprp2t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
			}

			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
		}
		else {
			localNoise = 0;
			currentPos = lastMainLabelPos = beginPos;
			forbidTimeStartT = -1;
			noiseNum = 0;
		}
		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;


		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp)
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
			}
			else {
				if (currentPos < currNTimestamp)
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;
		
		if (lastMainLabelPos < endPos) {
			continue;
		}

		//if ((isEdgeTypeFixed && fixLabel.find(idToLabel[mainLabel[i]]) == fixLabelEnd)) continue;
		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			/*if (Test::testingMode == 3 && maxIntv[i].second != -1 && lab[timestampPosForEMaxIntvlEndT] != mainLabel && intvB <= maxIntv[i].second) {
				cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
			}*/
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			selectedNum++;
			savePos = lastMainLabelPos - endPos;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				/*if (Test::testingMode == 3) {
					if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
						cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
					}
					else {
						cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
					}
				}*/
				if (maxEMaxIntvlEndT < currentT) {
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
		}
		else {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}
			
			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
		}
	}
}

void TGraphUDEL::edgeFilterOpt1MidR(int intvB, int intvE, vec(int)*& edgeSetsRAdd,/**/
	vec(int)*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;

	int timestampPosForEMaxIntvlEndT;
	for (int i = 0; i < nEdge; i++, labelPosForEdge += numOfLabel) {//O(|E|)
		//auto beginTest = std::chrono::steady_clock::now();

		mainLabel = lab[posUsedForEdgeFilter + i];//mainL = L^intvB(e)
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		mainLabelPos = labelPosForEdge + mainLabel;//position of <e,mainL>

		auto& now = vioT[mainLabelPos];//tabuT[e,mainL] before updated
		currentPos = scanT[mainLabelPos];//scanT[e,mainL] before updated
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];//maxIntv[e,mainL] 
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		if (currentPos != 0) {

			if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {// L^intvB(e) = L^(intvB-1)(e)
				if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid R sets for e
					continue;
				}
				else {
					auto intvItem = now.first;

					bool shrink = false;
					while (intvItem != nullptr) {//check tabuT[e,mainL]

						if (intvItem->item.first > maxIntv[i].second) break;//no noise
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							savePos = tempEndT - intvE;
							if (savePos >= 0) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}

						int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
						if (lab[tempEWPos] == mainLabel) {
							intvItem->item.second++;
						}
						else {
							int labelSum = intvItem->item.second - intvB + 2;
							if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
								intvItem->item.second++;
							}
						}

						if (intvItem->next != nullptr) {
							if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
								intvItem->item.second = intvItem->next->item.second;
								now.deleteNextNode(intvItem);
							}
							else intvItem = intvItem->next;
						}
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {//maintain maxIntv[e,mainL]
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							}
							break;
						}
						else intvItem = intvItem->next;
					}

					if (!shrink) {
						savePos = maxIntv[i].second - intvE;
						if (savePos >= 0) {//R set does not change
							for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
								temp++) {
								auto& intv = preEMaxIntvlPtr.q[temp];
								if (intv.second < intvE) {//no overlap
									preEMaxIntvlPtr.swapToTop(temp);
									preEMaxIntvlPtr.pop();
								}
								//else break;
							}
							selectedNum++;
							if (rightEndpoint < savePos) rightEndpoint = savePos;
							edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
						}
					}
				}

				//Test::comprp1n++;
				//Test::comprp1t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
				now.removeNodeLessThan(intvE);
				continue;
			}
			else {
				//Test::comprp2n++;
				//now.vioT.removeNodeLessThan(intvE);
				now.removeNodeLessThan(intvE);

				if (maxIntv[i].second == -1) {
					EMaxIntvlEndT = -1;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					EMaxIntvlEndT = maxIntv[i].second;
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
				}
				lastMainLabelT = max(EMaxIntvlEndT, intvB);
				if (currentPos < beginPos) { //does not scan
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					noiseNum = 0;
				}
				else {
					auto forbidIntv = now.first;
					if (forbidIntv == nullptr) {//does not need update
						if (currentPos != currNTimestamp) {
							int timestampPosForEL = currentPos * nEdge + i;
							edgeType = lab[timestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[timestampPosForEL], 1);
								int tempTPos = timestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = currentPos - eLvalue;
								localNoise = eLvalue;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								localNoise = 0;
								noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
							}
						}
					}
					else {
						//update vioT

						//get noiseNum, localNoise
						int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
						int tempPos = tempT - startT, intvLen = tempT - intvB;
						labelsSum = tempT - intvB;
						int tempTimestampPosForEL = tempPos * nEdge + i;
						edgeType = lab[tempTimestampPosForEL];
						if (edgeType != mainLabel) {
							int eLvalue = max(bef[tempTimestampPosForEL], 1);
							int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
							int tempPosForLabels = tempPos - eLvalue;
							localNoise = eLvalue - 1;
							while (lab[tempTPos] != mainLabel) {
								eLvalue = max(bef[tempTPos], 1);
								tempPosForLabels -= eLvalue;
								localNoise += eLvalue;
								tempTPos -= eLvalue * nEdge;
							}
							noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
						}
						else {
							noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
							localNoise = 0;
						}

						auto tempIntv = forbidIntv;
						while (forbidIntv != nullptr) {

							forbidTimeStartT = -1;
							lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
							int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
							while (tempT <= stopT) {

								intvLen = nextLabT - tempT;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType == mainLabel) {
									localNoise = 0;
									minNoise = Setting::delta * (labelsSum + 1);
									if (LESSEQ(noiseNum, minNoise)) {
										if (forbidTimeStartT != -1) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											forbidTimeStartT = -1;
											tempIntv = tempIntv->next;
										}
									}
									else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
									}
									else {
										checkNum = (noiseNum - minNoise) / Setting::delta;
										if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
										else noiseT = (int)checkNum + tempT;
										if (forbidTimeStartT == -1) {
											now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
										}
										else {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
											forbidTimeStartT = -1;
										}
									}
								}
								else {
									if (forbidTimeStartT == -1) {
										forbidTimeStartT = tempT + startT;
									}
									localNoise += intvLen;
									noiseNum += intvLen;
								}
								labelsSum += intvLen;

								tempT = nextLabT;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
									if (tempT <= stopT) {
										lazyUpdate(tempT - 1, tempTimestampPosForEL - nEdge, i);//update aft
										nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
									}
								}
								else {
									if (tempT <= stopT)	{
										lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
										nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
									}
								}
							}
							if (forbidTimeStartT != -1) {
								if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
									now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
									tempIntv = tempIntv->next;

									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;

								}
								else {
									tempIntv = forbidIntv = forbidIntv->next;
								}
							}
							else {
								forbidIntv = forbidIntv->pre;
								tempIntv = tempIntv->next;
								if (forbidIntv == nullptr) {
									now.deleteFirstNode();
								}
								else now.deleteNextNode(forbidIntv);
								forbidIntv = tempIntv;
							}

							if (forbidIntv != nullptr) {
								tempT = forbidIntv->item.first;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								labelsSum = tempT - intvB;
								localNoise = 0;
								stopT = forbidIntv->item.second;
							}
						}

						//get the last time where the edge has main label
						auto lastIntv = now.tail;
						int updateEndT = min(currentPos + startT, endT);
						if (lastIntv != nullptr) {
							if (lastIntv->item.second != updateEndT) {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
							else {
								lastMainLabelT = lastIntv->item.first - 1;
							}
						}
						else {
							lastMainLabelT = updateEndT;
							localNoise = 0;
						}
					}
					if (currentPos + 1 == currNTimestamp/* || MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
						int checkPos = (lastMainLabelT - startT)*nEdge + i;
						lazyUpdate(lastMainLabelT, checkPos, i);//update aft
						int labelEnd = lastMainLabelT - startT + max(aft[checkPos], 1) - 1;
						//dynamic
						if (localNoise <= Setting::c) {
							checkPos = labelEnd * nEdge + i;
							lazyUpdate(labelEnd, checkPos, i);//update aft
							if (aft[checkPos] != -MYINFINITE) {
								if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
									newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
									MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
									newEIntR->emplace_back(rd4);
									newValidMidResult.emplace_back(true);
								}
							}
							else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
								newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
								MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
								newEIntR->emplace_back(rd4);
								newValidMidResult.emplace_back(true);
							}
						}
						else {
							if (newPosInEIntR[mainLabelPos] >= 0) {
								newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
								newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
							}
						}

						if (lastMainLabelT >= intvE) {
							if (lastMainLabelT == maxIntv[i].second) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value

							}
							else if (lastMainLabelT < maxIntv[i].second) {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}

								if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
									if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value

							}
							else {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;

										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}
								if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
									preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
								}

								if (EMaxIntvlEndT < lastMainLabelT) {
									maxIntv[i].first = intvB;
									maxIntv[i].second = lastMainLabelT;
								}
								else {
									maxIntv[i].first = EMaxIntvlStartT;
									maxIntv[i].second = EMaxIntvlEndT;
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
							}
						}

						continue;
					}
					else {
						currentPos = max(currentPos + 1, beginPos);
						lastMainLabelPos = lastMainLabelT - startT;
					}
				}


				//Test::comprp2t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
			}

			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
		}
		else {
			localNoise = 0;
			currentPos = lastMainLabelPos = beginPos;
			forbidTimeStartT = -1;
			noiseNum = 0;
		}
		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;


		lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos - 1, tempTimestampPosForIT - nEdge, i);//update aft
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
				}
			}
			else {
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
				}
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;

		if (lastMainLabelPos < endPos) {
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			continue;
		}

		//if ((isEdgeTypeFixed && fixLabel.find(idToLabel[mainLabel[i]]) == fixLabelEnd)) continue;
		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			selectedNum++;
			savePos = lastMainLabelPos - endPos;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
		}
		else {
			//dynamic
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}

			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}
			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
		}
	}
}

void TGraphUDEL::edgeFilterOpt1MidRForDYN(int intvB, int intvE, int oriEndTE,
	vec(int)*& edgeSetsRAdd,
	vec(int)*& edgeSetsR, int& selectedNum,
	int& rightEndpoint,
	int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	//int eLen, eLvalue;
	//unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int beginPos = intvB - startT;
	int endPos = intvE - startT, checkEndPos = oriEndTE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;

	//bool addItem;
	//int oriEndT, oriStartT;
	int timestampPosForEMaxIntvlEndT;
	auto fromMidREnd = edgesInEIntR->end(), fromMidRIter = edgesInEIntR->begin();

	for (; fromMidRIter != fromMidREnd; ++fromMidRIter) {

		int i = *fromMidRIter;
		labelPosForEdge = i * numOfLabel;

		mainLabel = lab[posUsedForEdgeFilter + i];//mainL = L^intvB(e)
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		mainLabelPos = labelPosForEdge + mainLabel;//position of <e,mainL>
		auto& now = vioT[mainLabelPos];//tabuT[e,mainL] before updated
		currentPos = scanT[mainLabelPos];//scanT[e,mainL] before updated
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];//maxIntv[e,mainL] 
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;

		if (posInEIntR[mainLabelPos] == EMaxIntvlChange::INTVINIT) {//initial 
			continue;
		}
		else if (posInEIntR[mainLabelPos] >= 0) {//continue to scan
			if (currentPos == 0) continue; //initial

			//recover localNoise , noiseNum and lastMainLabelPos
			int timestampPosForEL = currentPos * nEdge + i;
			edgeType = lab[timestampPosForEL];
			if (edgeType != mainLabel) {
				int eLvalue = max(bef[timestampPosForEL], 1);
				int tempTPos = timestampPosForEL - eLvalue * nEdge;
				int tempPosForLabels = currentPos - eLvalue;
				localNoise = eLvalue;
				while (lab[tempTPos] != mainLabel) {
					eLvalue = max(bef[tempTPos], 1);
					tempPosForLabels -= eLvalue;
					localNoise += eLvalue;
					tempTPos -= eLvalue * nEdge;
				}
				noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
			}
			else {
				localNoise = 0;
				noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
			}

			if (maxIntv[i].second == -1) {
				lastMainLabelPos = beginPos;
			}
			else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
				lastMainLabelPos = max(maxIntv[i].second, beginPos);
			}
			else {
				EMaxIntvlEndT = -1;
				int last = preEMaxIntvlPtr.rear;
				if (last != preEMaxIntvlPtr.front) {
					int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
					for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
						EMaxIntvlEndT = preEMaxIntvlPtr.q[temp].second;
						if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
					}
					if (temp == preEMaxIntvlPtr.front) {
						EMaxIntvlEndT = preEMaxIntvlPtr.q[temp].second;
						if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
					}
				}
				lastMainLabelPos = max(EMaxIntvlEndT, beginPos);
			}

			//recover vioT
			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				lastMainLabelPos = forbidTimeStartT - 1;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
			currentPos = max(currentPos + 1, beginPos);
		}
		else {
			if (currentPos != 0) {

				if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {// L^intvB(e) = L^(intvB-1)(e)
					if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid R sets for e
						continue;
					}
					else {
						auto intvItem = now.first;

						bool shrink = false;
						while (intvItem != nullptr) {//check tabuT[e,mainL]

							if (intvItem->item.first > maxIntv[i].second) break;//no noise
							else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
								shrink = true;
								int tempEndT = intvItem->item.first - 1;
								savePos = tempEndT - intvE;
								if (savePos >= 0) {
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
										temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										//else break;
									}
									savePos = tempEndT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
								break;
							}

							int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
							if (lab[tempEWPos] == mainLabel) {
								intvItem->item.second++;
							}
							else {
								int labelSum = intvItem->item.second - intvB + 2;
								if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
									intvItem->item.second++;
								}
							}

							if (intvItem->next != nullptr) {
								if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
									intvItem->item.second = intvItem->next->item.second;
									now.deleteNextNode(intvItem);
								}
								else intvItem = intvItem->next;
							}
							else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
								shrink = true;
								int tempEndT = intvItem->item.first - 1;
								if (tempEndT >= intvE) {
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
										temp++) {//maintain maxIntv[e,mainL]
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
									}
									savePos = tempEndT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
									}
								}
								break;
							}
							else intvItem = intvItem->next;
						}

						if (!shrink) {
							if (maxIntv[i].second >= intvE) {//R set does not change
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < intvE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								savePos = maxIntv[i].second - oriEndTE;
								if (savePos >= 0) {
									selectedNum++;
									if (rightEndpoint < savePos) rightEndpoint = savePos;
									edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								}
							}
						}
					}

					now.removeNodeLessThan(intvE);
					continue;
				}
				else {
					//now.vioT.removeNodeLessThan(intvE);
					now.removeNodeLessThan(intvE);

					if (maxIntv[i].second == -1) {
						EMaxIntvlEndT = -1;
					}
					else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
						EMaxIntvlEndT = maxIntv[i].second;
					}
					else {
						EMaxIntvlEndT = -1;
						int last = preEMaxIntvlPtr.rear;
						if (last != preEMaxIntvlPtr.front) {
							int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
							for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
								tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
								if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
							}
							if (temp == preEMaxIntvlPtr.front) {
								tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
								if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
							}
						}
					}
					lastMainLabelT = max(EMaxIntvlEndT, intvB);
					if (currentPos < beginPos) { //does not scan
						localNoise = 0;
						currentPos = lastMainLabelPos = beginPos;
						noiseNum = 0;
					}
					else {
						auto forbidIntv = now.first;
						if (forbidIntv == nullptr) {//does not need update
							if (currentPos != currNTimestamp) {
								int timestampPosForEL = currentPos * nEdge + i;
								edgeType = lab[timestampPosForEL];
								if (edgeType != mainLabel) {
									int eLvalue = max(bef[timestampPosForEL], 1);
									int tempTPos = timestampPosForEL - eLvalue * nEdge;
									int tempPosForLabels = currentPos - eLvalue;
									localNoise = eLvalue;
									while (lab[tempTPos] != mainLabel) {
										eLvalue = max(bef[tempTPos], 1);
										tempPosForLabels -= eLvalue;
										localNoise += eLvalue;
										tempTPos -= eLvalue * nEdge;
									}
									noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
								}
								else {
									localNoise = 0;
									noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
								}
							}
						}
						else {
							//update vioT

							//get noiseNum, localNoise
							int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
							int tempPos = tempT - startT, intvLen = tempT - intvB;
							labelsSum = tempT - intvB;
							int tempTimestampPosForEL = tempPos * nEdge + i;
							edgeType = lab[tempTimestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[tempTimestampPosForEL], 1);
								int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = tempPos - eLvalue;
								localNoise = eLvalue - 1;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
								localNoise = 0;
							}

							auto tempIntv = forbidIntv;
							while (forbidIntv != nullptr) {

								forbidTimeStartT = -1;
								lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
								int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
								while (tempT <= stopT) {

									intvLen = nextLabT - tempT;
									edgeType = lab[tempTimestampPosForEL];
									if (edgeType == mainLabel) {
										localNoise = 0;
										minNoise = Setting::delta * (labelsSum + 1);
										if (LESSEQ(noiseNum, minNoise)) {
											if (forbidTimeStartT != -1) {
												now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
												forbidTimeStartT = -1;
												tempIntv = tempIntv->next;
											}
										}
										else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
											if (forbidTimeStartT == -1) {
												forbidTimeStartT = tempT + startT;
											}
										}
										else {
											//remove = 0;
											//lastMainLabelPos = nextLabPos - 1;
											checkNum = (noiseNum - minNoise) / Setting::delta;
											if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
											else noiseT = (int)checkNum + tempT;
											if (forbidTimeStartT == -1) {
												now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
												tempIntv = tempIntv->next;
											}
											else {
												now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
												tempIntv = tempIntv->next;
												forbidTimeStartT = -1;
											}
										}
									}
									else {
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
										localNoise += intvLen;
										noiseNum += intvLen;
									}
									labelsSum += intvLen;

									tempT = nextLabT;
									tempTimestampPosForEL = (tempT - startT) * nEdge + i;
									if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
										if (tempT <= stopT) {
											lazyUpdate(tempT - 1, tempTimestampPosForEL - nEdge, i);//update aft
											nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
										}
									}
									else {
										if (tempT <= stopT) {
											lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
											nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
										}
									}
								}
								if (forbidTimeStartT != -1) {
									if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
										now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
										tempIntv = tempIntv->next;

										forbidIntv = forbidIntv->pre;
										tempIntv = tempIntv->next;
										if (forbidIntv == nullptr) {
											now.deleteFirstNode();
										}
										else now.deleteNextNode(forbidIntv);
										forbidIntv = tempIntv;

									}
									else {
										tempIntv = forbidIntv = forbidIntv->next;
									}
								}
								else {
									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;
								}

								if (forbidIntv != nullptr) {
									tempT = forbidIntv->item.first;
									tempTimestampPosForEL = (tempT - startT) * nEdge + i;
									labelsSum = tempT - intvB;
									localNoise = 0;
									stopT = forbidIntv->item.second;
								}
							}

							//get the last time where the edge has main label
							auto lastIntv = now.tail;
							int updateEndT = min(currentPos + startT, endT);
							if (lastIntv != nullptr) {
								if (lastIntv->item.second != updateEndT) {
									lastMainLabelT = updateEndT;
									localNoise = 0;
								}
								else {
									lastMainLabelT = lastIntv->item.first - 1;
								}
							}
							else {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
						}
						if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
							int checkPos = (lastMainLabelT - startT)*nEdge + i;
							lazyUpdate(lastMainLabelT, checkPos, i);//update aft
							int labelEnd = lastMainLabelT - startT + max(aft[checkPos], 1) - 1;
							//dynamic
							if (localNoise <= Setting::c) {
								checkPos = labelEnd * nEdge + i;
								lazyUpdate(labelEnd, checkPos, i);//update aft
								if (aft[checkPos] != -MYINFINITE) {
									if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
										newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
										MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
										newEIntR->emplace_back(rd4);
										newValidMidResult.emplace_back(true);
									}
								}
								else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
									newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
									MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
									newEIntR->emplace_back(rd4);
									newValidMidResult.emplace_back(true);
								}
							}
							else {
								if (newPosInEIntR[mainLabelPos] >= 0) {
									newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
									newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
								}
							}

							if (lastMainLabelT >= intvE) {
								if (lastMainLabelT >= oriEndTE)
									posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R# for edge i is changed
								else if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R# for edge i is unchanged 

								if (lastMainLabelT == maxIntv[i].second) {
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
									}
									savePos = lastMainLabelT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
										edgeSetsRAdd[savePos].emplace_back(i);// add current row value
									}
								}
								else if (lastMainLabelT < maxIntv[i].second) {
									int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										else {
											int tempTimestampPos = (intv.second - startT)*nEdge + i;
											if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
												preMaxEMaxIntvlEndT = intv.second;
												preMaxEMaxIntvlStartT = intv.first;
											}
										}
									}

									if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
										if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
										if (EMaxIntvlEndT < lastMainLabelT) {
											maxIntv[i].first = intvB;
											maxIntv[i].second = lastMainLabelT;
										}
										else {
											maxIntv[i].first = EMaxIntvlStartT;
											maxIntv[i].second = EMaxIntvlEndT;
										}
									}

									savePos = lastMainLabelT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
										edgeSetsRAdd[savePos].emplace_back(i);// add current row value
									}
								}
								else {
									int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < intvE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										else {
											int tempTimestampPos = (intv.second - startT)*nEdge + i;

											if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
												preMaxEMaxIntvlEndT = intv.second;
												preMaxEMaxIntvlStartT = intv.first;
											}
										}
									}
									if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
										preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									}

									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}

									savePos = lastMainLabelT - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
										edgeSetsRAdd[savePos].emplace_back(i);// add current row value
									}
								}
							}

							continue;
						}
						else {
							currentPos = max(currentPos + 1, beginPos);
							lastMainLabelPos = lastMainLabelT - startT;
						}
					}

				}

				forbidTimeStartT = -1;
				auto tail = now.tail;
				if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
					forbidTimeStartT = tail->item.first;
					if (tail->pre != nullptr)
						now.deleteNextNode(tail->pre);
					else now.deleteFirstNode();
				}
			}
			else {
				localNoise = 0;
				currentPos = lastMainLabelPos = beginPos;
				forbidTimeStartT = -1;
				noiseNum = 0;
			}
		}
		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;


		lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos - 1, tempTimestampPosForIT - nEdge, i);//update aft
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
				}
			}
			else {
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
				}
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;

		if (lastMainLabelPos >= checkEndPos)
			posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R# for edge i is changed
		else if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R# for edge i is unchanged 
		if (lastMainLabelPos < endPos) {
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			continue;
		}

		//if ((isEdgeTypeFixed && fixLabel.find(idToLabel[mainLabel[i]]) == fixLabelEnd)) continue;
		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
				edgeSetsRAdd[savePos].emplace_back(i);// add current row value
			}
		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= intvE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}
			//else now.addItemAtLast(make_pair(currentT + 1, maxIntv[i].second));

			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
				edgeSetsRAdd[savePos].emplace_back(i);// add current row value
			}
		}
		else {
			//dynamic
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}

			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < intvE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= intvE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}
			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
				edgeSetsRAdd[savePos].emplace_back(i);// add current row value
			}
		}

	}
}
#pragma endregion


#pragma region FRTMPLUS
void TGraphUDEL::edgeFilterPlus(int intvB, int intvE, int filterE, int limited, bool*&hasE, int*& newE, int& newENum, vec(int)*& edgeSetsRAdd,
	vec(int)*& edgeSetsR, int& selectedNum, vec(int)*& edgeSetsRShortIntv, int& selectedNumShortIntv, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT;
	int timePos = filterE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen, intvLenNoNoise = filterE - intvB + 1;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;
	bool saveNoNoise;

	int timestampPosForEMaxIntvlEndT;
	for (int i = 0; i < nEdge; i++, labelPosForEdge += numOfLabel) {//O(|E|)
		
		mainLabel = lab[posUsedForEdgeFilter + i];//mainL = L^intvB(e)
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;

		saveNoNoise = false;
		edgeType = lab[TGraph::posUsedForEdgeFilterShortIntv + i];//label at intvE
		if (max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) >= intvLenNoNoise) {
			int intvStart = filterE - max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) + 1;
			int j = timePos + 1;
			if (maxIntvShortIntv[i].second >= timePos && maxIntvShortIntv[i].first <= beginPos) {
				int check = min(maxIntvShortIntv[i].second - timePos, limited);
				selectedNumShortIntv++;
				edgeSetsRShortIntv[check].emplace_back(i);
				saveNoNoise = true;
			}
			else if (maxIntvShortIntv[i].second >= intvB - startT && maxIntvShortIntv[i].first <= timePos) {
			}
			else {
				maxIntvShortIntv[i].first = intvStart;
				selectedNumShortIntv++;

				int check = max(aft[TGraph::posUsedForEdgeFilterShortIntv + i], 1) - 1;
				maxIntvShortIntv[i].second = timePos + check;
				check = min(check, limited);
				edgeSetsRShortIntv[check].emplace_back(i);
				saveNoNoise = true;
			}
		}
		
		mainLabelPos = labelPosForEdge + mainLabel;//position of <e,mainL>

		auto& now = vioT[mainLabelPos];//tabuT[e,mainL] before updated
		currentPos = scanT[mainLabelPos];//scanT[e,mainL] before updated
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];//maxIntv[e,mainL] 
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;

		if (currentPos != 0) {
			if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {// L^intvB(e) = L^(intvB-1)(e)
				if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid R sets for e
					if (saveNoNoise) {
						newE[newENum++] = i;
					}
					continue;
				}
				else {
					auto intvItem = now.first;

					bool shrink = false;
					while (intvItem != nullptr) {//check tabuT[e,mainL]

						if (intvItem->item.first > maxIntv[i].second) break;//no noise
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							savePos = tempEndT - intvE;
							if (savePos >= 0) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									//else break;
								}
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i);
								hasE[i] = true;
							}
							break;
						}

						int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
						if (lab[tempEWPos] == mainLabel) {
							intvItem->item.second++;
						}
						else {
							int labelSum = intvItem->item.second - intvB + 2;
							if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
								intvItem->item.second++;
							}
						}

						if (intvItem->next != nullptr) {
							if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
								intvItem->item.second = intvItem->next->item.second;
								now.deleteNextNode(intvItem);
							}
							else intvItem = intvItem->next;
						}
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {//maintain maxIntv[e,mainL]
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i);
								hasE[i] = true;
							}
							break;
						}
						else intvItem = intvItem->next;
					}

					if (!shrink) {
						savePos = maxIntv[i].second - intvE;
						if (savePos >= 0) {//R set does not change
							for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
								temp++) {
								auto& intv = preEMaxIntvlPtr.q[temp];
								if (intv.second < filterE) {//no overlap
									preEMaxIntvlPtr.swapToTop(temp);
									preEMaxIntvlPtr.pop();
								}
								//else break;
							}
							selectedNum++;
							if (rightEndpoint < savePos) rightEndpoint = savePos;
							edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
							hasE[i] = true;
						}
					}
				}
				now.removeNodeLessThan(intvE);

				if (!hasE[i] && saveNoNoise) {
					newE[newENum++] = i;
				}
				continue;
			}
			else {
				now.removeNodeLessThan(intvE);

				if (maxIntv[i].second == -1) {
					EMaxIntvlEndT = -1;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					EMaxIntvlEndT = maxIntv[i].second;
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
				}
				lastMainLabelT = max(EMaxIntvlEndT, intvB);
				if (currentPos < beginPos) { //does not scan
					//now.localNoise = 0;
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					noiseNum = 0;
				}
				else {
					auto forbidIntv = now.first;
					if (forbidIntv == nullptr) {//does not need update
						if (currentPos != currNTimestamp) {
							int timestampPosForEL = currentPos * nEdge + i;
							edgeType = lab[timestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[timestampPosForEL], 1);
								int tempTPos = timestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = currentPos - eLvalue;
								localNoise = eLvalue;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								localNoise = 0;
								noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
							}
						}
					}
					else {
						//update vioT

						//get noiseNum, localNoise
						int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
						int tempPos = tempT - startT, intvLen = tempT - intvB;
						labelsSum = tempT - intvB;
						int tempTimestampPosForEL = tempPos * nEdge + i;
						edgeType = lab[tempTimestampPosForEL];
						if (edgeType != mainLabel) {
							int eLvalue = max(bef[tempTimestampPosForEL], 1);
							int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
							int tempPosForLabels = tempPos - eLvalue;
							localNoise = eLvalue - 1;
							while (lab[tempTPos] != mainLabel) {
								eLvalue = max(bef[tempTPos], 1);
								tempPosForLabels -= eLvalue;
								localNoise += eLvalue;
								tempTPos -= eLvalue * nEdge;
							}
							noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
						}
						else {
							noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
							localNoise = 0;
						}

						auto tempIntv = forbidIntv;
						while (forbidIntv != nullptr) {

							forbidTimeStartT = -1;
							int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
							while (tempT <= stopT) {

								intvLen = nextLabT - tempT;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType == mainLabel) {
									localNoise = 0;
									minNoise = Setting::delta * (labelsSum + 1);
									if (LESSEQ(noiseNum, minNoise)) {
										if (forbidTimeStartT != -1) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											forbidTimeStartT = -1;
											tempIntv = tempIntv->next;
										}
									}
									else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
									}
									else {
										//remove = 0;
										checkNum = (noiseNum - minNoise) / Setting::delta;
										if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
										else noiseT = (int)checkNum + tempT;
										if (forbidTimeStartT == -1) {
											now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
										}
										else {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
											forbidTimeStartT = -1;
										}
									}
								}
								else {
									if (forbidTimeStartT == -1) {
										forbidTimeStartT = tempT + startT;
									}
									localNoise += intvLen;
									noiseNum += intvLen;
								}
								labelsSum += intvLen;

								tempT = nextLabT;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
									if (tempT <= stopT)
										nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
								}
								else {
									if (tempT <= stopT)
										nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
								}
							}
							if (forbidTimeStartT != -1) {
								if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
									now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
									tempIntv = tempIntv->next;

									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;
								}
								else {
									tempIntv = forbidIntv = forbidIntv->next;
								}
							}
							else {
								forbidIntv = forbidIntv->pre;
								tempIntv = tempIntv->next;
								if (forbidIntv == nullptr) {
									now.deleteFirstNode();
								}
								else now.deleteNextNode(forbidIntv);
								forbidIntv = tempIntv;
							}

							if (forbidIntv != nullptr) {
								tempT = forbidIntv->item.first;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								labelsSum = tempT - intvB;
								localNoise = 0;
								stopT = forbidIntv->item.second;
							}
						}

						//get the last time where the edge has main label
						auto lastIntv = now.tail;
						int updateEndT = min(currentPos + startT, endT);
						if (lastIntv != nullptr) {
							if (lastIntv->item.second != updateEndT) {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
							else {
								lastMainLabelT = lastIntv->item.first - 1;
							}
						}
						else {
							lastMainLabelT = updateEndT;
							localNoise = 0;
						}
					}
					if (currentPos + 1 == currNTimestamp/* || MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
						if (lastMainLabelT >= intvE) {
							if (lastMainLabelT == maxIntv[i].second) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);
								hasE[i] = true;
							}
							else if (lastMainLabelT < maxIntv[i].second) {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}

								if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
									if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= filterE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);
								hasE[i] = true;
							}
							else {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;

										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}
								if (maxIntv[i].second >= filterE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
									preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
								}

								if (EMaxIntvlEndT < lastMainLabelT) {
									maxIntv[i].first = intvB;
									maxIntv[i].second = lastMainLabelT;
								}
								else {
									maxIntv[i].first = EMaxIntvlStartT;
									maxIntv[i].second = EMaxIntvlEndT;
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);
								hasE[i] = true;
							}
						}

						if (!hasE[i] && saveNoNoise) {
							newE[newENum++] = i;
						}
						continue;
					}
					else {
						currentPos = max(currentPos + 1, beginPos);
						lastMainLabelPos = lastMainLabelT - startT;
					}
				}
			}

			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
		}
		else {
			localNoise = 0;
			currentPos = lastMainLabelPos = beginPos;
			forbidTimeStartT = -1;
			noiseNum = 0;
		}
		labelsSum = currentPos - beginPos;

        auto a = std::chrono::steady_clock::now(); 
		int tempTimestampPosForIT = currentPos * nEdge + i;

		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(intvE, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(intvE, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(intvE, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp)
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
			}
			else {
				if (currentPos < currNTimestamp)
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;

		//if (intvB == 0) {
		//	long long b = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count();
		//	//cout << b << endl;
		//	Test::counter += b;
		//}
		if (lastMainLabelPos < endPos) {
			if (saveNoNoise) {
				newE[newENum++] = i;
			}
			continue;
		}

		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);
			hasE[i] = true;
		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= filterE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << idToLabel[lab[timestampPosForEMaxIntvlEndT]] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}
			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);
			hasE[i] = true;
		}
		else {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= filterE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);
			hasE[i] = true;
		}

		if (!hasE[i] && saveNoNoise) {
			newE[newENum++] = i;
		}
	}
}

void TGraphUDEL::edgeFilterPlusMidR(int intvB, int intvE, int filterE, int limited, bool*&/*iSet&*/hasE, int*& newE, int& newENum, vec(int)*& edgeSetsRAdd,
	vec(int)*& edgeSetsR, int& selectedNum, vec(int)*& edgeSetsRShortIntv, int& selectedNumShortIntv, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT;
	int timePos = filterE - startT;
	int currentPos, mainLabelPos, currentT;
	int labelsSum, intvLen, intvLenNoNoise = filterE - intvB + 1;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge = 0;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;
	bool saveNoNoise;

	int timestampPosForEMaxIntvlEndT;
	for (int i = 0; i < nEdge; i++, labelPosForEdge += numOfLabel) {//O(|E|)
		mainLabel = lab[posUsedForEdgeFilter + i];//mainL = L^intvB(e)
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;
		
		saveNoNoise = false;
		edgeType = lab[TGraph::posUsedForEdgeFilterShortIntv + i];//label at intvE
		if (max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) >= intvLenNoNoise) {
			int intvStart = filterE - max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) + 1;
			int j = timePos + 1;
			if (maxIntvShortIntv[i].second >= timePos && maxIntvShortIntv[i].first <= beginPos) {
				int check = min(maxIntvShortIntv[i].second - timePos, limited);
				selectedNumShortIntv++;
				edgeSetsRShortIntv[check].emplace_back(i);
				saveNoNoise = true;
			}
			else if (maxIntvShortIntv[i].second >= intvB - startT && maxIntvShortIntv[i].first <= timePos) {
			}
			else {
				maxIntvShortIntv[i].first = intvStart;
				selectedNumShortIntv++;
				
				lazyUpdate(filterE, TGraph::posUsedForEdgeFilterShortIntv + i, i);//update aft
				int check = max(aft[TGraph::posUsedForEdgeFilterShortIntv + i], 1) - 1;
				maxIntvShortIntv[i].second = timePos + check;
				check = min(check, limited);
				edgeSetsRShortIntv[check].emplace_back(i);
				saveNoNoise = true;
			}
		}

		mainLabelPos = labelPosForEdge + mainLabel;//position of <e,mainL>
		
		auto& now = vioT[mainLabelPos];//tabuT[e,mainL] before updated
		currentPos = scanT[mainLabelPos];//scanT[e,mainL] before updated
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];//maxIntv[e,mainL] 
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		
		if (currentPos != 0) {

			if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {// L^intvB(e) = L^(intvB-1)(e)
				if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid R sets for e
					if (saveNoNoise) {
						newE[newENum++] = i;
					}
					continue;
				}
				else {
					auto intvItem = now.first;

					bool shrink = false;
					while (intvItem != nullptr) {//check tabuT[e,mainL]

						if (intvItem->item.first > maxIntv[i].second) break;//no noise
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							savePos = tempEndT - intvE;
							if (savePos >= 0) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i);
								hasE[i] = true;
							}
							break;
						}

						int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
						if (lab[tempEWPos] == mainLabel) {
							intvItem->item.second++;
						}
						else {
							int labelSum = intvItem->item.second - intvB + 2;
							if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
								intvItem->item.second++;
							}
						}

						if (intvItem->next != nullptr) {
							if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
								intvItem->item.second = intvItem->next->item.second;
								now.deleteNextNode(intvItem);
							}
							else intvItem = intvItem->next;
						}
						else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
							shrink = true;
							int tempEndT = intvItem->item.first - 1;
							if (tempEndT >= intvE) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
									temp++) {//maintain maxIntv[e,mainL]
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								savePos = tempEndT - intvE;
								selectedNum++;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i);
								hasE[i] = true;
							}
							break;
						}
						else intvItem = intvItem->next;
					}

					if (!shrink) {
						savePos = maxIntv[i].second - intvE;
						if (savePos >= 0) {//R set does not change
							for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
								temp++) {
								auto& intv = preEMaxIntvlPtr.q[temp];
								if (intv.second < filterE) {//no overlap
									preEMaxIntvlPtr.swapToTop(temp);
									preEMaxIntvlPtr.pop();
								}
							}
							selectedNum++;
							if (rightEndpoint < savePos) rightEndpoint = savePos;
							edgeSetsR[savePos].emplace_back(i);
							hasE[i] = true;
						}
					}
				}
				if (!hasE[i] && saveNoNoise) {
					newE[newENum++] = i;
				}
				now.removeNodeLessThan(intvE);
				continue;
			}
			else {
				//Test::comprp2n++;
				now.removeNodeLessThan(intvE);

				if (maxIntv[i].second == -1) {
					EMaxIntvlEndT = -1;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					EMaxIntvlEndT = maxIntv[i].second;
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
				}
				lastMainLabelT = max(EMaxIntvlEndT, intvB);
				if (currentPos < beginPos) { //does not scan
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					noiseNum = 0;
				}
				else {
					auto forbidIntv = now.first;
					if (forbidIntv == nullptr) {//does not need update
						if (currentPos != currNTimestamp) {
							int timestampPosForEL = currentPos * nEdge + i;
							edgeType = lab[timestampPosForEL];
							if (edgeType != mainLabel) {
								int eLvalue = max(bef[timestampPosForEL], 1);
								int tempTPos = timestampPosForEL - eLvalue * nEdge;
								int tempPosForLabels = currentPos - eLvalue;
								localNoise = eLvalue;
								while (lab[tempTPos] != mainLabel) {
									eLvalue = max(bef[tempTPos], 1);
									tempPosForLabels -= eLvalue;
									localNoise += eLvalue;
									tempTPos -= eLvalue * nEdge;
								}
								noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
							}
							else {
								localNoise = 0;
								noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
							}
						}
					}
					else {
						//update vioT

						//get noiseNum, localNoise
						int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
						int tempPos = tempT - startT, intvLen = tempT - intvB;
						labelsSum = tempT - intvB;
						int tempTimestampPosForEL = tempPos * nEdge + i;
						edgeType = lab[tempTimestampPosForEL];
						if (edgeType != mainLabel) {
							int eLvalue = max(bef[tempTimestampPosForEL], 1);
							int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
							int tempPosForLabels = tempPos - eLvalue;
							localNoise = eLvalue - 1;
							while (lab[tempTPos] != mainLabel) {
								eLvalue = max(bef[tempTPos], 1);
								tempPosForLabels -= eLvalue;
								localNoise += eLvalue;
								tempTPos -= eLvalue * nEdge;
							}
							noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
						}
						else {
							noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
							localNoise = 0;
						}

						auto tempIntv = forbidIntv;
						while (forbidIntv != nullptr) {

							forbidTimeStartT = -1;
							lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
							int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
							while (tempT <= stopT) {

								intvLen = nextLabT - tempT;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType == mainLabel) {
									localNoise = 0;
									minNoise = Setting::delta * (labelsSum + 1);
									if (LESSEQ(noiseNum, minNoise)) {
										if (forbidTimeStartT != -1) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											forbidTimeStartT = -1;
											tempIntv = tempIntv->next;
										}
									}
									else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
										if (forbidTimeStartT == -1) {
											forbidTimeStartT = tempT + startT;
										}
									}
									else {
										checkNum = (noiseNum - minNoise) / Setting::delta;
										if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
										else noiseT = (int)checkNum + tempT;
										if (forbidTimeStartT == -1) {
											now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
										}
										else {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
											tempIntv = tempIntv->next;
											forbidTimeStartT = -1;
										}
									}
								}
								else {
									if (forbidTimeStartT == -1) {
										forbidTimeStartT = tempT + startT;
									}
									localNoise += intvLen;
									noiseNum += intvLen;
								}
								labelsSum += intvLen;

								tempT = nextLabT;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
									if (tempT <= stopT) {
										lazyUpdate(tempT - 1, tempTimestampPosForEL - nEdge, i);//update aft
										nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
									}
								}
								else {
									if (tempT <= stopT) {
										lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
										nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
									}
								}
							}
							if (forbidTimeStartT != -1) {
								if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
									now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
									tempIntv = tempIntv->next;

									forbidIntv = forbidIntv->pre;
									tempIntv = tempIntv->next;
									if (forbidIntv == nullptr) {
										now.deleteFirstNode();
									}
									else now.deleteNextNode(forbidIntv);
									forbidIntv = tempIntv;

								}
								else {
									tempIntv = forbidIntv = forbidIntv->next;
								}
							}
							else {
								forbidIntv = forbidIntv->pre;
								tempIntv = tempIntv->next;
								if (forbidIntv == nullptr) {
									now.deleteFirstNode();
								}
								else now.deleteNextNode(forbidIntv);
								forbidIntv = tempIntv;
							}

							if (forbidIntv != nullptr) {
								tempT = forbidIntv->item.first;
								tempTimestampPosForEL = (tempT - startT) * nEdge + i;
								labelsSum = tempT - intvB;
								localNoise = 0;
								stopT = forbidIntv->item.second;
							}
						}

						//get the last time where the edge has main label
						auto lastIntv = now.tail;
						int updateEndT = min(currentPos + startT, endT);
						if (lastIntv != nullptr) {
							if (lastIntv->item.second != updateEndT) {
								lastMainLabelT = updateEndT;
								localNoise = 0;
							}
							else {
								lastMainLabelT = lastIntv->item.first - 1;
							}
						}
						else {
							lastMainLabelT = updateEndT;
							localNoise = 0;
						}
					}
					if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
						int checkPos = (lastMainLabelT - startT)*nEdge + i;
						lazyUpdate(lastMainLabelT, checkPos, i);//update aft
						int labelEnd = lastMainLabelT - startT + max(aft[checkPos], 1) - 1;
						//dynamic
						if (localNoise <= Setting::c) {
							checkPos = labelEnd * nEdge + i;
							lazyUpdate(labelEnd, checkPos, i);//update aft
							if (aft[checkPos] != -MYINFINITE) {
								if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
									newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
									MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
									newEIntR->emplace_back(rd4);
									newValidMidResult.emplace_back(true);
								}
							}
							else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
								newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
								MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
								newEIntR->emplace_back(rd4);
								newValidMidResult.emplace_back(true);
							}
						}
						else {
							if (newPosInEIntR[mainLabelPos] >= 0) {
								newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
								newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
							}
						}

						if (lastMainLabelT >= intvE) {
							if (lastMainLabelT == maxIntv[i].second) {
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
								}
								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
								hasE[i] = true;

							}
							else if (lastMainLabelT < maxIntv[i].second) {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;
										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}

								if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
									if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= filterE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
									if (EMaxIntvlEndT < lastMainLabelT) {
										maxIntv[i].first = intvB;
										maxIntv[i].second = lastMainLabelT;
									}
									else {
										maxIntv[i].first = EMaxIntvlStartT;
										maxIntv[i].second = EMaxIntvlEndT;
									}
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
								hasE[i] = true;
							}
							else {
								int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
								for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
									auto& intv = preEMaxIntvlPtr.q[temp];
									if (intv.second < filterE) {//no overlap
										preEMaxIntvlPtr.swapToTop(temp);
										preEMaxIntvlPtr.pop();
									}
									else {
										int tempTimestampPos = (intv.second - startT)*nEdge + i;

										if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
											preMaxEMaxIntvlEndT = intv.second;
											preMaxEMaxIntvlStartT = intv.first;
										}
									}
								}
								if (maxIntv[i].second >= filterE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
									preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
								}

								if (EMaxIntvlEndT < lastMainLabelT) {
									maxIntv[i].first = intvB;
									maxIntv[i].second = lastMainLabelT;
								}
								else {
									maxIntv[i].first = EMaxIntvlStartT;
									maxIntv[i].second = EMaxIntvlEndT;
								}

								selectedNum++;
								savePos = lastMainLabelT - intvE;
								if (rightEndpoint < savePos) rightEndpoint = savePos;
								edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
								edgeSetsRAdd[savePos].emplace_back(i);// add current row value
								hasE[i] = true;
							}
						}
						if (!hasE[i] && saveNoNoise) {
							newE[newENum++] = i;
						}
						continue;
					}
					else {
						currentPos = max(currentPos + 1, beginPos);
						lastMainLabelPos = lastMainLabelT - startT;
					}
				}
			}

			forbidTimeStartT = -1;
			auto tail = now.tail;
			if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
				forbidTimeStartT = tail->item.first;
				if (tail->pre != nullptr)
					now.deleteNextNode(tail->pre);
				else now.deleteFirstNode();
			}
		}
		else {
			localNoise = 0;
			currentPos = lastMainLabelPos = beginPos;
			forbidTimeStartT = -1;
			noiseNum = 0;
		}
		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;


		lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
				//if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos - 1, tempTimestampPosForIT - nEdge, i);//update aft
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
				}
			}
			else {
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
				}
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;


		if (lastMainLabelPos < endPos) {
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			if (saveNoNoise) {
				newE[newENum++] = i;
			}
			continue;
		}

		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			selectedNum++;
			savePos = lastMainLabelPos - endPos;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
			hasE[i] = true;
		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= filterE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}

			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
			hasE[i] = true;
		}
		else {
			//dynamic
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= filterE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}
			
			selectedNum++;
			savePos = currentT - intvE;
			if (rightEndpoint < savePos) rightEndpoint = savePos;
			edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
			edgeSetsRAdd[savePos].emplace_back(i);// add current row value
			hasE[i] = true;
		}
		if (!hasE[i] && saveNoNoise) {
			newE[newENum++] = i;
		}
	}

}

void TGraphUDEL::edgeFilterPlusMidRForDYN(int intvB, int intvE, int filterE, int limited, int oriEndTE, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, vec(int)*& edgeSetsRAdd,
	vec(int)*& edgeSetsR, int& selectedNum, vec(int)*& edgeSetsRShortIntv, int& selectedNumShortIntv, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) {

	int edgeType, mainLabel;
	int beginPos = intvB - startT;
	int endPos = intvE - startT, timePos = filterE - startT, checkEndPos = oriEndTE - startT;
	int intvLen, intvLenNoNoise = filterE - intvB + 1;
	int currentPos, mainLabelPos, currentT;
	int labelsSum;
	int savePos;
	int graphLastPos = endT - startT;
	int allLen = graphLastPos - beginPos + 1;
	int nextLabPos, lastMainLabelPos, lastMainLabelT;
	int forbidTimeStartT;
	int noiseNum;
	int checkPos = beginPos - 1;
	int labelPosForEdge;
	int EMaxIntvlStartT, EMaxIntvlEndT;
	int localNoise;
	int noiseT;
	double minNoise, checkNum;
	bool saveNoNoise;

	int timestampPosForEMaxIntvlEndT;
	auto fromMidREnd = edgesInEIntR->end(), fromMidRIter = edgesInEIntR->begin();
	for (; fromMidRIter != fromMidREnd; ++fromMidRIter) {

		int i = *fromMidRIter;
		
		labelPosForEdge = i * numOfLabel;

		mainLabel = lab[posUsedForEdgeFilter + i];//mainL = L^intvB(e)
		if (isEdgeTypeFixed && !fixLabel[mainLabel]) continue;

		mainLabelPos = labelPosForEdge + mainLabel;//position of <e,mainL>

		auto& now = vioT[mainLabelPos];//tabuT[e,mainL] before updated
		currentPos = scanT[mainLabelPos];//scanT[e,mainL] before updated
		CircularQueue<NVIntv>& preEMaxIntvlPtr = preMaxIntv[i];//maxIntv[e,mainL] 
		timestampPosForEMaxIntvlEndT = (maxIntv[i].second - startT) * nEdge + i;
		if (posInEIntR[mainLabelPos] == EMaxIntvlChange::INTVINIT) {//initial 
			continue;
		}
		else {
			saveNoNoise = false;
			edgeType = lab[TGraph::posUsedForEdgeFilterShortIntv + i];
			/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
				put here for less time*/
			if (max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) >= intvLenNoNoise) {
				int intvStart = filterE - max(bef[TGraph::posUsedForEdgeFilterShortIntv + i], 1) + 1;
				int j = timePos + 1;
				if (maxIntvShortIntv[i].second >= timePos && maxIntvShortIntv[i].first <= beginPos) {
					if (maxIntvShortIntv[i].second >= oriEndTE) {
						int check = min(maxIntvShortIntv[i].second - oriEndTE, limited);
						selectedNumShortIntv++;
						edgeSetsRShortIntv[check].emplace_back(i);
						saveNoNoise = true;
					}
				}
				else if (maxIntvShortIntv[i].second >= beginPos && maxIntvShortIntv[i].first <= timePos) {
				}
				else {//case 1 and 2
					maxIntvShortIntv[i].first = intvStart;
					selectedNumShortIntv++;

					lazyUpdate(filterE, TGraph::posUsedForEdgeFilterShortIntv + i, i);//update aft
					int check = max(aft[TGraph::posUsedForEdgeFilterShortIntv + i], 1) - 1;
					maxIntvShortIntv[i].second = timePos + check + startT;

					if (maxIntvShortIntv[i].second >= oriEndTE) {
						check = min(maxIntvShortIntv[i].second - oriEndTE, limited);
						edgeSetsRShortIntv[check].emplace_back(i);
						saveNoNoise = true;
					}
				}
			}
			
			if (posInEIntR[mainLabelPos] >= 0) {//continue to scan

				if (currentPos == 0) continue; //initial

				//recover localNoise , noiseNum and lastMainLabelPos
				int timestampPosForEL = currentPos * nEdge + i;
				edgeType = lab[timestampPosForEL];
				if (edgeType != mainLabel) {
					int eLvalue = max(bef[timestampPosForEL], 1);
					int tempTPos = timestampPosForEL - eLvalue * nEdge;
					int tempPosForLabels = currentPos - eLvalue;
					localNoise = eLvalue;
					while (lab[tempTPos] != mainLabel) {
						eLvalue = max(bef[tempTPos], 1);
						tempPosForLabels -= eLvalue;
						localNoise += eLvalue;
						tempTPos -= eLvalue * nEdge;
					}
					noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
				}
				else {
					localNoise = 0;
					noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
				}

				if (maxIntv[i].second == -1) {
					lastMainLabelPos = beginPos;
				}
				else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
					lastMainLabelPos = max(maxIntv[i].second, beginPos);
				}
				else {
					EMaxIntvlEndT = -1;
					int last = preEMaxIntvlPtr.rear;
					if (last != preEMaxIntvlPtr.front) {
						int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
						for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
							EMaxIntvlEndT = preEMaxIntvlPtr.q[temp].second;
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
						}
						if (temp == preEMaxIntvlPtr.front) {
							EMaxIntvlEndT = preEMaxIntvlPtr.q[temp].second;
							if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
						}
					}
					lastMainLabelPos = max(EMaxIntvlEndT, beginPos);
				}

				//recover vioT
				forbidTimeStartT = -1;
				auto tail = now.tail;
				if (tail != nullptr && tail->item.second == currentPos + startT) {
					forbidTimeStartT = tail->item.first;
					lastMainLabelPos = forbidTimeStartT - 1;
					if (tail->pre != nullptr)
						now.deleteNextNode(tail->pre);
					else now.deleteFirstNode();
				}
				currentPos = max(currentPos + 1, beginPos);
			}
			else {
				if (currentPos != 0) {

					if (mainLabel == lab[posUsedForEdgeFilter - nEdge + i]) {// L^intvB(e) = L^(intvB-1)(e)
						if (maxIntv[i].second == -1 || maxIntv[i].second < intvE || lab[timestampPosForEMaxIntvlEndT] != mainLabel) {//no valid R sets for e
							if (saveNoNoise) {
								newE[newENum++] = i; 
								posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R for edge i is changed
							}
							else {
								if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
							}
							continue;
						}
						else {
							auto intvItem = now.first;

							bool shrink = false;
							while (intvItem != nullptr) {//check tabuT[e,mainL]

								if (intvItem->item.first > maxIntv[i].second) break;//no noise
								else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
									shrink = true;
									int tempEndT = intvItem->item.first - 1;
									savePos = tempEndT - intvE;
									if (savePos >= 0) {
										for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
											temp++) {
											auto& intv = preEMaxIntvlPtr.q[temp];
											if (intv.second < filterE) {//no overlap
												preEMaxIntvlPtr.swapToTop(temp);
												preEMaxIntvlPtr.pop();
											}
											//else break;
										}
										savePos = tempEndT - oriEndTE;
										if (savePos >= 0) {
											selectedNum++;
											if (rightEndpoint < savePos) rightEndpoint = savePos;
											edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
											hasE[i] = true;
										}
									}
									break;
								}

								int tempEWPos = (intvItem->item.second - startT)*nEdge + i;
								if (lab[tempEWPos] == mainLabel) {
									intvItem->item.second++;
								}
								else {
									int labelSum = intvItem->item.second - intvB + 2;
									if (MORE(dif[tempEWPos + nEdge] - dif[posUsedForEdgeFilter + i], Setting::delta * labelSum)) {
										intvItem->item.second++;
									}
								}

								if (intvItem->next != nullptr) {
									if (intvItem->item.second + 1 == intvItem->next->item.first) {//combine
										intvItem->item.second = intvItem->next->item.second;
										now.deleteNextNode(intvItem);
									}
									else intvItem = intvItem->next;
								}
								else if (intvItem->item.second >= maxIntv[i].second) {//R set to which e belongs shrinks
									shrink = true;
									int tempEndT = intvItem->item.first - 1;
									if (tempEndT >= intvE) {
										for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
											temp++) {//maintain maxIntv[e,mainL]
											auto& intv = preEMaxIntvlPtr.q[temp];
											if (intv.second < filterE) {//no overlap
												preEMaxIntvlPtr.swapToTop(temp);
												preEMaxIntvlPtr.pop();
											}
										}
										savePos = tempEndT - oriEndTE;
										if (savePos >= 0) {
											selectedNum++;
											if (rightEndpoint < savePos) rightEndpoint = savePos;
											edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
											hasE[i] = true;
										}
									}
									break;
								}
								else intvItem = intvItem->next;
							}

							if (!shrink) {
								if (maxIntv[i].second >= intvE) {//R set does not change
									for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear;
										temp++) {
										auto& intv = preEMaxIntvlPtr.q[temp];
										if (intv.second < filterE) {//no overlap
											preEMaxIntvlPtr.swapToTop(temp);
											preEMaxIntvlPtr.pop();
										}
										//else break;
									}
									savePos = maxIntv[i].second - oriEndTE;
									if (savePos >= 0) {
										selectedNum++;
										if (rightEndpoint < savePos) rightEndpoint = savePos;
										edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
										hasE[i] = true;
									}
								}
							}
						}
						if (saveNoNoise) {
							if (!hasE[i]) { 
								newE[newENum++] = i; 
								posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R for edge i is changed
							}
						}
						else {
							if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
						}
						now.removeNodeLessThan(intvE);
						continue;
					}
					else {
						//Test::comprp2n++;
						now.removeNodeLessThan(intvE);

						if (maxIntv[i].second == -1) {
							EMaxIntvlEndT = -1;
						}
						else if (lab[timestampPosForEMaxIntvlEndT] == mainLabel) {
							EMaxIntvlEndT = maxIntv[i].second;
						}
						else {
							EMaxIntvlEndT = -1;
							int last = preEMaxIntvlPtr.rear;
							if (last != preEMaxIntvlPtr.front) {
								int temp = preEMaxIntvlPtr.rear - 1 /*+ CircularQueue<NVIntv>::queueSize) % CircularQueue<NVIntv>::queueSize*/;
								for (; temp != preEMaxIntvlPtr.front; temp--/*) % CircularQueue<NVIntv>::queueSize*/) {
									tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
									if (lab[(EMaxIntvlEndT - startT)*nEdge + i] == mainLabel) break;
								}
								if (temp == preEMaxIntvlPtr.front) {
									tie(EMaxIntvlStartT, EMaxIntvlEndT) = preEMaxIntvlPtr.q[temp];
									if (lab[(EMaxIntvlEndT - startT)*nEdge + i] != mainLabel) EMaxIntvlEndT = -1;
								}
							}
						}
						lastMainLabelT = max(EMaxIntvlEndT, intvB);
						if (currentPos < beginPos) { //does not scan
							//now.localNoise = 0;
							localNoise = 0;
							currentPos = lastMainLabelPos = beginPos;
							noiseNum = 0;
						}
						else {
							auto forbidIntv = now.first;
							if (forbidIntv == nullptr) {//does not need update
								if (currentPos != currNTimestamp) {
									int timestampPosForEL = currentPos * nEdge + i;
									edgeType = lab[timestampPosForEL];
									if (edgeType != mainLabel) {
										int eLvalue = max(bef[timestampPosForEL], 1);
										int tempTPos = timestampPosForEL - eLvalue * nEdge;
										int tempPosForLabels = currentPos - eLvalue;
										localNoise = eLvalue;
										while (lab[tempTPos] != mainLabel) {
											eLvalue = max(bef[tempTPos], 1);
											tempPosForLabels -= eLvalue;
											localNoise += eLvalue;
											tempTPos -= eLvalue * nEdge;
										}
										noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
									}
									else {
										localNoise = 0;
										noiseNum = dif[timestampPosForEL] - dif[posUsedForEdgeFilter + i];
									}
								}
							}
							else {
								//update vioT

								//get noiseNum, localNoise
								int tempT = forbidIntv->item.first, stopT = forbidIntv->item.second;
								int tempPos = tempT - startT, intvLen = tempT - intvB;
								labelsSum = tempT - intvB;
								int tempTimestampPosForEL = tempPos * nEdge + i;
								edgeType = lab[tempTimestampPosForEL];
								if (edgeType != mainLabel) {
									int eLvalue = max(bef[tempTimestampPosForEL], 1);
									int tempTPos = tempTimestampPosForEL - eLvalue * nEdge;
									int tempPosForLabels = tempPos - eLvalue;
									localNoise = eLvalue - 1;
									while (lab[tempTPos] != mainLabel) {
										eLvalue = max(bef[tempTPos], 1);
										tempPosForLabels -= eLvalue;
										localNoise += eLvalue;
										tempTPos -= eLvalue * nEdge;
									}
									noiseNum = dif[tempTPos] - dif[posUsedForEdgeFilter + i] + localNoise;
								}
								else {
									noiseNum = dif[tempTimestampPosForEL] - dif[posUsedForEdgeFilter + i];
									localNoise = 0;
								}

								auto tempIntv = forbidIntv;
								while (forbidIntv != nullptr) {

									forbidTimeStartT = -1;
									lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
									int nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
									while (tempT <= stopT) {

										intvLen = nextLabT - tempT;
										edgeType = lab[tempTimestampPosForEL];
										if (edgeType == mainLabel) {
											localNoise = 0;
											minNoise = Setting::delta * (labelsSum + 1);
											if (LESSEQ(noiseNum, minNoise)) {
												//lastMainLabelPos = nextLabPos - 1;
												//remove = 0;
												if (forbidTimeStartT != -1) {
													now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
													forbidTimeStartT = -1;
													tempIntv = tempIntv->next;
												}
											}
											else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
												//remove += intvLen;
												if (forbidTimeStartT == -1) {
													forbidTimeStartT = tempT + startT;
												}
											}
											else {
												//remove = 0;
												//lastMainLabelPos = nextLabPos - 1;
												checkNum = (noiseNum - minNoise) / Setting::delta;
												if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + tempT;
												else noiseT = (int)checkNum + tempT;
												if (forbidTimeStartT == -1) {
													now.addNodeAft(tempIntv, make_pair(tempT, min(noiseT, stopT)));
													tempIntv = tempIntv->next;
												}
												else {
													now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, min(noiseT, stopT)));
													tempIntv = tempIntv->next;
													forbidTimeStartT = -1;
												}
											}
										}
										else {
											if (forbidTimeStartT == -1) {
												forbidTimeStartT = tempT + startT;
											}
											localNoise += intvLen;
											noiseNum += intvLen;
										}
										labelsSum += intvLen;

										tempT = nextLabT;
										tempTimestampPosForEL = (tempT - startT) * nEdge + i;
										if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
											if (tempT <= stopT) {
												lazyUpdate(tempT - 1, tempTimestampPosForEL - nEdge, i);//update aft
												nextLabT = min(nextLabT - 2 - aft[tempTimestampPosForEL - nEdge], stopT) + 1;
											}
										}
										else {
											if (tempT <= stopT) {
												lazyUpdate(tempT, tempTimestampPosForEL, i);//update aft
												nextLabT = min(tempT + max(aft[tempTimestampPosForEL], 1) - 1, stopT) + 1;
											}
										}
									}
									if (forbidTimeStartT != -1) {
										if (forbidTimeStartT != forbidIntv->item.first || tempT - 1 != stopT) {
											now.addNodeAft(tempIntv, make_pair(forbidTimeStartT, tempT - 1));
											tempIntv = tempIntv->next;

											forbidIntv = forbidIntv->pre;
											tempIntv = tempIntv->next;
											if (forbidIntv == nullptr) {
												now.deleteFirstNode();
											}
											else now.deleteNextNode(forbidIntv);
											forbidIntv = tempIntv;

										}
										else {
											tempIntv = forbidIntv = forbidIntv->next;
										}
									}
									else {
										forbidIntv = forbidIntv->pre;
										tempIntv = tempIntv->next;
										if (forbidIntv == nullptr) {
											now.deleteFirstNode();
										}
										else now.deleteNextNode(forbidIntv);
										forbidIntv = tempIntv;
									}

									if (forbidIntv != nullptr) {
										tempT = forbidIntv->item.first;
										tempTimestampPosForEL = (tempT - startT) * nEdge + i;
										labelsSum = tempT - intvB;
										localNoise = 0;
										stopT = forbidIntv->item.second;
									}
								}

								//get the last time where the edge has main label
								auto lastIntv = now.tail;
								int updateEndT = min(currentPos + startT, endT);
								if (lastIntv != nullptr) {
									//tie(intvStartT, intvEndT/*, std::ignore*/) = lastIntv->item;
									if (lastIntv->item.second != updateEndT) {
										lastMainLabelT = updateEndT;
										localNoise = 0;
									}
									else {
										lastMainLabelT = lastIntv->item.first - 1;
									}
								}
								else {
									lastMainLabelT = updateEndT;
									localNoise = 0;
								}
							}
							if (currentPos + 1 == currNTimestamp /*|| MORE(noiseNum, Setting::delta * allLen)*/ || localNoise > Setting::c) {
								int checkPos = (lastMainLabelT - startT)*nEdge + i;
								lazyUpdate(lastMainLabelT, checkPos, i);//update aft
								int labelEnd = lastMainLabelT - startT + max(aft[checkPos], 1) - 1;
								//dynamic
								if (localNoise <= Setting::c) {
									checkPos = labelEnd * nEdge + i;
									lazyUpdate(labelEnd, checkPos, i);//update aft
									if (aft[checkPos] != -MYINFINITE) {
										if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
											newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
											MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
											newEIntR->emplace_back(rd4);
											newValidMidResult.emplace_back(true);
										}
									}
									else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
										newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
										MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
										newEIntR->emplace_back(rd4);
										newValidMidResult.emplace_back(true);
									}
								}
								else {
									if (newPosInEIntR[mainLabelPos] >= 0) {
										newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
										newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
									}
								}


								if (lastMainLabelT >= intvE) {
									if (lastMainLabelT >= oriEndTE)
										posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R# for edge i is changed
									else if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R# for edge i is unchanged 

									if (lastMainLabelT == maxIntv[i].second) {
										for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
											auto& intv = preEMaxIntvlPtr.q[temp];
											if (intv.second < filterE) {//no overlap
												preEMaxIntvlPtr.swapToTop(temp);
												preEMaxIntvlPtr.pop();
											}
										}
										savePos = lastMainLabelT - oriEndTE;
										if (savePos >= 0) {
											selectedNum++;
											if (rightEndpoint < savePos) rightEndpoint = savePos;
											edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
											edgeSetsRAdd[savePos].emplace_back(i);// add current row value
											hasE[i] = true;
										}
									}
									else if (lastMainLabelT < maxIntv[i].second) {
										int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
										for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
											auto& intv = preEMaxIntvlPtr.q[temp];
											if (intv.second < filterE) {//no overlap
												preEMaxIntvlPtr.swapToTop(temp);
												preEMaxIntvlPtr.pop();
											}
											else {
												int tempTimestampPos = (intv.second - startT)*nEdge + i;
												if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
													preMaxEMaxIntvlEndT = intv.second;
													preMaxEMaxIntvlStartT = intv.first;
												}
											}
										}

										if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
											if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= filterE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
											if (EMaxIntvlEndT < lastMainLabelT) {
												maxIntv[i].first = intvB;
												maxIntv[i].second = lastMainLabelT;
											}
											else {
												maxIntv[i].first = EMaxIntvlStartT;
												maxIntv[i].second = EMaxIntvlEndT;
											}
										}

										savePos = lastMainLabelT - oriEndTE;
										if (savePos >= 0) {
											selectedNum++;
											if (rightEndpoint < savePos) rightEndpoint = savePos;
											edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
											edgeSetsRAdd[savePos].emplace_back(i);// add current row value
											hasE[i] = true;
										}
									}
									else {
										int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
										for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
											auto& intv = preEMaxIntvlPtr.q[temp];
											if (intv.second < filterE) {//no overlap
												preEMaxIntvlPtr.swapToTop(temp);
												preEMaxIntvlPtr.pop();
											}
											else {
												int tempTimestampPos = (intv.second - startT)*nEdge + i;

												if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
													preMaxEMaxIntvlEndT = intv.second;
													preMaxEMaxIntvlStartT = intv.first;
												}
											}
										}
										if (maxIntv[i].second >= filterE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
											preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
										}

										if (EMaxIntvlEndT < lastMainLabelT) {
											maxIntv[i].first = intvB;
											maxIntv[i].second = lastMainLabelT;
										}
										else {
											maxIntv[i].first = EMaxIntvlStartT;
											maxIntv[i].second = EMaxIntvlEndT;
										}

										savePos = lastMainLabelT - oriEndTE;
										if (savePos >= 0) {
											selectedNum++;
											if (rightEndpoint < savePos) rightEndpoint = savePos;
											edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
											edgeSetsRAdd[savePos].emplace_back(i);// add current row value
											hasE[i] = true;
										}
									}
								}
								if (saveNoNoise) {
									if (!hasE[i]) { 
										newE[newENum++] = i;
										posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R for edge i is changed
									}
								}
								else {
									if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
								}
								continue;
							}
							else {
								currentPos = max(currentPos + 1, beginPos);
								lastMainLabelPos = lastMainLabelT - startT;
							}
						}


						//Test::comprp2t += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - beginTest).count();
					}

					forbidTimeStartT = -1;
					auto tail = now.tail;
					if (tail != nullptr && tail->item.second + 1 == currentPos + startT) {
						forbidTimeStartT = tail->item.first;
						if (tail->pre != nullptr)
							now.deleteNextNode(tail->pre);
						else now.deleteFirstNode();
					}
				}
				else {
					localNoise = 0;
					currentPos = lastMainLabelPos = beginPos;
					forbidTimeStartT = -1;
					noiseNum = 0;
				}
			}
		}
		labelsSum = currentPos - beginPos;

		int tempTimestampPosForIT = currentPos * nEdge + i;


		lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
		nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
		while (currentPos < currNTimestamp) {
			intvLen = nextLabPos - currentPos;
			edgeType = lab[tempTimestampPosForIT];
			if (edgeType == mainLabel) {
				localNoise = 0;
				minNoise = Setting::delta * (labelsSum + 1);
				if (LESSEQ(noiseNum, minNoise)) {
					lastMainLabelPos = nextLabPos - 1;
					//remove = 0;
					if (forbidTimeStartT != -1) {
						now.addItemAtLast(make_pair(forbidTimeStartT, currentPos + startT - 1));
						forbidTimeStartT = -1;
					}
				}
				else if (MORE(noiseNum, Setting::delta * (labelsSum + intvLen))) {
					if (forbidTimeStartT == -1) {
						forbidTimeStartT = max(endPos, currentPos) + startT;
					}
				}
				else {
					lastMainLabelPos = nextLabPos - 1;
					checkNum = (noiseNum - minNoise) / Setting::delta;
					if (EQ(round(checkNum), checkNum)) noiseT = (int)round(checkNum) - 1 + currentPos + startT;
					else noiseT = (int)checkNum + currentPos + startT;
					if (forbidTimeStartT == -1) {
						now.addItemAtLast(make_pair(max(endPos, currentPos) + startT, noiseT));
					}
					else {
						now.addItemAtLast(make_pair(forbidTimeStartT, noiseT));
						forbidTimeStartT = -1;
					}
				}
			}
			else {
				if (forbidTimeStartT == -1) {
					forbidTimeStartT = max(endPos, currentPos) + startT;
				}
				localNoise += intvLen;
				noiseNum += intvLen;
				if (localNoise > Setting::c) break;
			//	if (MORE(noiseNum, Setting::delta * allLen)) break;
			}
			labelsSum += intvLen;

			currentPos = nextLabPos;
			tempTimestampPosForIT = currentPos * nEdge + i;
			if (edgeType == mainLabel) {//Lab[currentPos] != mainLabel
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos - 1, tempTimestampPosForIT - nEdge, i);//update aft
					nextLabPos = min(nextLabPos - 2 - aft[tempTimestampPosForIT - nEdge], endT - startT) + 1;
				}
			}
			else {
				if (currentPos < currNTimestamp) {
					lazyUpdate(currentPos, tempTimestampPosForIT, i);//update aft
					nextLabPos = min(currentPos + max(aft[tempTimestampPosForIT], 1) - 1, endT - startT) + 1;
				}
			}
		}
		if (forbidTimeStartT != -1) {
			now.addItemAtLast(make_pair(forbidTimeStartT, min(nextLabPos + startT - 1, endT)));
		}
		scanT[mainLabelPos] = nextLabPos - 1;

		if (lastMainLabelPos >= checkEndPos)
			posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R# for edge i is changed
		else if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;//R# for edge i is unchanged 
		if (lastMainLabelPos < endPos) {
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			//dynamic
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}


			if (saveNoNoise) {
				newE[newENum++] = i;
				posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R for edge i is changed
			}
			else {
				if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
			}
			continue;
		}

		//if ((isEdgeTypeFixed && fixLabel.find(idToLabel[mainLabel[i]]) == fixLabelEnd)) continue;
		currentT = lastMainLabelPos + startT;
		if (currentT == maxIntv[i].second) {
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
			}
			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
				edgeSetsRAdd[savePos].emplace_back(i);// add current row value
				hasE[i] = true;
			}
		}
		else if (currentT < maxIntv[i].second) {
			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}

			if (lab[timestampPosForEMaxIntvlEndT] != mainLabel) {
				if (maxIntv[i].second > preMaxEMaxIntvlEndT && maxIntv[i].second >= filterE) preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
				if (maxEMaxIntvlEndT < currentT) {
					/*if (Test::testingMode == 3) {
						Test::containment++;
						if (maxIntv[i].second > preMaxEMaxIntvlEndT) {
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << maxIntv[i].first << "," << maxIntv[i].second << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
						else{
							cout << edgeList[i].first << "," << edgeList[i].second << ":[" << preMaxEMaxIntvlStartT << "," << preMaxEMaxIntvlEndT << "," << lab[timestampPosForEMaxIntvlEndT + i] << "] | [" << intvB << "," << currentT << "," << idToLabel[mainLabel] << "]" << endl;
						}
					}*/
					maxIntv[i].first = intvB;
					maxIntv[i].second = currentT;
				}
				else {
					maxIntv[i].first = maxEMaxIntvlStartT;
					maxIntv[i].second = maxEMaxIntvlEndT;
					if (maxEMaxIntvlEndT > currentT) {
						now.addItemAtLast(make_pair(currentT + 1, maxEMaxIntvlEndT));
					}
				}
			}

			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
				edgeSetsRAdd[savePos].emplace_back(i);// add current row value
				hasE[i] = true;
			}
		}
		else {
			//dynamic
			int checkPos = lastMainLabelPos * nEdge + i;
			lazyUpdate(lastMainLabelPos, checkPos, i);//update aft
			int labelEnd = lastMainLabelPos + max(aft[checkPos], 1) - 1;
			if (localNoise <= Setting::c) {
				checkPos = labelEnd * nEdge + i;
				lazyUpdate(labelEnd, checkPos, i);//update aft
				if (aft[checkPos] != -MYINFINITE) {
					if (-aft[checkPos] - 1 <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
						newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
						MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
						newEIntR->emplace_back(rd4);
						newValidMidResult.emplace_back(true);
					}
				}
				else if (endT - startT - labelEnd <= Setting::c && newPosInEIntR[mainLabelPos] < 0) {
					newPosInEIntR[mainLabelPos] = (int)newEIntR->size();
					MidResult* rd4 = DBG_NEW MidResult(intvB, i/*, lastMainLabelPos*/, maxIntv[i], preMaxIntv[i], scanT[mainLabelPos], now);
					newEIntR->emplace_back(rd4);
					newValidMidResult.emplace_back(true);
				}
			}
			else {
				if (newPosInEIntR[mainLabelPos] >= 0) {
					newValidMidResult[newPosInEIntR[mainLabelPos]] = false;
					newPosInEIntR[mainLabelPos] = EMaxIntvlChange::INTVINIT;
				}
			}


			int maxEMaxIntvlEndT = -1, maxEMaxIntvlStartT;
			int preMaxEMaxIntvlEndT = -1, preMaxEMaxIntvlStartT;
			for (int temp = preEMaxIntvlPtr.front; temp != preEMaxIntvlPtr.rear; temp++) {
				auto& intv = preEMaxIntvlPtr.q[temp];
				if (intv.second < filterE) {//no overlap
					preEMaxIntvlPtr.swapToTop(temp);
					preEMaxIntvlPtr.pop();
				}
				else {
					int tempTimestampPos = (intv.second - startT)*nEdge + i;
					if (lab[tempTimestampPos] == mainLabel && intv.second > maxEMaxIntvlEndT) {
						maxEMaxIntvlEndT = intv.second;
						maxEMaxIntvlStartT = intv.first;
					}
					if (lab[tempTimestampPos] == lab[timestampPosForEMaxIntvlEndT] && intv.second > preMaxEMaxIntvlEndT) {
						preMaxEMaxIntvlEndT = intv.second;
						preMaxEMaxIntvlStartT = intv.first;
					}
				}
			}
			if (maxIntv[i].second >= filterE && maxIntv[i].second > preMaxEMaxIntvlEndT) {
				preEMaxIntvlPtr.push(make_pair(maxIntv[i].first, maxIntv[i].second));
			}

			if (maxEMaxIntvlEndT < lastMainLabelPos + startT) {
				maxIntv[i].first = intvB;
				maxIntv[i].second = lastMainLabelPos + startT;
			}
			else {
				maxIntv[i].first = maxEMaxIntvlStartT;
				maxIntv[i].second = maxEMaxIntvlEndT;
				if (maxEMaxIntvlEndT > lastMainLabelPos + startT) {
					now.addItemAtLast(make_pair(lastMainLabelPos + startT + 1, maxEMaxIntvlEndT));
				}
			}

			savePos = currentT - oriEndTE;
			if (savePos >= 0) {
				selectedNum++;
				if (rightEndpoint < savePos) rightEndpoint = savePos;
				edgeSetsR[savePos].emplace_back(i/*, mainLabel*/);
				edgeSetsRAdd[savePos].emplace_back(i);// add current row value
				hasE[i] = true;
			}
		}

		if (saveNoNoise) {
			if (!hasE[i]) {
				newE[newENum++] = i;
				posInEIntR[mainLabelPos] = EMaxIntvlChange::CHANGED;//R for edge i is changed
			}
		}
		else {
			if (posInEIntR[mainLabelPos] != -1) posInEIntR[mainLabelPos] = EMaxIntvlChange::UNCHANGED;
		}
	}
}

#pragma endregion


#pragma region expCheck
bool TGraphUDEL::bothFitDefAndSameLabelChangeStartTimePos(vec(int)& edges, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) {
	auto edgeEnd = edges.end();
	int edgeId;

	int sum, field = startTimePos2 - startTimePos1 + 1;
	int label;
	double minNoise, noiseLimit, checkNum;
	int next, p, noiseT, noiseNum/*, timestampPosForLabels*/;
	for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter) {
		edgeId = *edgeIter;
		//timestampPosForLabels = edgeId * allNTimestamp;
		label = lab[mainLabelPos*nEdge + edgeId];
		int tempPosForEndTime = endTimePos * nEdge + edgeId;
		if (label != lab[tempPosForEndTime]) {
			//cout << mainLabelPos << " " << endTimePos << " " << edgeId << "%%" << endl;
			return false;
		}
		p = startTimePos2;
		int tempstampPosForEL = p * nEdge + edgeId;
		next = p - max(bef[tempstampPosForEL], 1);
		while (p >= startTimePos1) {
			if (label != lab[tempstampPosForEL]) {
				int stopP = max(next + 1, startTimePos1);
				for (int tabu = p; tabu >= stopP; tabu--) {
					if (expandMask[tabu - startTimePos1]) {
						field--;
						expandMask[tabu - startTimePos1] = false;
						if (field == 0)break;
					}
				}
			}
			else {
				sum = endTimePos - p;
				noiseNum = dif[tempPosForEndTime] - dif[tempstampPosForEL];
				noiseLimit = Setting::delta * sum;

				if (MORE(noiseNum, noiseLimit)) {
					minNoise = noiseLimit + Setting::delta;
					if (noiseNum >= minNoise) {
						checkNum = (noiseNum - minNoise) / Setting::delta;
						if (EQ(round(checkNum), checkNum)) noiseT = p - ((int)round(checkNum) - 1);
						else noiseT = p - (int)checkNum;


						for (int tabu = p; tabu >= startTimePos1 && tabu >= noiseT; tabu--) {
							if (expandMask[tabu - startTimePos1]) {
								field--;
								expandMask[tabu - startTimePos1] = false;
								if (field == 0)break;
							}
						}
					}
				}
			}
			if (field == 0)break;
			p = next;
			tempstampPosForEL = p * nEdge + edgeId;
			if (p >= 0) next = p - max(bef[tempstampPosForEL],1);
		}
		//cout << field << " " << edgeId << endl;
		if (field == 0)break;
	}
	return field != 0;
}
bool TGraphUDEL::bothFitDefAndSameLabelChangeStartTimePos(vec(int)& edges, int*&subCCId, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) {
	auto edgeEnd = edges.end();
	auto subCCIter = &subCCId[0];
	int edgeId;

	int sum, field = startTimePos2 - startTimePos1 + 1;
	int label;
	double minNoise, noiseLimit, checkNum;
	int next, p, noiseT, noiseNum/*, timestampPosForLabels*/;
	for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter,++subCCIter) {
		if (*subCCIter != -1) {
			edgeId = *edgeIter;
			//timestampPosForLabels = edgeId * allNTimestamp;
			label = lab[mainLabelPos*nEdge + edgeId];
			int tempPosForEndTime = endTimePos * nEdge + edgeId;
			if (label != lab[tempPosForEndTime]) {
				return false;
			}
			p = startTimePos2;
			int tempstampPosForEL = p * nEdge + edgeId;
			next = p - max(bef[tempstampPosForEL], 1);
			while (p >= startTimePos1) {
				if (label != lab[tempstampPosForEL]) {
					int stopP = max(next + 1, startTimePos1);
					for (int tabu = p; tabu >= stopP; tabu--) {
						if (expandMask[tabu - startTimePos1]) {
							field--;
							expandMask[tabu - startTimePos1] = false;
							if (field == 0)break;
						}
					}
				}
				else {
					sum = endTimePos - p;
					noiseNum = dif[tempPosForEndTime] - dif[tempstampPosForEL];
					noiseLimit = Setting::delta * sum;

					if (MORE(noiseNum, noiseLimit)) {
						minNoise = noiseLimit + Setting::delta;
						if (noiseNum >= minNoise) {
							checkNum = (noiseNum - minNoise) / Setting::delta;
							if (EQ(round(checkNum), checkNum)) noiseT = p - ((int)round(checkNum) - 1);
							else noiseT = p - (int)checkNum;


							for (int tabu = p; tabu >= startTimePos1 && tabu >= noiseT; tabu--) {
								if (expandMask[tabu - startTimePos1]) {
									field--;
									expandMask[tabu - startTimePos1] = false;
									if (field == 0)break;
								}
							}
						}
					}
				}
				if (field == 0)break;
				p = next;
				tempstampPosForEL = p * nEdge + edgeId;
				if (p >= 0) next = p - max(bef[tempstampPosForEL], 1);
			}
			//cout << field << " " << edgeId << endl;
			if (field == 0)break;
		}
	}
	return field != 0;
}
bool TGraphUDEL::bothFitDefAndSameLabelChangeStartTimePos(vec(int)& subCCs, vec(int)& edges, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) {
	auto edgeEnd = edges.end();
	auto subCCEnd = subCCs.end();
	int edgeId;

	int sum, field = startTimePos2 - startTimePos1 + 1;
	int label;
	double minNoise, noiseLimit, checkNum;
	int next, p, noiseT, noiseNum/*, timestampPosForLabels*/;
	for (auto subCCIter = subCCs.begin(); subCCIter != subCCEnd; ++subCCIter) {
		edgeId = edges[*subCCIter];
		//timestampPosForLabels = edgeId * allNTimestamp;
		label = lab[mainLabelPos*nEdge + edgeId];
		int tempPosForEndTime = endTimePos * nEdge + edgeId;
		if (label != lab[tempPosForEndTime]) {
			return false;
		}
		p = startTimePos2;
		int tempstampPosForEL = p * nEdge + edgeId;
		next = p - max(bef[tempstampPosForEL],1);
		while (p >= startTimePos1) {
			if (label != lab[tempstampPosForEL]) {
				int stopP = max(next + 1, startTimePos1);
				for (int tabu = p; tabu >= stopP; tabu--) {
					if (expandMask[tabu - startTimePos1]) {
						field--;
						expandMask[tabu - startTimePos1] = false;
						if (field == 0)break;
					}
				}
			}
			else {
				sum = endTimePos - p;
				noiseNum = dif[tempPosForEndTime] - dif[tempstampPosForEL];
				noiseLimit = Setting::delta * sum;

				if (MORE(noiseNum, noiseLimit)) {
					minNoise = noiseLimit + Setting::delta;
					if (noiseNum >= minNoise) {
						checkNum = (noiseNum - minNoise) / Setting::delta;
						if (EQ(round(checkNum), checkNum)) noiseT = p - ((int)round(checkNum) - 1);
						else noiseT = p - (int)checkNum;


						for (int tabu = p; tabu >= startTimePos1 && tabu >= noiseT; tabu--) {
							if (expandMask[tabu - startTimePos1]) {
								field--;
								expandMask[tabu - startTimePos1] = false;
								if (field == 0)break;
							}
						}
					}
				}
			}
			if (field == 0)break;
			p = next;
			tempstampPosForEL = p * nEdge + edgeId;
			if (p >= 0) next = p - max(bef[tempstampPosForEL],1);
		}
		//cout << field << " " << edgeId << endl;
		if (field == 0)break;
	}
	return field != 0;
}
bool TGraphUDEL::bothFitDefAndSameLabelChangeEndTimePos(vec(int)& subCCs, vec(int)& edges, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) {
	auto edgeEnd = edges.end();
	auto subCCEnd = subCCs.end();
	int edgeId;

	int sum, field = endTimePos2 - endTimePos1 + 1;
	int label;
	double minNoise, noiseLimit, checkNum;
	int next, p, noiseT, noiseNum/*, timestampPosForLabels*/;
	for (auto subCCIter = subCCs.begin(); subCCIter != subCCEnd; ++subCCIter) {
		edgeId = edges[*subCCIter];
		//timestampPosForLabels = edgeId * allNTimestamp;
		label = lab[mainLabelPos*nEdge + edgeId];
		int tempPosForStartTime = startTimePos * nEdge + edgeId;
		if (label != lab[tempPosForStartTime]) {
			//cout << mainLabelPos << " "<< startTimePos<<" "<<edgeId << "%"<< endl;
			return false;
		}
		p = endTimePos2;
		int tempstampPosForEL = p * nEdge + edgeId;
		next = p - max(bef[tempstampPosForEL],1);
		while (p >= endTimePos1) {
			
			if (label != lab[tempstampPosForEL]) {
				int stopP = max(next + 1, endTimePos1);
				for (int tabu = p; tabu >= stopP; tabu--) {
					if (expandMask[tabu - endTimePos1]) {
						field--;
						expandMask[tabu - endTimePos1] = false;
						if (field == 0)break;
					}
				}
			}
			else {
				
				sum = p - startTimePos + 1;
				noiseNum = dif[tempstampPosForEL] - dif[tempPosForStartTime] ;
				noiseLimit = Setting::delta * (sum - max(bef[tempstampPosForEL],1));
				
				if (MORE(noiseNum, noiseLimit)) {
					minNoise = noiseLimit + Setting::delta;
					if (noiseNum >= minNoise) {
						checkNum = (noiseNum - minNoise) / Setting::delta;
						if (EQ(round(checkNum), checkNum)) noiseT = next + (int)round(checkNum);
						else noiseT = next + (int)checkNum + 1;

						for (int tabu = min(noiseT, endTimePos2); tabu >= endTimePos1 && tabu >= next; tabu--) {
							if (expandMask[tabu - endTimePos1]) {
								field--;
								expandMask[tabu - endTimePos1] = false;
								if (field == 0)break;
							}
						}
					}
				}
			}
			if (field == 0)break;
			p = next;
			tempstampPosForEL = p * nEdge + edgeId;
			if (p >= 0) next = p - max(bef[tempstampPosForEL],1);
		}
		//if(Test::counter == -1) cout << field << " " << edgeId << endl;
		if (field == 0)break;
	}
	//if (Test::counter == -1) cout << field << endl;
	return field != 0;
}
bool TGraphUDEL::bothFitDefAndSameLabelChangeEndTimePos(vec(int)& edges, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) {
	auto edgeEnd = edges.end();
	int edgeId;

	int sum, field = endTimePos2 - endTimePos1 + 1;
	int label;
	double minNoise, noiseLimit, checkNum;
	int next, p, noiseT, noiseNum/*, timestampPosForLabels*/;
	for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter) {
		edgeId = *edgeIter;
		//timestampPosForLabels = edgeId * allNTimestamp;
		label = lab[mainLabelPos*nEdge + edgeId];
		int tempPosForStartTime = startTimePos * nEdge + edgeId;
		if (label != lab[tempPosForStartTime]) {
			return false;
		}
		p = endTimePos2;
		int tempstampPosForEL = p * nEdge + edgeId;
		next = p - max(bef[tempstampPosForEL],1);
		while (p >= endTimePos1) {
			if (label != lab[tempstampPosForEL]) {
				
				int stopP = max(next + 1, endTimePos1);
				
				for (int tabu = p; tabu >= stopP; tabu--) {
					if (expandMask[tabu - endTimePos1]) {
						field--;
						/*if(Test::counter == -1 && tabu == 26273) {
							cout << startTimePos << " "<<edgeId << "%%" << endl;
						}*/
						expandMask[tabu - endTimePos1] = false;
						if (field == 0)break;
					}
				}
			}
			else {
				sum = p - startTimePos + 1;
				noiseNum = dif[tempstampPosForEL] - dif[tempPosForStartTime];
				noiseLimit = Setting::delta * (sum - max(bef[tempstampPosForEL],1));

				if (MORE(noiseNum, noiseLimit)) {
					minNoise = noiseLimit + Setting::delta;
					if (noiseNum >= minNoise) {
						checkNum = (noiseNum - minNoise) / Setting::delta;
						if (EQ(round(checkNum), checkNum)) noiseT = next + (int)round(checkNum);
						else noiseT = next + (int)checkNum + 1;

						for (int tabu = min(noiseT, endTimePos2); tabu >= endTimePos1 && tabu >= next; tabu--) {
							if (expandMask[tabu - endTimePos1]) {
								field--;
								/*if (Test::counter == -1 && tabu == 26273) {
									cout << startTimePos << " " << edgeId <<"%%"<<endl;
								}*/
								expandMask[tabu - endTimePos1] = false;
								if (field == 0)break;
							}
						}
					}
				}
			}
			if (field == 0)break;
			p = next;
			tempstampPosForEL = p * nEdge + edgeId;
			if (p >= 0) next = p - max(bef[tempstampPosForEL],1);
		}
		//cout << field << " " << edgeId << endl;
		if (field == 0)break;

	}

	return field != 0;
}
bool TGraphUDEL::bothFitDefAndSameLabelChangeEndTimePos(vec(int)& edges, int*&subCCId, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) {
	auto edgeEnd = edges.end();
	auto subCCIter = &subCCId[0];
	int edgeId;
	
	int sum, field = endTimePos2 - endTimePos1 + 1;
	int label;
	double minNoise, noiseLimit, checkNum;
	int next, p, noiseT, noiseNum/*, timestampPosForLabels*/;
	for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter, ++subCCIter) {
		if (*subCCIter != -1) {
			edgeId = *edgeIter;
			//timestampPosForLabels = edgeId * allNTimestamp;
			label = lab[mainLabelPos*nEdge + edgeId];
			int tempPosForStartTime = startTimePos * nEdge + edgeId;
			if (label != lab[tempPosForStartTime]) {
				return false;
			}
			p = endTimePos2;
			int tempstampPosForEL = p * nEdge + edgeId;
			next = p - max(bef[tempstampPosForEL], 1);
			while (p >= endTimePos1) {
				if (label != lab[tempstampPosForEL]) {

					int stopP = max(next + 1, endTimePos1);

					for (int tabu = p; tabu >= stopP; tabu--) {
						if (expandMask[tabu - endTimePos1]) {
							field--;
							expandMask[tabu - endTimePos1] = false;
							if (field == 0)break;
						}
					}
				}
				else {
					sum = p - startTimePos + 1;
					noiseNum = dif[tempstampPosForEL] - dif[tempPosForStartTime];
					noiseLimit = Setting::delta * (sum - max(bef[tempstampPosForEL], 1));

					if (MORE(noiseNum, noiseLimit)) {
						minNoise = noiseLimit + Setting::delta;
						if (noiseNum >= minNoise) {
							checkNum = (noiseNum - minNoise) / Setting::delta;
							if (EQ(round(checkNum), checkNum)) noiseT = next + (int)round(checkNum);
							else noiseT = next + (int)checkNum + 1;

							for (int tabu = min(noiseT, endTimePos2); tabu >= endTimePos1 && tabu >= next; tabu--) {
								if (expandMask[tabu - endTimePos1]) {
									field--;
									expandMask[tabu - endTimePos1] = false;
									if (field == 0)break;
								}
							}
						}
					}
				}
				if (field == 0)break;
				p = next;
				tempstampPosForEL = p * nEdge + edgeId;
				if (p >= 0) next = p - max(bef[tempstampPosForEL], 1);
			}
			if (field == 0)break;
		}
	}

	return field != 0;
}
#pragma endregion


/*test whether a motif have all edges with same label in endpoints of interval*/
void TGraphUDEL::checkMotifEndpointsD5(TMotifII*& motif) {
	vec(int)* motifEdge = motif->getMotifEdge();

	veciter(int) listEnd = motifEdge->end();
	int startP = motif->getStartT() - startT;
	int endP = motif->getEndT() - startT;
	int startLab, endLab;
	//vec(int)* maskEdge = motif->getMaskEdge();
	//bool finded;
	int noiseNum = 0, tempNoise, recordNoise;
	int maxLeft = 0, minRight = 0x7fffffff;
	for (auto iter = motifEdge->begin();
		iter != listEnd; ++iter) {
		//finded = Util::findItem(maskEdge, iter->id);

		//check endpoint
		startLab = getEdgeLabel(*iter, startP);
		endLab = getEdgeLabel(*iter, endP);
		if (startLab != endLab) {
			cout << "motif startT=" << motif->getStartT() << " motif endT=" << motif->getEndT() << endl;
			cout << "endpoint error: edge " << (*iter) << endl;
			return;
		}
		//check local noise and global noise
		noiseNum = 0;
		tempNoise = 0;
		for (int i = startP + 1; i < endP; ++i) {
			if (getEdgeLabel(*iter, i) != startLab) {
				tempNoise++;
				if (tempNoise > Setting::c) {
					cout << "motif startT=" << motif->getStartT() << " motif endT=" << motif->getEndT() << endl;
					cout << "local noise number error: edge " << (*iter) << endl;
					return;
				}
				noiseNum++;
			}
			else {
				tempNoise = 0;
			}
		}
		if (tempNoise > Setting::c) return;
		if (MORE(noiseNum, Setting::delta*(endP - startP + 1))) {
			cout << "motif startT=" << motif->getStartT() << " motif endT=" << motif->getEndT() << endl;
			cout << "global noise number error: edge " << (*iter) << endl;
			return;
		}

		tempNoise = 0;
		recordNoise = noiseNum;
		int p = startP - 1;
		int mainLabelP = startP;
		for (; p >= 0; --p) {
			if (getEdgeLabel(*iter, p) != startLab) {
				tempNoise++;
				if (tempNoise > Setting::c) break;
			}
			else {
				tempNoise = 0;
				mainLabelP = p;
			}
		}
		if (maxLeft < mainLabelP)maxLeft = mainLabelP;
		tempNoise = 0;
		recordNoise = noiseNum;
		mainLabelP = endP;
		p = endP + 1;
		for (; p < currNTimestamp; ++p) {
			if (getEdgeLabel(*iter, p) != startLab) {
				tempNoise++;
				if (tempNoise > Setting::c) break;
			}
			else {
				tempNoise = 0;
				mainLabelP = p;
			}
		}
		if (minRight > mainLabelP)minRight = mainLabelP;
	}

	//check expandable
	for (int left = maxLeft; left < startP; left++) {
		auto iter = motifEdge->begin();
		for (; iter != listEnd; ++iter) {
			if (!isSameLabel(*iter, left, startP)) {
				break;
			}
		}
		if (iter == listEnd) {
			for (int right = endP + 1; right <= minRight; right++) {
				auto iter = motifEdge->begin();
				for (; iter != listEnd; ++iter) {
					if (isSameLabel(*iter, right, startP)) {
						if (MORE(dif[right * nEdge + *iter] - dif[left * nEdge + *iter],
							Setting::delta*(right - left + 1))) {
							break;
						}
					}
					else break;
				}
				if (iter == listEnd) {
					cout << "motif startT=" << motif->getStartT() << " motif endT=" << motif->getEndT() << endl;
					cout << "expandable error: " << left << " " << right << ", edges: ";
					for (auto iter = motifEdge->begin(); iter != listEnd; ++iter)
						cout << (*iter) << endl;
					return;
				}
			}
		}
	}
}