#include "TGraphUDEL.h"
#include "stdafx.h"

TGraphUDEL::TGraphUDEL(const TGraphUDEL& ances):TGraph(ances) {
	int size = nEdge * allNTimestamp;
	lab = DBG_NEW int[size];
	bef = DBG_NEW int[size];
	aft = DBG_NEW int[size];
	for (int i = 0; i < size; i++) {
		lab[i] = ances.lab[i];
		bef[i] = ances.bef[i];
		aft[i] = ances.aft[i];
	}
}

void TGraphUDEL::loadInfomation(const char* src, int fixedE, int k) {
	int size = allNTimestamp * nEdge, timestampPosForEL, timestampPosForELBre;

	lab = DBG_NEW int[size];
	bef = DBG_NEW int[size];
	aft = DBG_NEW int[size];
	for (int i = 0; i < size; i++) {
		lab[i] = 0x7fffffff;
	}
	
	int nodeNum = 0;//node num
	int flagT = -1;//check the same time as the first input data
	int u, v, t;
	int w;

	FILE* file;
	file = fopen(src, "r+");
	if (!file) {
		FILE_N_E
			exit(0);
	}
	char line[LINE_LENGTH];
	CLEARALL(line, 0, LINE_LENGTH, char);
	int sep1, sep2, sep3;//separator pos
	long ind = 0;
	//int vertexN = 0;
	numOfLabel = 0;
	this->edgeList = DBG_NEW NodePair[this->nEdge];
	
	edge2ind = DBG_NEW set<Edge>();
	vec(int) u_arr, v_arr, t_arr, w_arr;//allocate the memory
	u_arr.reserve(ALLOC_MEM);
	v_arr.reserve(ALLOC_MEM);
	t_arr.reserve(ALLOC_MEM);
	w_arr.reserve(ALLOC_MEM);
	int filterT = fixedE > 0 ? fixedE : allNTimestamp - 1;
	try {
		while (fgets(line, LINE_LENGTH, file)) {
			if (strlen(line) == 0) continue;
			sep1 = (int)(find(line, line + LINE_LENGTH, SEP_CHAR) - line);
			sep2 = (int)(find(line + sep1 + 1, line + LINE_LENGTH, SEP_CHAR) - line);
			sep3 = (int)(find(line + sep2 + 1, line + LINE_LENGTH, SEP_CHAR) - line);
			u = STR2INT(line);
			v = STR2INT(line + sep1 + 1);
			t = STR2INT(line + sep2 + 1);
			w = STR2INT(line + sep3 + 1);
			if (labelToId.find(w) == labelToId.end()) {
				labelToId[w] = numOfLabel++;
			}
			if (flagT == -1) flagT = t;

			Edge e(u, v);
			if (flagT == t) {//save static graph structure
				//edgeSet->insert(e1);
				set<Edge>::iterator e2IdIter = edge2ind->find(e);
				if (e2IdIter == edge2ind->end()) {
					e.id = ind;
					edgeList[ind] = NodePair(e.s, e.t);
					++ind;
					edge2ind->emplace(e);//save the edge of graph
				}
				else continue;//multiple edges

			}
			
			if (t <= filterT) {
				w = labelToId[w];
				u_arr.emplace_back(u);
				v_arr.emplace_back(v);
				t_arr.emplace_back(t);
				w_arr.emplace_back(w);
			}
			
			//cout << t << endl;
		}
		fclose(file);
	}
	catch (exception&) {
		LOAD_ERROR
	}

	BEGIN_TIMER(b)
	size_t length = u_arr.size();
	for (size_t i = 0; i < length; i++) {
		u = u_arr[i];
		v = v_arr[i];
		t = t_arr[i];
		w = w_arr[i];
		Edge e(u, v);
		int edgeInd = edge2ind->find(e)->id;
		timestampPosForEL = t * nEdge + edgeInd;
		if (lab[timestampPosForEL] == 0x7fffffff) {
			timestampPosForELBre = timestampPosForEL - nEdge;
			lab[timestampPosForEL] = w;
			if (t == 0) bef[edgeInd] = 1;
			else if (w == lab[timestampPosForELBre]) {
				bef[timestampPosForEL] = bef[timestampPosForELBre] + 1;
			}
			else bef[timestampPosForEL] = 1;
		}
	}

	this->startT = 0;
	this->endT = fixedE > 0 ? fixedE : (allNTimestamp - 1);
	this->currNTimestamp = fixedE > 0 ? (fixedE + 1) : allNTimestamp;
	int intvE, intvS;
	int eLabelSize = nEdge * numOfLabel;
	tail = DBG_NEW int[eLabelSize];
	CLEARALL(tail, -1, eLabelSize, int);
	int posForELabel = 0;
	dif = DBG_NEW int[size];
	CLEARALL(dif, 0, size, int);
	int* tempLabelsNum = DBG_NEW int[numOfLabel];
	int timestampPos = 0;
	for (int id = 0; id < nEdge; id++, posForELabel+= numOfLabel) {//O(E)
		intvE = currNTimestamp - 1;
		while (intvE >= 0) {
			intvS = intvE - bef[intvE*nEdge + id] + 1;
			for (int intvP = intvS; intvP <= intvE; intvP++) {
				aft[intvP*nEdge + id] = intvE - intvP + 1;
			}
			intvE = intvS - 1;
		}
		intvS = 0;
		while (intvS < currNTimestamp) {
			int lab = getEdgeLabel(id, intvS);
			intvE = intvS + aft[intvS*nEdge + id] - 1;
			if (tail[posForELabel + lab] == -1) {//first timestamp for label
				bef[intvS*nEdge + id] = -MYINFINITE;
			}
			else {
				bef[intvS*nEdge + id] = aft[tail[posForELabel + lab] * nEdge + id] = tail[posForELabel + lab] - intvS;
			}
			tail[posForELabel + lab] = intvE;
			aft[intvE * nEdge + id] = -MYINFINITE;
			intvS = intvE + 1;
		}

		CLEARALL(tempLabelsNum, 0, numOfLabel, int);
		timestampPos = 0;
		for (int j = 0; j < currNTimestamp; j++, timestampPos += nEdge) {
			int labelId = lab[timestampPos + id];
			tempLabelsNum[labelId]++;
			dif[timestampPos + id] = j + 1 - tempLabelsNum[labelId];
		}
	}
	
	auto iterEnd = this->labelToId.end();
	this->idToLabel = DBG_NEW int[this->labelToId.size()];
	for (auto iter = this->labelToId.begin(); iter != iterEnd; ++iter) {
		this->idToLabel[iter->second] = iter->first;
	}
	createStructForDefType();
	createCommonStruct();
	
	delete[] tempLabelsNum;
	
#ifdef _USENSTIMER
	cout << "preprocess: " << END_TIMER(b) << "ns" << endl;
#else
	cout << "preprocess: " << END_TIMER(b) << "ms" << endl;
#endif // _USENSTIMER

	//test the number of maximal intervals for each edge, where the edge has the same labels
	long long intvNum = 0, maxIntvNum = 0, eachIntvNum;
	//int avg = 0;
	//double all = 0;
	for (int id = 0; id < nEdge; id++, posForELabel += numOfLabel) {//O(E)
		int t = 0;
		eachIntvNum = 0;
		//avg = 0;
		while (t < currNTimestamp) {
			//avg += max(aft[t*nEdge + id], 1);
			t = t + max(aft[t*nEdge + id],1);
			eachIntvNum++;
		}
		//all += avg * 1.0 / eachIntvNum;
		intvNum += eachIntvNum;
		maxIntvNum = max(maxIntvNum, eachIntvNum);
	}
	cout << "maximum interval number: " << maxIntvNum << ", average interval number: " << intvNum * 1.0 / nEdge << endl;
	//cout<< all / nEdge << endl;
}

void TGraphUDEL::updateDS(const char* src, int fixedE, int newFixedE) {

	int u, v, t;
	int w;

	FILE* file;
	file = fopen(src, "r+");
	if (!file) {
		FILE_N_E
			exit(0);
	}
	char line[LINE_LENGTH];
	CLEARALL(line, 0, LINE_LENGTH, char);
	int sep1, sep2, sep3;//separator pos
	long ind = 0;
	vec(int) u_arr, v_arr, t_arr, w_arr;//allocate the memory
	u_arr.reserve(ALLOC_MEM);
	v_arr.reserve(ALLOC_MEM);
	t_arr.reserve(ALLOC_MEM);
	w_arr.reserve(ALLOC_MEM);
	try {
		while (fgets(line, LINE_LENGTH, file)) {
			if (strlen(line) == 0) continue;
			sep1 = (int)(find(line, line + LINE_LENGTH, SEP_CHAR) - line);
			sep2 = (int)(find(line + sep1 + 1, line + LINE_LENGTH, SEP_CHAR) - line);
			sep3 = (int)(find(line + sep2 + 1, line + LINE_LENGTH, SEP_CHAR) - line);
			u = STR2INT(line);
			v = STR2INT(line + sep1 + 1);
			t = STR2INT(line + sep2 + 1);
			w = STR2INT(line + sep3 + 1);

			if (t > fixedE && t <= newFixedE) {
				w = labelToId[w];
				u_arr.emplace_back(u);
				v_arr.emplace_back(v);
				t_arr.emplace_back(t);
				w_arr.emplace_back(w);
			}
		}
		fclose(file);
	}
	catch (exception&) {
		LOAD_ERROR
	}

	BEGIN_TIMER(a)
	size_t length = u_arr.size();
	for (size_t i = 0; i < length; i++) {
		u = u_arr[i];
		v = v_arr[i];
		t = t_arr[i];
		w = w_arr[i];
		Edge e(u, v);
		int edgeInd = edge2ind->find(e)->id;
		int timestampPosForEL = t * nEdge + edgeInd;
		if (lab[timestampPosForEL] == 0x7fffffff) {
			int timestampPosForELBre = timestampPosForEL - nEdge;
			lab[timestampPosForEL] = w; //update lab_t for t in [T+1,T+delta T]
			if (w == lab[timestampPosForELBre]) {  //update positive value of bef_t for t in [T+1,T+delta T]
				bef[timestampPosForEL] = max(bef[timestampPosForELBre],1) + 1;
			}
			else bef[timestampPosForEL] = 1;
		}
	}
	
	int intvE, intvS;
	int posForELabel = 0;
	int* tempLabelsNum = DBG_NEW int[numOfLabel];
	int timestampPos = 0;
	for (int id = 0; id < nEdge; id++, posForELabel += numOfLabel) {//O(E)
		
		intvE = currNTimestamp - 1;
		while (intvE > fixedE) { //update positive value of aft_t for t in [T+1,T+delta T]
			intvS = intvE - bef[intvE*nEdge + id] + 1;
			/*for (int intvP = intvS; intvP <= intvE; intvP++) {
				aft[intvP*nEdge + id] = intvE - intvP + 1;
			}*/
			for (int intvP = fixedE + 1; intvP <= intvE; intvP++) {
				aft[intvP*nEdge + id] = intvE - intvP + 1;
				//if(id == 1)cout << intvP << " " << aft[intvP*nEdge + id] << endl;
			}
			if (intvS <= fixedE) { //at most |L| times
				aft[intvS*nEdge + id] = intvE - intvS + 1; 
				//if (id == 1)cout << intvS << " " << aft[intvS*nEdge + id] << endl;
			}
			intvE = intvS - 1;
		}
		
		CLEARALL(tempLabelsNum, 0, numOfLabel, int);
		for (int lab = 0; lab < numOfLabel; lab++) {
			int lastP = tail[posForELabel + lab];
			if (lastP >= 0)
				tempLabelsNum[lab] = lastP - startT + 1 - dif[lastP * nEdge + id];
			else
				tempLabelsNum[lab] = 0;
		}
		
		while (intvS < currNTimestamp) { //update tail_lab and negative values of bef_t and aft_t (for aft_t, t may <= T, at most |L| times)
			int lab = getEdgeLabel(id, intvS);
			intvE = intvS + aft[intvS*nEdge + id] - 1;
			if (intvS > fixedE) {
				if (tail[posForELabel + lab] == -1) {//first timestamp for label
					bef[intvS*nEdge + id] = -MYINFINITE;
				}
				else {
					//if (id == 1)cout << tail[posForELabel + lab] << " ! " << tail[posForELabel + lab] - intvS << endl;
					bef[intvS*nEdge + id] = aft[tail[posForELabel + lab] * nEdge + id] = tail[posForELabel + lab] - intvS;
				}
			}
			tail[posForELabel + lab] = intvE;
			aft[intvE * nEdge + id] = -MYINFINITE;
			//if (id == 1)cout << intvE << " @ " << -MYINFINITE << endl;
			intvS = intvE + 1;
		}

		timestampPos = (fixedE + 1) * nEdge;
		for (int j = fixedE + 1; j < currNTimestamp; j++, timestampPos += nEdge) { //update dif_t for t in [T+1,T+delta T]
			int labelId = lab[timestampPos + id];
			tempLabelsNum[labelId]++;
			dif[timestampPos + id] = j + 1 - tempLabelsNum[labelId];
		}
	}
	
	delete[] tempLabelsNum;
	//delete[] maxNoiseNum;
#ifdef _USENSTIMER
	cout << "updatetime: " << END_TIMER(a) << "ns" << endl;
#else
	cout << "updatetime: " << END_TIMER(a) << "ms" << endl;
#endif // _USENSTIMER
}

/*compRES for exact label matches*/
void TGraphUDEL::edgeFilter(int intvB, int intvE,
	SAVEINFO_Vec*& edgeSetsR, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed) {
	int timePos = intvE - startT;
	int check;
	int edgeType;
	//unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int intvLen = intvE - intvB + 1;
	int endPos = endT - startT;
	int intvStart;
	for (int i = 0; i < nEdge; i++) {
		edgeType = lab[TGraph::posUsedForEdgeFilter + i];
		/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
			put here for less time*/
		if (max(bef[TGraph::posUsedForEdgeFilter + i],1) < intvLen /*||
			(isEdgeTypeFixed && fixLabel.find(edgeType)
				== fixLabelEnd)*/) {
			continue;
		}
		intvStart = intvE - max(bef[TGraph::posUsedForEdgeFilter + i],1) + 1;
		//intvStart = startIdx[i][eLPos] - startT;
		int j = timePos + 1;
		//Test::edgesNum++;
		if (maxIntv[i].second >= timePos && maxIntv[i].first <= intvB - startT) {//case 4
			check = maxIntv[i].second - timePos;
			selectedNum++;
			edgeSetsR[check].emplace_back(i,
				intvStart, edgeType);
			//Test::allCost2 += 1;
			//Test::allCost1 += check;
		}
		else if (maxIntv[i].second >= intvB - startT && maxIntv[i].first <= timePos) {//case 3
			continue;
		}
		else {//case 1 and 2
			maxIntv[i].first = intvStart;
			selectedNum++;
			
			check = max(aft[TGraph::posUsedForEdgeFilter + i],1) - 1;
			maxIntv[i].second = timePos + check + startT;
			edgeSetsR[check].emplace_back(i,
				intvStart, edgeType);
			
			//Test::maxCost = max(Test::maxCost, check);
			//Test::allCost2 += check;
			//Test::allCost1 += check;
		}
	}
}
