#pragma once
#include "stdafx.h"
#include "TGraph.h"


/*use Index table*/
class TGraphUDEL: public TGraph {
public:
	
	TGraphUDEL() :bef(nullptr), lab(nullptr), aft(nullptr), dif(nullptr), tail(nullptr) {
	}
	// copy-constructed
	TGraphUDEL(const TGraphUDEL& ances);
	
	//output the information of TGraph conveniently
	friend ostream& operator<<(ostream& output, const TGraphUDEL& e) {
		output << "Input temporal graph information:" << endl;
		output << "number of node: " << e.nNode << "\tnumber of edge: " << e.nEdge << endl;
		output << "time length: " << e.currNTimestamp <<
			"\tstart time: " << e.startT << "\tend time: " << e.endT << endl;
		return output;
	}

	virtual ~TGraphUDEL(){
		delete[]bef;
		delete[]lab; 
		delete[]aft; 
		delete[] dif;
		delete[] tail;
	}

	/*get the label of edge with edgeId at the timePos (time-startT)*/
	inline int getEdgeLabel(int edgeId, int timePos) {
		return lab[timePos*nEdge + edgeId];
	}

	//has same label at timePos1 and timePos2
	bool isSameLabel(int edgeId, int timePos1, int timePos2) {
		return lab[timePos1*nEdge + edgeId] == lab[timePos2*nEdge + edgeId];
	}
	bool bothSameLabel(vec(int)& edges, int timePos1, int timePos2) {
		auto edgeEnd = edges.end();
		for (auto edgeIter = edges.begin(); edgeIter != edgeEnd; ++edgeIter) {
			if (lab[timePos1*nEdge + *edgeIter] != lab[timePos2*nEdge + *edgeIter]) {
				return false;
			}
		}
		return true;
	}
	bool bothSameLabel(vec(int)& subCCs, vec(int)& edges, int timePos1, int timePos2) {
		auto edgeEnd = edges.end();
		auto subCCEnd = subCCs.end();
		int edgeId;
		for (auto subCCIter = subCCs.begin(); subCCIter != subCCEnd; ++subCCIter) {
			edgeId = edges[*subCCIter];
			if (lab[timePos1*nEdge + edgeId] != lab[timePos2*nEdge + edgeId]) {
				return false;
			}
		}
		return true;
	}


	/*test whether a motif have all edges with same label in endpoints of interval*/
	void checkMotifEndpointsD5(TMotifII*& motif);
	

	#pragma region construct and update the temporal graph
	void loadInfomation(const char* src, int fixedE, int k);
	//update the temporal graph
	void updateDS(const char* src, int fixedE, int newFixedE);
	#pragma endregion

	inline void lazyUpdate(int t, int pos, int edgeId) {
		int p = t - max(bef[pos], 1) + 1;
		if (t == p) return;
		int checkPos = p * nEdge + edgeId;
		int realAft = max(aft[pos], 1);
		if (realAft == aft[checkPos] - t + p) return;
		aft[pos] = aft[checkPos] - t + p;
	}
	/*edgeFilter for FRTMExact*/
	void edgeFilter(int intvB, int intvE,
		SAVEINFO_Vec*& edgeSetsR, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed);

	void edgeFilterFRTM(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed, bool dynamicMode, int* edgeSetsRAdd = nullptr);
	//for m <= T - k + 1
	void edgeFilterFRTMMidRForDYN(int intvB, int intvE, int oriEndTE,
		vec(int)*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed, int* edgeSetsRAdd = nullptr);


	//for short intervals
	void edgeFilterShortIntv(int intvB, int intvE, int limited,
		vec(int)*& edgeSetsR, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed);//the same as edgeFilter but limit at [intvB,limitedE]
	//for short intervals m > T - k + 1
	void edgeFilterShortIntvMidR(int intvB, int intvE, int limited, int lastTimeNoNoise,
		vec(int)*& edgeSetsR, int& selectedNum, int choiceEndT, bool*& fixLabel, bool isEdgeTypeFixed);//the same as edgeFilter but limit at [intvB,limitedE]
	//for short intervals m <= T - k + 1
	void edgeFilterShortIntvMidRForDYN(int intvB, int intvE, int limited, bool*& hasE, int*& newE, int& newENum, int oriEndT, int lastTimeNoNoise,
		vec(int)*& edgeSetsR, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed);//the same as edgeFilter but limit at [intvB,limitedE]


	void edgeFilterPlus(int intvB, int intvE, int filterE, int limited, bool*& hasE, int*& newE, int& newENum, int*& edgeSetsRAdd,
		vec(int)*& selectedEdge, int& selectedNum, vec(int)*& selectedEdgeNoEdge, int& selectedNumNoEdge,
		int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed, bool dynamicMode);
	void edgeFilterPlusMidRForDYN(int intvB, int intvE, int filterE, int limited, int oriEndTE, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, int*& edgeSetsRAdd,
		vec(int)*& selectedEdge, int& selectedNum, vec(int)*& selectedEdgeNoEdge, int& selectedNumNoEdge, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed);

private:
	inline void removeNonOverlapPreEMaxIntv(CircularQueue<NVIntv>& preEMaxIntvlPtr, int checkT);
	inline void removeNonOverlapPreEMaxIntvWithMem(CircularQueue<NVIntv>& preEMaxIntvlPtr, int checkT, int edgeId, int timestampPosForEMaxIntvlEndT, int& preMaxEMaxIntvlEndT/*, int& preMaxEMaxIntvlStartT*/);
	inline void removeNonOverlapPreEMaxIntvWithMem(CircularQueue<NVIntv>& preEMaxIntvlPtr, int checkT, int edgeId, int mainLabel, int timestampPosForEMaxIntvlEndT, int& preMaxEMaxIntvlEndT/*, int& preMaxEMaxIntvlStartT*/, int& maxEMaxIntvlEndT, int& maxEMaxIntvlStartT);
	template<class T>
	inline void saveToEdgeSetsR(int savePos, int edgeId, vec(T)*& edgeSetsR, ForbidPairList* noiseRecord, int& selectedNum, int& rightEndpoint);
	inline void saveToMidResult(int lastMainLabelPos, int edgeId, int localNoise, int mainLabelPos, int intvB, ForbidPairList& now);

	
	inline void getCurrentNoise(int currentPos, int mainLabel, int eid, int& localNoise, int& noiseNum);
	inline void updateVioT(int intvB, ForbidPairList& now, ForbidPairNode*& forbidIntv, int mainLabel, int eid, bool dynamicMode,
		int& localNoise, int& noiseNum);
	inline void scan(int currentPos, int beginPos, int endPos, ForbidPairList& now,
		int mainLabel, int mainLabelPos, int eid, int forbidTimeStartT, bool dynamicMode, int& localNoise,
		int& noiseNum, int& lastMainLabelPos);
	inline void fetchPreEMaxIntvl(CircularQueue<NVIntv>& preEMaxIntvlPtr, int eid, int mainLabel,
		int timestampPosForEMaxIntvlEndT, int& EMaxIntvlStartT, int& EMaxIntvlEndT);

	//fit definitions in [timePos1+startT, timePos2+startT]
	bool bothFitDefAndSameLabelChangeStartTimePos(vec(int)& edges, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask);
	bool bothFitDefAndSameLabelChangeEndTimePos(vec(int)& edges, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask);
protected:


	#pragma region DEL Table  
		int *lab; // lab_t:temporal graph label   lab[t+edgeId*T]  O(ET)
		int *bef;// len_t:the times of edges keeping their label fixed  bef[t+edgeId*T]  O(ET)
		int *aft;// aft_t: the times of edges keeping their label fixed in the time after t   aft[t+edgeId*T] O(ET)
		int *tail;//the last interval with different labels  aft[t+edgeId*T]  O(EL) 
		int* dif;//dif[e][t1]-dif[e][t2] = noise label number of lab[t1][e] in [t1,t2] if lab[t1][e]=lab[t2][e]
#pragma endregion  
};

