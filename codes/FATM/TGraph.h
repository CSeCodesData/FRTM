#pragma once
#include "stdafx.h"
#include "TMotif.h"
#include "Util.h"
#include "Edge.h"
#include "LinkedList.h"
#include "DisjointSet.h"
#include "Strategy.h"
#include "CircularQueue.h"
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>

#define ALLOC_MEM 20000000
/*map the interval [Ts,Te] to the one dimension position of result*/
#define resultPos(Ts,Te,startT,endT,k) ((Ts - startT)*((endT<<1) - (startT<<1) - Ts - (k<<1) + 5) >>1) + Te - (Ts + k - 1)
#define getPos(d1,d2,d1size) d1*d1size+d2

using NVIntv = Pair;
int CircularQueue<NVIntv>::queueSize = 0;
int CircularQueue<int>::queueSize = 0;

using CCTYPEOPT1 = CComponentsFRTM;
#define ComponentsTypeOPT1 ComponentsType::CCFRTM

struct MidResult {
	CircularQueue<NVIntv> preMaxIntv;
	int row;
	NVIntv eMaxIntvl;
	int scanT;
	ForbidPairList vioT;
	//int scanNext;
	int edgeId;
	MidResult() = default;
	MidResult(int r, int e, NVIntv& intv/*, int next*/, CircularQueue<NVIntv>& preEMaxIntv, int scanTime, ForbidPairList& forbitT) {
		row = r;
		edgeId = e;
		eMaxIntvl = intv;
		//scanNext = next;
		preEMaxIntv.copyTo(preMaxIntv);
		scanT = scanTime;
		vioT.copyTo(forbitT);
	}
	~MidResult() {}
};


using MotifPos = pair<int, int>;//position of motif in result, id of motif in one position (interval)




/* Temporal Graph: the nodes and edges are fixed
   and the weights vary with time*/
class TGraph {
	
public:

	TGraph() {
		nNode = Setting::nodes;
		nEdge = Setting::edges;
		allNTimestamp = Setting::allNTimestamp;
		//ufset = nullptr;
	}
	// copy-constructed
	TGraph(const TGraph& ances);

	inline int getCurrNTimestamp() const { return currNTimestamp; }
	inline void setCurrNTimestamp(int cnt)  { currNTimestamp = cnt; }

	inline int getAllNTimestamp() const { return allNTimestamp; }

	inline int getNNode() const { return nNode; }

	inline int getNEdge() const { return nEdge; }

	inline int getStartT() const { return startT; }

	inline int getEndT() const { return endT; }

	inline void setStartT(int startT) { this->startT=startT; }

	inline void setEndT(int endT) { this->endT=endT; }

	inline void getEdgeList(NodePair*& s) { s = edgeList; }

	inline int getNumOfLabel() const { return numOfLabel; }

	inline int getLabelId(int label) { return labelToId[label]; }

	inline int getMaxIntervalLength() { return maxIntervalLength; }

	//inline void getDisjointSet(DisjointSet*& ufset) { ufset = ufsetKConnect; }

	/*testing for definition 1*/
	virtual void testLabel(TMotifI*& motif) {}

	//output the information of TGraph conveniently
	friend ostream& operator<<(ostream& output, const TGraph& e) {
		output << "Input temporal graph information:" << endl;
		output << "number of node: " << e.nNode << "\tnumber of edge: " << e.nEdge << endl;
		output << "time length: " << e.currNTimestamp <<
			"\tstart time: " << e.startT << "\tend time: " << e.endT << endl;
		return output;
	}

	/*load the setting of fixed labels */
	void readFixedLabel(bool*& fixLabel,
		const char* file, bool& isEdgeTypeFixed) {
		FILE* f;
		f = fopen(file, "r+");
		if (!f) {
			isEdgeTypeFixed = false;
			return;
		}
		else {
			char line[LINE_LENGTH];
			CLEARALL(line, 0, LINE_LENGTH, char);
			while (fgets(line, LINE_LENGTH, f)) {
				if (strlen(line) == 0) continue;
				int label = STR2INT(line);
				int labelId = labelToId[label];
				if (labelId < numOfLabel) {
					cout << label << " ";
					fixLabel[labelId] = true;
				}
				else {
					cout << "fixed labels are different from temporal graph" << endl;
				}
			}
			cout << endl;
			fclose(f);
		}
	}

	#pragma region construct and update the temporal graph
	void createCommonStruct();
	void releaseCommonStruct();

	void resetStruct();

	void createStructForDefType();
	void releaseStructForDefType();

	virtual ~TGraph() {
		releaseCommonStruct();
		releaseStructForDefType();
		delete[] edgeList;
		delete edge2ind;
		delete[] idToLabel;
		delete[] ind2node;
	}
	
	#pragma endregion

	/*print motif with specific edge*/
	void printMotif(TMotifII*& motif, int motifId, int node1 = -1, int node2 = -1);
	void printMotif2(TMotifII*& motif, int motifId);
	
	/*test whether a motif have all edges with same label in endpoints of interval*/
	virtual void checkMotifEndpointsD5(TMotifII*& motif) = 0;

		
	#pragma region edgeFilter
		#pragma region exact label matches
		virtual void edgeFilter(int intvB, int intvE,
			SAVEINFO_Vec*& selectedEdge, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		#pragma endregion

		#pragma region FRTM
		virtual void edgeFilterFRTM(int intvB, int intvE,
			vec(int)*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed, bool dynamicMode, int* edgeSetsRAdd = nullptr) = 0;
		virtual void edgeFilterFRTMMidRForDYN(int intvB, int intvE, int oriEndTE,
			vec(int)/*SAVESIMINFO_Vec*/*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed, int* edgeSetsRAdd = nullptr) = 0;
		#pragma endregion

		#pragma region FRTMOpt1
		virtual void edgeFilterShortIntvMidRForDYN(int intvB, int intvE, int limited, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, int oriEndTe, int lastTimeNoNoise,
			vec(int)*& selectedEdge, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed) = 0;//the same as edgeFilter but limit at [intvB,limitedE]
		virtual void edgeFilterShortIntv(int intvB, int intvE, int limited, vec(int)*& selectedEdge, int& selectedNum,
			bool*& fixLabel, bool isEdgeTypeFixed) = 0;//the same as edgeFilter but limit at [intvB,limitedE]
		virtual void edgeFilterShortIntvMidR(int intvB, int intvE, int limited, int lastTimeNoNoise,
			vec(int)*& selectedEdge, int& selectedNum, int choiceEndT, bool*& fixLabel, bool isEdgeTypeFixed) = 0;//the same as edgeFilter but limit at [intvB,limitedE]
		
		virtual void edgeFilterPlus(int intvB, int intvE, int filterE, int limited, bool*& hasE, int*& newE, int& newENum, int*& edgeSetsRAdd,
			vec(int)*& selectedEdge, int& selectedNum, vec(int)*& selectedEdgeNoEdge, int& selectedNumNoEdge,
			int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed, bool dynamicMode) = 0;
		virtual void edgeFilterPlusMidRForDYN(int intvB, int intvE, int filterE, int limited, int oriEndTE, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, int*& edgeSetsRAdd,
			vec(int)*& selectedEdge, int& selectedNum, vec(int)*& selectedEdgeNoEdge, int& selectedNumNoEdge, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		#pragma endregion
	#pragma endregion

	/*get the label of edge with edgeId at the timePos (time-startT)*/
	virtual int getEdgeLabel(int edgeId, int timePos) = 0;
	//has same label at timePos1 and timePos2
	inline virtual bool isSameLabel(int edgeId, int timePos1, int timePos2) = 0;
	virtual bool bothSameLabel(vec(int)& edges, int timePos1, int timePos2) = 0;
	virtual bool bothSameLabel(vec(int)& subCCs, vec(int)& edges, int timePos1, int timePos2) = 0;

	#pragma region maximal and non-expandable property checking
		#pragma region FRTMExact
		/*generate motifs in one interval for FRTMExact*/
		void genMotifInOneIntv(SAVEINFO_VIter& iterStart, SAVEINFO_VIter& iterEnd,
			i2iHMap& vertex2Pos, DisjointSet*& connectivity, int& vertexNum,
			vec(int)& combineCCPos,/* int& realMotifNum,*/
			i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
			vec(int)& saveCCPos, int motifStartT, int motifEndT,
			vec(TMotifI*)*& result, int pos, long long& motifNumber);

		/*fetch new edges from one R edge set and insert into the union-find set which maintains connectivity*/
		void maintainConnectivity(SAVEINFO_VIter& infoBegin,
			SAVEINFO_VIter& infoEnd, i2iHMap& vertex2Pos, DisjointSet*& connectivity, int&vertexNum, vec(int)& combineCCPos);
			
		/*combine two ccs*/
		void combineCC(CComponents*& origin, CComponents*& now);
		
		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponents(vec(CComponents*)& tempComponents,
			i2iHMap& vertex2Pos, /*DisjointSet*& disjointSet,*/DisjointSet*& connectivity,
			i2iHMap& root2Comp, /*int& tempComponentsSize,int& realMotifNum,*/
			vec(int)& combineCCPos);

		/*add edges into ccs*/
		void updateNewEdgeInfo(
			SAVEINFO_VIter& infoBegin, SAVEINFO_VIter& infoEnd,
			vec(CComponents*)& tempComponents,
			i2iHMap& vertex2Pos, DisjointSet*& connectivity,
			i2iHMap& root2Comp, /*int& tempComponentsSize, int& realMotifNum,*/
			vec(int)& saveCCPos, int startTime, int endTime);
		#pragma endregion	
		
		#pragma region maintain connectivity in maxCheck for FRTM and FRTMOPT1
		/*fetch new edges from one R edge set and insert into the union-find set which maintains connectivity*/
		void maintainConnectivity(veciter(int)& infoBegin,
			veciter(int)& infoEnd, /*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivity,/* int&vertexNum,*/ vec(int)& combineCCPos);
		
		#pragma endregion

		#pragma region maxCheck for FRTM

		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponentsFRTM(vec(CComponentsII*)& tempComponents,
			DisjointSet*& connectivity,
			i2iHMap& root2Comp, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
			/*int& tempComponentsSize,int& realMotifNum,*/
			vec(int)& combineCCPos);
		void combineComponentsFRTMOPT1(vec(CComponentsII*)& tempComponents,
			DisjointSet*& connectivity,
			i2iHMap& root2Comp, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
			/*int& tempComponentsSize,int& realMotifNum,*/
			vec(int)& combineCCPos, int*& haveNewEdgeCC);
		

		/*recompute ccs*/
		//void maintainTempCCForUF(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
		//	i2iHMap& vertex2Pos, DisjointSet* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTM* generatedCC);// union find set
		void maintainTempCCForUFWithTime(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
			DisjointSet* connectivity, i2iHMap& root2Comp, /*i2iHMap* subCCMap, */int*& subCCId, i2iHMap& root2Id, CComponentsFRTM* generatedCC, int currentTime);// union find set


		/*add edges into ccs*/
		void updateCCFRTM(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp,
			LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, 
			veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsType type);

		/*add edges into ccs and recompute ccs*/
		void updateNewEdgeInfoFRTM(
			veciter(int)& infoBegin, veciter(int)& infoEnd,
			vec(CComponentsII*)& tempComponents,
			/*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivity, DisjointSet* tempDisjointSet,
			i2iHMap& root2Comp, /*i2iHMap* subCCMap,*/ int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue, /*int& realMotifNum,*/
			vec(int)& saveCCPos, int startTime, int endTime, int k, ComponentsType type);
		#pragma endregion

		#pragma region maxCheck for FRTMOpt1 
			
		/*recompute ccs*/
		//bool maintainTempCCForUFOpt1(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue, unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,/* int*& lastEdgeSetR,*/DisjointSet* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTMMax* generatedCC, int endTime, int removeEdge, i2iHMap& haveNewEdgeCC, int*& remainEdges, int currentTime);// union find set
		
		template <class T>
		bool maintainTempCCForUFOpt1WithTime(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
			int*& edgesAdd, DisjointSet* connectivity, i2iHMap& root2Comp,  int*& subCCId, i2iHMap& root2Id, T*& generatedCC,
			int startTime, int endTime, int removeEdge, int*& haveNewEdgeCC, int currentTime);// union find set
		
		/*add edges into ccs*/
		void updateCCOpt1(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp, 
			/*unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,*/int*& edgesAdd, int*& haveNewEdgeCC, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsType ccstype);
		
		/*add edges into ccs and recompute ccs*/
		void updateNewEdgeInfoOpt1(
			veciter(int)& infoBegin, veciter(int)& infoEnd,/* unordered_set<int>& edgesAdd,*/int*& edgesAdd,
			vec(CComponentsII*)& tempComponents,
			/*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivity, DisjointSet* tempDisjointSet,
			i2iHMap& root2Comp, int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue,
			vec(int)& saveCCPos, int startTime, int endTime, int k, int*& haveNewEdgeCC, int*& tempHaveNewEdgeCC, ComponentsType ccstype);
#pragma endregion
			
		#pragma region maxCheck for FRTMPLUS
		/*fetch new edges from one R edge set and insert into the union-find set which maintains connectivity*/
		void maintainConnectivityPlus(veciter(int)& infoBegin,
					veciter(int)& infoEnd, int*& addE, int& addENum, /*i2iHMap& vertex2Pos,*/ DisjointSet*& connectivityNoNoise, /*int&vertexNum,*/
					bool*&/*iSet&*/ hasE, /*iSet&*/bool*& saveR, DisjointSet*& connectivity, vec(int)& combineCCPos);

		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponentsPlus(vec(CComponentsShortIntv*)& tempComponents,
					DisjointSet*& connectivity,
					i2iHMap& root2Comp, /*int& tempComponentsSize,int& realMotifNum,*/
					vec(int)& combineCCPos);

		/*add edges into ccs and recompute ccs*/
		void updateNewEdgeInfoPlus(
					veciter(int)& infoBegin,
					int*& addE, int addENum,
					vec(CComponentsShortIntv*)& tempComponents,
					/*i2iHMap& vertex2PosNoNoise,*/ DisjointSet*& connectivity,
					i2iHMap& root2Comp, bool*&/*iSet&*/ hasE, vec(int)& saveCCPos, int startTime, int endTime, ComponentsType type);

		/*record roots of ccs for the expandable property checking*/
		void recordSavedRoot(vec(CComponentsII*)& tempComponents, unordered_map<int, LinkedNode<int>*>& hasChecked, int*& newE, int& newENum, /*iSet&*/bool*& saveR, DisjointSet*& connectivity/*, i2iHMap& vertex2Pos*/);
		
		#pragma endregion


		#pragma region expCheck
		
		
		//checking each edge of the cc when checking the expandble property
		virtual bool bothFitDefAndSameLabelChangeStartTimePos(vec(int)& edges, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) = 0;
		virtual bool bothFitDefAndSameLabelChangeEndTimePos(vec(int)& edges, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) = 0;

		//expandable checking for ccs in each interval
		using CheckExpandableFRTM = void (TGraph::*)(int savePos, CComponentsFRTM* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask, vec(TMotifII*)*& result, long long& motifNumber);
		void expCheck(vec(int)&saveCCPos,
			vec(CComponents*)& tempComponents,
			vec(TMotifI*)*& result, int k, int motifStartT, int motifEndT,
			long long& motifNumber);


		template<typename T, typename CCType>
		void expCheckFRTM(vec(int)& saveCCPos,
			vec(CComponentsII*)& tempComponents,
			vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
			long long& motifNumber, T checkExpandable, CCType*& emptyCC);
		template<typename CCType>
		void expCheckFRTMPlus(vec(int)&saveCCPos,
			vec(CComponentsShortIntv*)& tempComponentsNoNoise,
			vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
			long long& motifNumber, DisjointSet*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked,
			i2iHMap& root2Comp,/* i2iHMap& vertex2Pos,*/ vec(CComponentsII*)& tempComponents, CCType*& emptyCC);
		template<typename CCType>
		void expCheckFRTMPlusMidR(vec(int)&saveCCPos,
			vec(CComponentsShortIntv*)& tempComponentsNoNoise,
			vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
			long long& motifNumber, DisjointSet*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked,
			i2iHMap& root2Comp,/* i2iHMap& vertex2Pos,*/ vec(CComponentsII*)& tempComponents, CCType*& emptyCC);
		
		
		//expandable checking for each cc
		template<typename T>
		void checkExpandableFRTM(int savePos, T* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask, vec(TMotifII*)*& result, long long& motifNumber);

		template<typename T>
		void checkExpandableFRTMMidR(int savePos, T* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask,  /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/vec(TMotifII*)*& result, long long& motifNumber);
		
		//void checkExpandableOpt1(int savePos, CComponentsII* tempCC,
		//	int motifStartT, int motifEndT, bool*& expandMask,  /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber);
		//void checkExpandableOpt1MidR(int savePos, CComponentsII* tempCC,
		//	int motifStartT, int motifEndT, bool*& expandMask,  /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber);

		#pragma endregion

	#pragma endregion 
		
	#pragma region incremental algorithm
	/*row number<=T-k+1*/
	void findRTMotifsDynamic(int k, vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed, ComponentsType cctype);
	void findRTMotifsDynamicOpt1(int k, vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed, ComponentsType cctype);
	void findRTMotifsDynamicPlus(int k, vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed, ComponentsType cctype);
	
	/*update the immediate result*/
	void updateMidResult();
	void showMidResult(vec(TMotifII*)*& result);

	virtual void updateDS(const char* src, int fixedE, int newFixedE) = 0;
	#pragma endregion
	
	
	void getMotifEdges(vec(TMotifI*)*& res, TMotifI*& motif, vector<int>& edgesList, int k) {
		int motifStartT = motif->getStartT();
		TMotifI* tempMotif;
		queue<TMotifI*> queue;
		queue.push(motif);
		while (!queue.empty()) {
			tempMotif = queue.front();
			queue.pop();

			vec(TEdge)* motifEdge = tempMotif->getMotifEdge();
			veciter(TEdge) listEnd = motifEdge->end();
			for (auto iter = motifEdge->begin();
				iter != listEnd; ++iter) {
				edgesList.emplace_back(iter->id);
			}
			if (tempMotif->getEndT() != endT) {
				//reuse edges in other motifs
				vec(SaveCCInfo)* otherEdge = tempMotif->getOtherEdge();
				if (otherEdge != nullptr) {
					veciter(SaveCCInfo) otherListEnd = otherEdge->end();
					for (auto iter = otherEdge->begin();
						iter != otherListEnd; ++iter) {
						int resultP = resultPos(motifStartT, iter->saveEndTime, startT, endT, k);
						queue.push((TMotifI*)res[resultP][iter->savePos]);
					}
				}
			}
		}
	}

#pragma region construct the temporal graph
	/*
	construct the temporal graph
	graph file src: each line of the file describe an edge in a timestamp.
	each line has four number: u,v,t,w (separated by ',' with no other space)
	for weight w of edge (u,v) in time t. Ids of node u and v are not guaranteed to be
	continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
	first, and followed by edges in time 1, 2...
	*/
	void constructGraph(const char* src, int k, int fixedE = -1);

	virtual void loadInfomation(const char* src, int fixedE, int k) = 0;
#pragma endregion
	virtual inline void lazyUpdate(int t, int pos, int edgeId) = 0;

	template<typename T>
	void releaseCCs(vec(T*)& setCC);
	//static int* degree;

private:
	//compute the intersection of scope[e] for each edge in the cc
	void getIntersectionOfScope(vec(int)& edges, pair<int,int>& intersectIntv);
	
	inline void connected(NodePair* temp, int edgeId, DisjointSet*& connectivity, /*i2iHMap& vertex2Pos,*/ vec(int)& combineCCPos/*, int&vertexNum*/);

	inline bool checkLeftOrBothExpandable(vec(int)& edges, int motifStartT, int motifEndT,
		pair<int, int>& intersectIntv, bool*& expandMask);
	inline bool checkRightOrLeftOrBothExpandable(vec(int)& edges, int motifStartT, int motifEndT,
		pair<int, int>& intersectIntv, bool*& expandMask);
	inline bool checkUpdateInMIntR(vec(int)& edges, veciter(int)& edgeEnd, int motifStartT);

	template<typename T>
	inline void removeEdgesAndUpdateScopes(T*& generatedCC, int*& subCCId, vector<pair<int, int>>& tempRecordFromQueue, int endTime, int& removeEdge);

	void printMotifEdges(vec(int)*& edgesList, vec(int)* removed, int motifStartT);
	
protected:
	#pragma region auxiliary variable and structure 
		int nNode, nEdge; // number of nodes and edges
		int currNTimestamp;	// number of snapshots(now)
		int allNTimestamp;	// number of snapshots(all)
		int startT, endT; // beginning timestamp and ending timestamp
		
		NodePair* edgeList;// map edges' id to (v1,v2)
		#define printE(id) cout<<id<<":"<<ind2node[edgeList[id].first]<<","<<ind2node[edgeList[id].second]<<endl;
		#define printENoId(id) cout<<ind2node[edgeList[id].first]<<","<<ind2node[edgeList[id].second]<<endl;
		#define printEAft(id) cout<<id<<":"<<edgeList[id].first<<","<<edgeList[id].second<<endl;
		#define printEAftNoId(id) cout<<edgeList[id].first<<","<<edgeList[id].second<<endl;

		absl::flat_hash_set<Edge>* edge2ind;// map from an edge to its index
		absl::flat_hash_map<int, int> node2ind;
		int* ind2node;

		/*vec(int) addEdge;
		unordered_set<int> delEdge;
		bool* edgeBase;
		DisjointSet** disSet;
		vec(int) combine;
		int* tempE;
		vector<vec(int)>* tempCC;
		unordered_map<int, int>* tempRoot2Id;*/

		absl::flat_hash_map<int, int> labelToId; // map the label to its label id
		int* idToLabel; // map the label id to its corresponding label
		int numOfLabel;//the number of different labels
		
		/*used in exact label matches*/
		SaveCCInfo* motifSaveInfo; 
		int maxIntervalLength;//used for reducing memory of R edge set

		/*used in edgeFilter*/
		ForbidPairList* vioT = nullptr;
		int* scanT;
		NVIntv*maxIntv; // current maximum relaxed interval for each edge
		NVIntv*maxIntvShortIntv; // current maximum relaxed interval for each edge for short intervals
		CircularQueue<NVIntv> *preMaxIntv;// previous maximum relaxed interval for each edge

		//scope[e] for each edge
		NVIntv* scope;

		/*for incremental algorithms*/
		vec(int)* edgesInEIntR = nullptr;//edges in the EIntR(used for edgeFilter in incremental algorithms)
		bool *edgesFetched = nullptr;//whether edges have been fetched from EIntR
		vec(bool) validMidResult, newValidMidResult;//the item in EIntR may be updated
		vec(MidResult*)* EIntR = nullptr, *newEIntR = nullptr;
		int *posInEIntR = nullptr, *newPosInEIntR = nullptr;//the position of edge stored in the EIntR for each edge
		vec(MotifPos)* MIntR = nullptr, *newMIntR = nullptr;

		//pair<int,int>* edgeBef = nullptr;
public:

	//static int saveEdgesControl;//save the disjoint set if the number of edges in sub-connected components > saveEdgesControl

	static int posUsedForEdgeFilter;
	static int posUsedForEdgeFilterShortIntv;
#pragma endregion
};


template<typename T>
void TGraph::releaseCCs(vec(T*)& setCC) {
	auto listEnd = setCC.end();
	for (auto listIter = setCC.begin();
		listIter != listEnd; ++listIter) {
		if (*listIter != nullptr) {
			delete *listIter;
		}
	}
	setCC.clear();
}


template<typename T>
void TGraph::checkExpandableFRTM(int savePos, T* tempCC,
	int motifStartT, int motifEndT, bool*& expandMask, vec(TMotifII*)*& result, long long& motifNumber) {
	//newly generated connected components
	pair<int, int> intersectIntv;
	if (tempCC->subCCNum == 0) {//not edges removed  (must be non right expandable)
		if (motifStartT != startT) {
			intersectIntv.first = tempCC->maxEMaxIntvlStartTQueue.top().first;
			intersectIntv.second = tempCC->minEMaxIntvlEndT;
			if (intersectIntv.first == motifStartT || //non left expandable => non both expandable
				checkLeftOrBothExpandable(tempCC->edges, motifStartT, motifEndT, intersectIntv, expandMask)//left expandable or both expandable
				) {
#ifdef TestMode
				Test::gne12 += tempCC->edges.size();
				if (intersectIntv.first == motifStartT) Test::gnefield++;
				else {
					int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first + 1);
					Test::gnefield += tempField * tempCC->edges.size();
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);
				}
#endif
				tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
			}
		}
		else {
#ifdef TestMode
			Test::gne11 += tempCC->edges.size();
			Test::gnefield++;
#endif
			tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
		}
		//Test::counter +=END_NSTIMER;
	}
	else {
		//BEGIN_TIMER(b)
		int subCCNum = tempCC->subCCNum;
		//check each subCC
		for (int i = 0; i < subCCNum; ++i) {
			if (tempCC->subCCsaved[i]) {//exists e¡ÊR+[m,i] (must be non right expandable)
				//BEGIN_NSTIMER; 
				if (motifStartT != 0) {
					getIntersectionOfScope(tempCC->subCCs[i], intersectIntv);
					if (intersectIntv.first == motifStartT || //non left expandable => non both expandable
						checkLeftOrBothExpandable(tempCC->subCCs[i], motifStartT, motifEndT, intersectIntv, expandMask)//left expandable or both expandable
						) {
#ifdef TestMode
						Test::gne212 += tempCC->subCCs[i].size();
						if (intersectIntv.first == motifStartT) Test::gnefield++;
						else {
							int tempField = (intersectIntv.second - motifEndT + 1) * (motifStartT - intersectIntv.first + 1);
							Test::gnefield += tempField * tempCC->subCCs[i].size();
							Test::gnemaxfield = max(Test::gnemaxfield, tempField);
						}
#endif
						tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
					}
				}
				else {
#ifdef TestMode
					Test::gne211 += tempCC->subCCs[i].size();
					Test::gnefield++;
#endif
					tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
				}
				//Test::counter += END_NSTIMER;
			}
			else {
				//BEGIN_NSTIMER; 
				getIntersectionOfScope(tempCC->subCCs[i], intersectIntv);
				if (intersectIntv.first == motifStartT) {//only check right expandable
#ifdef TestMode
					Test::gne221 += tempCC->subCCs[i].size();
					int tempField;
					if (intersectIntv.second - motifEndT - 1 > 0) {
						tempField = (intersectIntv.second - motifEndT - 1) * (motifStartT - intersectIntv.first + 1);
					}
					else tempField = 0;
					Test::gnefield += tempField * tempCC->subCCs[i].size();
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);
#endif
					if (motifStartT == 0 && intersectIntv.second == endT) continue;//right expandable
					CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
					if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
						tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
					}
					//Test::counter += END_NSTIMER;
				}
				else {//left expandable, right expandable or both expandable
#ifdef TestMode
					Test::gne222 += tempCC->subCCs[i].size();
					int tempField = (intersectIntv.second - motifEndT) * (motifStartT - intersectIntv.first);
					Test::gnefield += tempField * tempCC->subCCs[i].size();
					Test::gnemaxfield = max(Test::gnemaxfield, tempField);
#endif
					if (checkRightOrLeftOrBothExpandable(tempCC->subCCs[i], motifStartT, motifEndT, intersectIntv, expandMask)) {
						tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
					}
				}
			}
		}
		//if (!Test::counter5) Test::counter8 += END_TIMER(b);
		//else Test::counter9 += END_TIMER(b);
	}
}


template<typename T>
void TGraph::checkExpandableFRTMMidR(int savePos, T* tempCC,
	int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/  vec(TMotifII*)*& result, long long& motifNumber) {
	//newly generated connected components
	pair<int, int> intersectIntv;
	if (tempCC->subCCNum == 0) {//not edges removed  (must be non right expandable)
		auto edgeEnd = tempCC->edges.end();
		if (motifStartT != startT) {
			intersectIntv.first = tempCC->maxEMaxIntvlStartTQueue.top().first;
			intersectIntv.second = tempCC->minEMaxIntvlEndT;

			if (intersectIntv.first == motifStartT ||  //non left expandable => non both expandable
				checkLeftOrBothExpandable(tempCC->edges, motifStartT, motifEndT, intersectIntv, expandMask) //check left/both expandable
				) {
				if (checkUpdateInMIntR(tempCC->edges, edgeEnd, motifStartT)) {
					newMIntR->emplace_back(savePos, (int)result[savePos].size());
				}
				tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
			}
		}
		else {
			if (checkUpdateInMIntR(tempCC->edges, edgeEnd, motifStartT)) {
				newMIntR->emplace_back(savePos, (int)result[savePos].size());
			}
			tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
		}
		//Test::counter +=END_NSTIMER;
	}
	else {
		int subCCNum = tempCC->subCCNum;
		//check each subCC
		for (int i = 0; i < subCCNum; ++i) {
			veciter(int) subCCEnd = tempCC->subCCs[i].end();
			if (tempCC->subCCsaved[i]) {//exists e¡ÊR+[m,i] (must be non right expandable)
				//BEGIN_NSTIMER; 
				if (motifStartT != 0) {
					getIntersectionOfScope(tempCC->subCCs[i], intersectIntv);
					if (intersectIntv.first == motifStartT || //non left expandable => non both expandable
						checkLeftOrBothExpandable(tempCC->subCCs[i], motifStartT, motifEndT, intersectIntv, expandMask) //check left/both expandable
						) {
						if (checkUpdateInMIntR(tempCC->subCCs[i], subCCEnd, motifStartT)) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}
						tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
					}
				}
				else {
					if (checkUpdateInMIntR(tempCC->subCCs[i], subCCEnd, motifStartT)) {
						newMIntR->emplace_back(savePos, (int)result[savePos].size());
					}
					tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
				}
				//Test::counter += END_NSTIMER;
			}
			else {
				//BEGIN_NSTIMER; 
				getIntersectionOfScope(tempCC->subCCs[i], intersectIntv);
				if (intersectIntv.first == motifStartT) {//only check right expandable
					if (motifStartT == 0 && intersectIntv.second == endT) continue;//right expandable
					CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
					if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->subCCs[i], motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
						if (checkUpdateInMIntR(tempCC->subCCs[i], subCCEnd, motifStartT)) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}
						tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
					}
					//Test::counter += END_NSTIMER;
				}
				else {//left expandable, right expandable or both expandable
					if (checkRightOrLeftOrBothExpandable(tempCC->subCCs[i], motifStartT, motifEndT, intersectIntv, expandMask)) {
						if (checkUpdateInMIntR(tempCC->subCCs[i], subCCEnd, motifStartT)) {
							newMIntR->emplace_back(savePos, (int)result[savePos].size());
						}
						tempCC->saveToResultSubCC(motifStartT, motifEndT, i, result, savePos, motifNumber);
					}
					//Test::counter2 += END_NSTIMER;
				}
			}
		}
	}
}

template<typename T>
void TGraph::removeEdgesAndUpdateScopes(T*& generatedCC, int*& subCCId, vector<pair<int, int>>& tempRecordFromQueue, int endTime, int& removeEdge) {
	removeEdge = 0;
	int queueSize = (int)generatedCC->noisePosQueue->size();
	if (queueSize != 0) {
		CLEARALL(subCCId, 0, generatedCC->edges.size(), int);
	}

	generatedCC->tabuTChangePos = -1;
	tempRecordFromQueue.clear();
	while (queueSize != 0) {//worst time O(Em log (delta Em))
		auto ccLastNoisePair = generatedCC->noisePosQueue->top();

		if (ccLastNoisePair.first < endTime)
			break;
		else {
			if (ccLastNoisePair.first == endTime) {//edge has noise in this time

				removeEdge++;
				subCCId[ccLastNoisePair.second.pos] = -1;//remove this edge temporarily
				tempRecordFromQueue.emplace_back(ccLastNoisePair.second.startTime, // noise end position
					ccLastNoisePair.second.pos);
				if (generatedCC->tabuTChangePos < ccLastNoisePair.second.startTime - 1) {
					generatedCC->tabuTChangePos = ccLastNoisePair.second.startTime - 1;
				}
				generatedCC->noisePosQueue->pop();
			}
			else {
				generatedCC->noisePosQueue->pop();
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
		queueSize = (int)generatedCC->noisePosQueue->size();
	}

	if (generatedCC->tabuTChangePos != -1) {
		if (generatedCC->noisePosQueue->size() != 0) {
			auto ccLastNoisePair = generatedCC->noisePosQueue->top();
			generatedCC->tabuTChangePos = max(ccLastNoisePair.first, generatedCC->tabuTChangePos);
		}

		int nextTime = endTime - 1;
		for (auto item : tempRecordFromQueue) {
			if (nextTime >= item.first) {
				generatedCC->noisePosQueue->emplace(
					nextTime, // noise next position
					NoisePos(
						item.first, // noise end position
						item.second // edge position in newCC->edges
					)
				);
			}
		}
	}
	while (generatedCC->preEMaxIntvlChangePos->size() != 0) {//worst time O(log (delta Em))
		auto preEMaxIntvlPair = generatedCC->preEMaxIntvlChangePos->top();

		if (preEMaxIntvlPair.endTime < endTime)
			break;
		else {
			int id = preEMaxIntvlPair.id;
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

			generatedCC->preEMaxIntvlChangePos->pop();
		}
	}
}


template <class T>
bool TGraph::maintainTempCCForUFOpt1WithTime(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>& tempRecordFromQueue,
	/*unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,*/ int*& edgesAdd, DisjointSet* connectivity,
	i2iHMap& root2Comp, /*i2iHMap* subCCMap,*/ int*& subCCId, i2iHMap& root2Id, T*& generatedCC,
	int startTime, int endTime, int removeEdge, int*& haveNewEdgeCC, /*int*& remainEdges,*/ int currentTime) {
	//BEGIN_NSTIMER(p);
	DisjointSetWithTime* ds = (DisjointSetWithTime*)connectivity;
	auto tempEdges = &generatedCC->edges;
	int* subCCIter;
	int size = (int)tempEdges->size();
	veciter(int) motifEdgesEnd = tempEdges->end();
	int edgePos = 0, ccid;
	int subCCNum = 0, realSubCCNum = 0;
	int remainEdgesSize = size - removeEdge;

	//if (remainEdgesSize == 1) {//one cc
	//	if (generatedCC->subCCs == nullptr) {
	//		generatedCC->subCCs = DBG_NEW vec(int)[1];
	//		generatedCC->subCCsaved = DBG_NEW bool[1];
	//	}
	//	else {
	//		generatedCC->subCCs[0].clear();
	//	}
	//	generatedCC->subCCNum = 1;
	//	//generatedCC->subCCHaveNewEdges[0] = false;
	//	subCCIter = &subCCId[0];
	//	for (auto iter = generatedCC->edges.begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
	//		if (*subCCIter != -1) {
	//			if (edgesAdd[*iter] == startTime) {
	//				generatedCC->subCCs[0].emplace_back(*iter);
	//				generatedCC->subCCsaved[0] = edgePos >= generatedCC->newInsert;//true: right non-expandable
	//				return true;
	//			}
	//			else {
	//				return false;
	//			}
	//		}
	//	}
	//	//Test::counter2 += END_NSTIMER(p);
	//	return false;
	//}
	//else if (remainEdgesSize == 2) {
	//	bool isNewEdges1, isNewEdges2;
	//	subCCIter = &subCCId[0];
	//	int edge1Pos = -1, edge2Pos = -1, edge1V, edge1U, edge2V, edge2U, edgeId1, edgeId2;
	//	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
	//		if (*subCCIter != -1) {
	//			if (edge1Pos == -1) {
	//				edgeId1 = *iter;
	//				edge1Pos = edgePos;
	//				tie(edge1V, edge1U) = edgeList[edgeId1];
	//				//isNewEdges1 = edgesAdd.find(edgeId1) != edgesAddEnd;
	//				isNewEdges1 = edgesAdd[edgeId1] == startTime;
	//			}
	//			else {
	//				edgeId2 = *iter;
	//				edge2Pos = edgePos;
	//				tie(edge2V, edge2U) = edgeList[edgeId2];
	//				//isNewEdges2 = edgesAdd.find(edgeId2) != edgesAddEnd;
	//				isNewEdges2 = edgesAdd[edgeId2] == startTime;
	//			}
	//		}
	//	}
	//	if (!isNewEdges1 && !isNewEdges2) return false;
	//	else {
	//		if (edge1V == edge2V || edge1V == edge2U || edge1U == edge2V || edge1U == edge2U) {
	//			generatedCC->subCCNum = 1;
	//			if (generatedCC->subCCs == nullptr) {
	//				generatedCC->subCCs = DBG_NEW vec(int)[1];
	//				generatedCC->subCCsaved = DBG_NEW bool[1];
	//			}
	//			else {
	//				generatedCC->subCCs[0].clear();
	//			}
	//			generatedCC->subCCs[0].emplace_back(edgeId1);
	//			generatedCC->subCCs[0].emplace_back(edgeId2);
	//			generatedCC->subCCsaved[0] = max(edge1Pos, edge2Pos) >= generatedCC->newInsert ? true : false;
	//			//Test::counter2 += END_NSTIMER(p);
	//		}
	//		else {
	//			realSubCCNum = isNewEdges1 && isNewEdges2 ? 2 : 1;
	//			if (generatedCC->subCCs == nullptr) {
	//				generatedCC->subCCs = DBG_NEW vec(int)[realSubCCNum];
	//				generatedCC->subCCsaved = DBG_NEW bool[realSubCCNum];
	//			}
	//			else {
	//				generatedCC->subCCNum = (int)_msize(generatedCC->subCCsaved) / sizeof(bool);
	//				if (generatedCC->subCCNum < realSubCCNum) {
	//					delete[] generatedCC->subCCs;
	//					delete[] generatedCC->subCCsaved;
	//					generatedCC->subCCs = DBG_NEW vec(int)[realSubCCNum];
	//					generatedCC->subCCsaved = DBG_NEW bool[realSubCCNum];
	//				}
	//			}
	//			generatedCC->subCCNum = 0;
	//			if (isNewEdges1) {
	//				generatedCC->subCCs[generatedCC->subCCNum].clear();
	//				generatedCC->subCCs[generatedCC->subCCNum].emplace_back(edgeId1);
	//				generatedCC->subCCsaved[generatedCC->subCCNum] = edge1Pos >= generatedCC->newInsert;
	//				generatedCC->subCCNum++;
	//			}
	//			if (isNewEdges2) {
	//				generatedCC->subCCs[generatedCC->subCCNum].clear();
	//				generatedCC->subCCs[generatedCC->subCCNum].emplace_back(edgeId2);
	//				generatedCC->subCCsaved[generatedCC->subCCNum] = edge2Pos >= generatedCC->newInsert;
	//				generatedCC->subCCNum++;
	//			}
	//			//Test::counter2 += END_NSTIMER(p);
	//		}
	//		return true;
	//	}
	//}

	//BEGIN_NSTIMER(a);
	root2Id.clear();
	//get the number of new connected components
	subCCIter = &subCCId[0];
	//int remainP = 0;
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter) {//O(Em)
		if (*subCCIter != -1) {

			//BEGIN_NSTIMER(a1);
			auto tempEdge = &edgeList[*iter];

			/*union-find operation means that sId and tId are connected*/
			ds->addE(tempEdge->first, tempEdge->second, currentTime);
		}
	}
	//Test::counter += END_NSTIMER(a);
	//BEGIN_NSTIMER(q);
	//get subCCNum
	subCCIter = &subCCId[0];
	for (auto iter = tempEdges->begin(); iter != motifEdgesEnd; ++iter, ++subCCIter) {//O(Em)
		if (*subCCIter != -1) {
			/*the root of the edge's vertex in the disjoint set*/
			int root = ds->findRoot(edgeList[*iter].first);
			auto ccIt = root2Id.find(root);
			if (ccIt == root2Id.end()) {//new subCC
				root2Id[root] = subCCNum;
				*subCCIter = subCCNum;
				haveNewEdgeCC[subCCNum++] = (edgesAdd[*iter] == startTime) ? realSubCCNum++ : -1;
			}
			else {
				*subCCIter = ccIt->second;
				if (haveNewEdgeCC[*subCCIter] == -1 && edgesAdd[*iter] == startTime) {
					haveNewEdgeCC[*subCCIter] = realSubCCNum++;
				}
			}
		}
	}
	/*if (subCCNum == 1) { 
		Test::counter3++;
		Test::counter4+= tempEdges->size();
	}*/
	if (realSubCCNum == 0) { 
		//if (subCCNum == 1) Test::counter3 += END_TIMER(a);
		return false;
	}

	//BEGIN_NSTIMER(r);
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
	for (auto iter = generatedCC->edges.begin(); iter != motifEdgesEnd; ++iter, ++subCCIter, ++edgePos) {//O(Em)
		ccid = *subCCIter;
		if (ccid != -1) {
			int realCCId = haveNewEdgeCC[ccid];
			if (realCCId != -1) {
				generatedCC->subCCs[realCCId].emplace_back(*iter);
				generatedCC->subCCsaved[realCCId] = edgePos >= generatedCC->newInsert;//true: right non-expandable
			}
		}
	}
	//Test::counter3 += END_NSTIMER(r);
	return true;
}



template<typename T, typename CCType>
void TGraph::expCheckFRTM(vec(int)& saveCCPos,
	vec(CComponentsII*)& setCC,
	vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, 
	long long& motifNumber, T checkExpandable, CCType*& emptyCC) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	CComponentsII* tempCC;
	for (auto saveMotifIter : saveCCPos) {//need check expandable
		tempCC = setCC[saveMotifIter];
		//cout << motifStartT << " " << motifEndT << " " << *saveMotifIter << endl;
		
		(this->*checkExpandable)(savePos, (CCType*)tempCC, motifStartT, motifEndT, expandMask, result, motifNumber);

		int edgeNum = (int)tempCC->edges.size();

		if (edgeNum > tempCC->newInsert) {
			tempCC->newInsert = edgeNum;
		}
	}
	saveCCPos.clear();
}

template<typename CCType>
void TGraph::expCheckFRTMPlus(vec(int)& saveCCPos,
	vec(CComponentsShortIntv*)& setCCShortIntv,
	vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask,
	long long& motifNumber, DisjointSet*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked, 
	i2iHMap& root2Comp, /*i2iHMap& vertex2Pos,*/ vec(CComponentsII*)& setCC, CCType*& emptyCC) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	veciter(int) saveMotifEnd = saveCCPos.end();
	pair<int, int> intersectIntv;
	auto mapEnd = root2Comp.end();
	auto checkedEnd = hasChecked.end();
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		auto tempCC = setCCShortIntv[*saveMotifIter];
		if (tempCC->startT == motifStartT) {//left unexpandable
			intersectIntv.first = tempCC->scopeL;
			intersectIntv.second = tempCC->scopeR;

			if (intersectIntv.second == -1) {//only in [motifStartT, motifStartT+limited)
#ifdef TestMode
				Test::gnenonoise += tempCC->edges.size();
#endif
				tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
			}
			else {
				bool unsavedBefore = false;
				if (!tempCC->haveNonExpand) {
					int firstENode = edgeList[*tempCC->edges.begin()].first;
					//int cc = root2Comp[connectivity->findRoot(vertex2Pos[firstENode])];
					int cc = root2Comp[connectivity->findRoot(firstENode)];
					auto checkedCC = (CCType*)setCC[cc];
					auto checkIter = hasChecked.find(cc);
					unsavedBefore = checkIter == checkedEnd || (checkedCC->tabuTChangePos != -1 && checkedCC->tabuTChangePos < motifEndT);
				}

				if (tempCC->haveNonExpand || !unsavedBefore) {//check expandable
					if (intersectIntv.first == motifStartT) {//only check right expandable
#ifdef TestMode
						Test::gnenonoise += tempCC->edges.size();
						int tempField;
						if (intersectIntv.second - motifEndT - 1 > 0) {
							tempField = intersectIntv.second - motifEndT - 1;
						}
						else tempField = 0;
						Test::gnefield += tempField * tempCC->edges.size();
						Test::gnemaxfield = max(Test::gnemaxfield, tempField);
#endif

						CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
						if (!bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, motifStartT - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
							tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
						}
					}
					else {//left expandable, right expandable or both expandable
#ifdef TestMode
						Test::gnenonoise += tempCC->edges.size();
						int tempField = (intersectIntv.second - motifEndT) * (motifStartT - intersectIntv.first);
						Test::gnefield += tempField * tempCC->edges.size();
						Test::gnemaxfield = max(Test::gnemaxfield, tempField);/**/
#endif

						if (checkRightOrLeftOrBothExpandable(tempCC->edges, motifStartT, motifEndT, intersectIntv, expandMask)) {
							tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
						}
					}
				}
			}
		}
	}
	saveCCPos.clear();
}

template<typename CCType>
void TGraph::expCheckFRTMPlusMidR(vec(int)& saveCCPos,
	vec(CComponentsShortIntv*)& setCCShortIntv,
	vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
	long long& motifNumber, DisjointSet*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked, 
	i2iHMap& root2Comp, /*i2iHMap& vertex2Pos,*/ vec(CComponentsII*)& setCC, CCType*& emptyCC) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	veciter(int) saveMotifEnd = saveCCPos.end();
	TMotifII* motif;
	pair<int, int> intersectIntv;
	bool needCheck;
	auto mapEnd = root2Comp.end();
	auto checkedEnd = hasChecked.end();
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		auto tempCC = setCCShortIntv[*saveMotifIter];
		auto edgeEnd = tempCC->edges.end();

		if (tempCC->startT == motifStartT) {//left unexpandable
			intersectIntv.first = tempCC->scopeL;
			intersectIntv.second = tempCC->scopeR;
			if (intersectIntv.second == -1) {//only in [motifStartT, motifStartT+limited)
				if (checkUpdateInMIntR(tempCC->edges, edgeEnd, motifStartT)) {
					newMIntR->emplace_back(savePos, (int)result[savePos].size());
				}
				tempCC->saveToResultCC(motifStartT, motifEndT, result, savePos, motifNumber);
			}
			else {
				bool unsavedBefore = false;
				if (!tempCC->haveNonExpand) {
					int firstENode = edgeList[*tempCC->edges.begin()].first;
					//int cc = root2Comp[connectivity->findRoot(vertex2Pos[firstENode])];
					int cc = root2Comp[connectivity->findRoot(firstENode)];
					auto checkedCC = (CCType*)setCC[cc];
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
								if (needCheck && newPosInEIntR[edgeId * numOfLabel + getEdgeLabel(edgeId, motifStartT)] == EMaxIntvlChange::INTVINIT) { needCheck = false; break; }
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
						//bool isSave = true;
						//if (motifStartT - intersectIntv.first + 1 >= intersectIntv.second - motifEndT) {
						//	int j = intersectIntv.second;
						//	for (; j > motifEndT; --j) {
						//		CLEARALL(expandMask, true, motifStartT - intersectIntv.first + 1, bool);
						//		if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - startT, j - startT, motifStartT - startT, expandMask)) {//O(tEs)
						//			isSave = false;
						//			break;
						//		}
						//	}
						//}
						//else {
						//	int j = intersectIntv.first;
						//	for (; j <= motifStartT; ++j) {
						//		CLEARALL(expandMask, true, intersectIntv.second - motifEndT, bool);
						//		if (bothFitDefAndSameLabelChangeEndTimePos(tempCC->edges, j - startT, motifEndT + 1 - startT, intersectIntv.second - startT, motifStartT - startT, expandMask)) {//O(tEs)
						//			isSave = false;
						//			break;
						//		}
						//	}
						//}
						//if (isSave) {
						//	CLEARALL(expandMask, true, motifStartT - intersectIntv.first, bool);
						//	if (bothFitDefAndSameLabelChangeStartTimePos(tempCC->edges, intersectIntv.first - startT, motifStartT - 1 - startT, motifEndT - startT, motifStartT - startT, expandMask)) {//O(tEs)
						//		isSave = false;
						//	}
						//	if (isSave) {}
						if (checkRightOrLeftOrBothExpandable(tempCC->edges, motifStartT, motifEndT, intersectIntv, expandMask)) {
							
							motif = DBG_NEW TMotifII(motifStartT, motifEndT);

							needCheck = true;
							for (auto edgeId : tempCC->edges) {//O(motif number)
								if (needCheck && newPosInEIntR[edgeId * numOfLabel + getEdgeLabel(edgeId, motifStartT)] == EMaxIntvlChange::INTVINIT) { needCheck = false; break; }
							}
							if (needCheck) {
								newMIntR->emplace_back(savePos, (int)result[savePos].size());
							}

							motif->copyEdges(tempCC->edges);
							result[savePos].emplace_back(motif);
							motifNumber++;
						}
								//}
					}
				}
			}
		}
	}
	saveCCPos.clear();
}
