#pragma once
#include "stdafx.h"
#include "TMotif.h"
#include "Util.h"
#include "Edge.h"
#include "LinkedList.h"
#include "DisjointSet.h"
#include "Strategy.h"
#include "CircularQueue.h"

#define ALLOC_MEM 20000000
/*map the interval [Ts,Te] to the one dimension position of result*/
#define resultPos(Ts,Te,startT,endT,k) ((Ts - startT)*((endT<<1) - (startT<<1) - Ts - (k<<1) + 5) >>1) + Te - (Ts + k - 1)
#define getPos(d1,d2,d1size) d1*d1size+d2

using NVIntv = pair<int, int>;
int CircularQueue<NVIntv>::queueSize = 0;
int CircularQueue<int>::queueSize = 0;

struct MidResult {
	CircularQueue<NVIntv> preMaxIntv;
	int row;
	pair<int, int> eMaxIntvl;
	int scanT;
	ForbidPairList vioT;
	//int scanNext;
	int edgeId;
	MidResult() = default;
	MidResult(int r, int e, pair<int, int>& intv/*, int next*/, CircularQueue<NVIntv>& preEMaxIntv, int scanTime, ForbidPairList& forbitT) {
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

	/*testing for definition 1*/
	virtual void testLabel(TMotifI*& motif) {}
	/*testing for definition 2*/
	//virtual void testOverlap(TMotifI*& motif) {}

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
			vec(int)/*SAVESIMINFO_Vec*/*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		virtual void edgeFilterFRTMMidR(int intvB, int intvE,
			vec(int)/*SAVESIMINFO_Vec*/*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		virtual void edgeFilterFRTMMidRForDYN(int intvB, int intvE, int oriEndTE,
			vec(int)/*SAVESIMINFO_Vec*/*& edgeSetsR, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		#pragma endregion

		#pragma region FRTMOPT1
		virtual void edgeFilterOpt1(int intvB, int intvE, vec(int)*& edgeSetsRAdd,
			vec(int)/*SAVESIMINFO_Vec*/*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		virtual void edgeFilterOpt1MidR(int intvB, int intvE,  vec(int)*& edgeSetsRAdd, 
			vec(int)/*SAVESIMINFO_Vec*/*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		virtual void edgeFilterOpt1MidRForDYN(int intvB, int intvE, int oriEndTE, vec(int)*& edgeSetsRAdd,
			vec(int)*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		#pragma endregion

		#pragma region FRTMOpt1
		virtual void edgeFilterShortIntvMidRForDYN(int intvB, int intvE, int limited, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, int oriEndTe, int lastTimeNoNoise,
			vec(int)*& selectedEdge, int& selectedNum, bool*& fixLabel, bool isEdgeTypeFixed) = 0;//the same as edgeFilter but limit at [intvB,limitedE]
		virtual void edgeFilterShortIntv(int intvB, int intvE, int limited, vec(int)*& selectedEdge, int& selectedNum,
			bool*& fixLabel, bool isEdgeTypeFixed) = 0;//the same as edgeFilter but limit at [intvB,limitedE]
		virtual void edgeFilterShortIntvMidR(int intvB, int intvE, int limited, int lastTimeNoNoise,
			vec(int)*& selectedEdge, int& selectedNum, int choiceEndT, bool*& fixLabel, bool isEdgeTypeFixed) = 0;//the same as edgeFilter but limit at [intvB,limitedE]
		
		virtual void edgeFilterPlus(int intvB, int intvE, int filterE, int limited, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, vec(int)*& edgeSetsRAdd,
			vec(int)*& selectedEdge, int& selectedNum, vec(int)*& selectedEdgeNoEdge, int& selectedNumNoEdge,
			int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		virtual void edgeFilterPlusMidR(int intvB, int intvE, int filterE, int limited, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, vec(int)*& edgeSetsRAdd,
			vec(int)*& edgeSetsR, int& selectedNum, vec(int)*& selectedEdgeNoEdge, int& selectedNumNoEdge, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed) = 0;
		virtual void edgeFilterPlusMidRForDYN(int intvB, int intvE, int filterE, int limited, int oriEndTE, bool*&/*iSet&*/ hasE, int*& newE, int& newENum, vec(int)*& edgeSetsRAdd,
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
			i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, int& vertexNum,
			vec(int)& combineCCPos,/* int& realMotifNum,*/
			i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
			vec(int)& saveCCPos, int motifStartT, int motifEndT,
			vec(TMotifI*)*& result, int pos, long long& motifNumber);

		/*fetch new edges from one R edge set and insert into the union-find set which maintains connectivity*/
		void maintainConnectivity(SAVEINFO_VIter& infoBegin,
			SAVEINFO_VIter& infoEnd, i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, int&vertexNum, vec(int)& combineCCPos);
			
		/*combine two ccs*/
		void combineCC(CComponents*& origin, CComponents*& now);
		
		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponents(vec(CComponents*)& tempComponents,
			i2iHMap& vertex2Pos, /*DisjointSet*& disjointSet,*/DynamicConnectivity*& connectivity,
			i2iHMap& root2Comp, /*int& tempComponentsSize,int& realMotifNum,*/
			vec(int)& combineCCPos);

		/*add edges into ccs*/
		void updateNewEdgeInfo(
			SAVEINFO_VIter& infoBegin, SAVEINFO_VIter& infoEnd,
			vec(CComponents*)& tempComponents,
			i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity,
			i2iHMap& root2Comp, /*int& tempComponentsSize, int& realMotifNum,*/
			vec(int)& saveCCPos, int startTime, int endTime);
		#pragma endregion	
		
		#pragma region maintain connectivity in maxCheck for FRTM and FRTMOPT1
		/*fetch new edges from one R edge set and insert into the union-find set which maintains connectivity*/
		void maintainConnectivity(veciter(int)& infoBegin,
			veciter(int)& infoEnd, i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, int&vertexNum, vec(int)& combineCCPos);
		#pragma endregion

		#pragma region maxCheck for FRTM

		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponentsFRTM(vec(CComponentsII*)& tempComponents,
			i2iHMap& vertex2Pos, /*DisjointSet*& disjointSet, */DynamicConnectivity*& connectivity,
			i2iHMap& root2Comp, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
			/*int& tempComponentsSize,int& realMotifNum,*/
			vec(int)& combineCCPos);
			
		/*recompute ccs*/
		using MaintainCC = void (TGraph::*)(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
			i2iHMap& vertex2Pos, DynamicConnectivity* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTM* generatedCC);
		void maintainTempCCForUF(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
			i2iHMap& vertex2Pos, DynamicConnectivity* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTM* generatedCC);// union find set
			
		/*add edges into ccs*/
		using UpdateComponent = void (TGraph::*)(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp,
			LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,/* int& realMotifNum,*/
			veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsD5Type type);
		void updateCCFRTM(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp,
			LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, 
			veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsD5Type type);

		/*add edges into ccs and recompute ccs*/
		void updateNewEdgeInfoFRTM(
			veciter(int)& infoBegin, veciter(int)& infoEnd,
			vec(CComponentsII*)& tempComponents,
			i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, DisjointSet* tempDisjointSet,
			i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, vector<pair<int, int>>&tempRecordFromQueue, /*int& realMotifNum,*/
			vec(int)& saveCCPos, int startTime, int endTime, int k, MaintainCC strategy, UpdateComponent updateStrategy, ComponentsD5Type type);
		#pragma endregion

		#pragma region maxCheck for FRTMOpt1 
		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponentsOpt1(vec(CComponentsII*)& tempComponents,
			i2iHMap& vertex2Pos, /*DisjointSet*& disjointSet, */DynamicConnectivity*& connectivity,
			i2iHMap& root2Comp, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,
			/*int& tempComponentsSize,int& realMotifNum,*/
			vec(int)& combineCCPos);
			
		/*recompute ccs*/
		using MaintainCCOPT1 = bool (TGraph::*)(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue, unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,/* int*& lastEdgeSetR,*/ i2iHMap& vertex2Pos, DynamicConnectivity* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTMOPT1* generatedCC, int endTime, int removeEdge, i2iHMap& haveNewEdgeCC, int*& remainEdges);
		bool maintainTempCCForUFOpt1(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue, unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,/* int*& lastEdgeSetR,*/ i2iHMap& vertex2Pos, DynamicConnectivity* connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, CComponentsFRTMOPT1* generatedCC, int endTime, int removeEdge, i2iHMap& haveNewEdgeCC, int*& remainEdges/*, int startTime*/);// union find set
			
		/*add edges into ccs*/
		void updateCCOpt1(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp, unordered_set<int>& edgesAdd, unordered_set<int>::iterator& edgesAddEnd,/* int*& lastEdgeSetR,*/
			LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, /*int& realMotifNum,*/
			veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k);
			
		/*add edges into ccs and recompute ccs*/
		void updateNewEdgeInfoOpt1(
			veciter(int)& infoBegin, veciter(int)& infoEnd, unordered_set<int>& edgesAdd,/* int*& lastEdgeSetR,*/
			vec(CComponentsII*)& tempComponents,
			i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, DisjointSet* tempDisjointSet,
			i2iHMap& root2Comp, i2iHMap* subCCMap, int*& subCCId, i2iHMap& root2Id, LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked, /*bool*& expandMask,
			int*& maskLen,int& tempComponentsSize,*/vector<pair<int, int>>&tempRecordFromQueue, /*int& realMotifNum,*/
			vec(int)& saveCCPos, int startTime, int endTime, int k, i2iHMap& haveNewEdgeCC, int*& remainEdges, MaintainCCOPT1 strategy);
		#pragma endregion
			
		#pragma region maxCheck for FRTMPLUS
		/*fetch new edges from one R edge set and insert into the union-find set which maintains connectivity*/
		void maintainConnectivityPlus(veciter(int)& infoBegin,
					veciter(int)& infoEnd, int*& addE, int& addENum, i2iHMap& vertex2Pos, DynamicConnectivity*& connectivityNoNoise, int&vertexNum,
					bool*&/*iSet&*/ hasE, /*iSet&*/bool*& saveR, DynamicConnectivity*& connectivity, vec(int)& combineCCPos);

		/*combine ccs which are connected after adding new edges from one R edge set*/
		void combineComponentsPlus(vec(CComponentsShortIntv*)& tempComponents,
					i2iHMap& vertex2Pos, /*DisjointSet*& disjointSet,*/DynamicConnectivity*& connectivity,
					i2iHMap& root2Comp, /*int& tempComponentsSize,int& realMotifNum,*/
					vec(int)& combineCCPos);

		/*add edges into ccs and recompute ccs*/
		void updateNewEdgeInfoPlus(
					veciter(int)& infoBegin,
					int*& addE, int addENum,
					vec(CComponentsShortIntv*)& tempComponents,
					i2iHMap& vertex2PosNoNoise, DynamicConnectivity*& connectivity,
					i2iHMap& root2Comp, bool*&/*iSet&*/ hasE, vec(int)& saveCCPos, int startTime, int endTime);

		/*record roots of ccs for the expandable property checking*/
		void recordSavedRoot(vec(CComponentsII*)& tempComponents, unordered_map<int, LinkedNode<int>*>& hasChecked, int*& newE, int& newENum, /*iSet&*/bool*& saveR, DynamicConnectivity*& connectivity, i2iHMap& vertex2Pos);

		#pragma endregion


		#pragma region expCheck
		
		
		//checking each edge of the cc when checking the expandble property
		virtual bool bothFitDefAndSameLabelChangeStartTimePos(vec(int)& edges, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) = 0;
		virtual bool bothFitDefAndSameLabelChangeStartTimePos(vec(int)& edges, int*& subCCId, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) = 0;
		virtual bool bothFitDefAndSameLabelChangeStartTimePos(vec(int)& subCCs, vec(int)& edges, int startTimePos1, int startTimePos2, int endTimePos, int mainLabelPos, bool*& expandMask) = 0;
		virtual bool bothFitDefAndSameLabelChangeEndTimePos(vec(int)& subCCs, vec(int)& edges, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) = 0;
		virtual bool bothFitDefAndSameLabelChangeEndTimePos(vec(int)& edges, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) = 0;
		virtual bool bothFitDefAndSameLabelChangeEndTimePos(vec(int)& edges, int*& subCCId, int startTimePos, int endTimePos1, int endTimePos2, int mainLabelPos, bool*& expandMask) = 0;

		//expandable checking for ccs in each interval
		using CheckExpandable = void (TGraph::*)(int savePos, CComponentsII* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos, */vec(TMotifII*)*& result, long long& motifNumber);
		void expCheck(vec(int)&saveCCPos,
			vec(CComponents*)& tempComponents,
			vec(TMotifI*)*& result, int k, int motifStartT, int motifEndT,
			long long& motifNumber);
		void expCheckFRTM(vec(int)& saveCCPos,
			vec(CComponentsII*)& tempComponents,
			vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
			long long& motifNumber, CheckExpandable checkExpandable);
		void expCheckFRTMPlus(vec(int)&saveCCPos,
			vec(CComponentsShortIntv*)& tempComponentsNoNoise,
			vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
			long long& motifNumber, DynamicConnectivity*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked, i2iHMap& root2Comp, i2iHMap& vertex2Pos,
			vec(CComponentsII*)& tempComponents);
		void expCheckFRTMPlusMidR(vec(int)&saveCCPos,
			vec(CComponentsShortIntv*)& tempComponentsNoNoise,
			vec(TMotifII*)*& result, int k, int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck,*/
			long long& motifNumber, DynamicConnectivity*& connectivity, unordered_map<int, LinkedNode<int>*>& hasChecked, i2iHMap& root2Comp, i2iHMap& vertex2Pos, vec(CComponentsII*)& tempComponents);
		
		//expandable checking for each cc
		void checkExpandableFRTM(int savePos, CComponentsII* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask, /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber);
		void checkExpandableFRTMMidR(int savePos, CComponentsII* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask,  /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/vec(TMotifII*)*& result, long long& motifNumber);
		void checkExpandableOpt1(int savePos, CComponentsII* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask,  /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber);
		void checkExpandableOpt1MidR(int savePos, CComponentsII* tempCC,
			int motifStartT, int motifEndT, bool*& expandMask,  /*i2tupHMap& expCheck, ibPairVec_Iter ccPos,*/ vec(TMotifII*)*& result, long long& motifNumber);

		#pragma endregion

	#pragma endregion 
		
	#pragma region incremental algorithm
	/*row number<=T-k+1*/
	void findRTMotifsDynamic(int k, vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed);
	void findRTMotifsDynamicOpt1(int k, vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed);
	void findRTMotifsDynamicPlus(int k, vec(TMotifII*)*& newResult, int oriEndT, long long& motifNumber, bool*& fixLabel, bool isEdgeTypeFixed);
	
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
private:
	//compute the intersection of scope[e] for each edge in the cc
	void getIntersectionOfScope(vec(int)& edges, pair<int,int>& intersectIntv);
	void getIntersectionOfScope(vec(int)& edges, int*& subCCId, pair<int,int>& intersectIntv);
	void getIntersectionOfScope(vec(int)& subCCs, vec(int)& edges, pair<int,int>& intersectIntv);
	
protected:
	#pragma region auxiliary variable and structure 
		int nNode, nEdge; // number of nodes and edges
		int currNTimestamp;	// number of snapshots(now)
		int allNTimestamp;	// number of snapshots(all)
		int startT, endT; // beginning timestamp and ending timestamp
		
		NodePair* edgeList;// map edges' id to (v1,v2)
		#define printE(id) cout<<edgeList[id].first<<" "<<edgeList[id].second<<endl;

		set<Edge>* edge2ind;// map from an edge to its index

		unordered_map<int, int> labelToId; // map the label to its label id
		int* idToLabel; // map the label id to its corresponding label
		int numOfLabel;//the number of different labels
		
		/*used in exact label matches*/
		SaveCCInfo* motifSaveInfo; 
		int maxIntervalLength;//used for reducing memory of R edge set

		/*used in edgeFilter*/
		ForbidPairList* vioT = nullptr;
		int* scanT;
		pair<int, int> *maxIntv; // current maximum relaxed interval for each edge
		pair<int, int> *maxIntvShortIntv; // current maximum relaxed interval for each edge for short intervals
		CircularQueue<NVIntv> *preMaxIntv;// previous maximum relaxed interval for each edge

		//scope[e] for each edge
		pair<int,int>* scope;

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

