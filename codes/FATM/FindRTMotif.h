#pragma once
#include "FindTMotif.h"

class FindRTMotif : public FindTMotif {

public:
	using ComputeRESForRTM = void (TGraph::*)(int intvB, int intvE,
		vec(int)*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed);
	using UpdateComponent = void (TGraph::*)(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp,
		LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,/* int& realMotifNum,*/
		veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsType type);
	using MaintainCC = void (TGraph::*)(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
		i2iHMap& vertex2Pos, DisjointSet*& connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, i2iHMap& root2Id, CComponentsFRTM* generatedCC);

	static void reorganizeResult(TGraph*& graph, vec(TMotifII*)*& resultTF, int graphEndT, int choiceEndT);

#pragma region static algorithm
	static void FRTM(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT, bool dynamicMode, ComponentsType cctype);
	static void FRTMOpt1(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT, bool dynamicMode, ComponentsType cctype);
	static void FRTMPlus(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT, bool dynamicMode, ComponentsType cctype);
#pragma endregion 

#pragma region incremental algorithm
	static void DFRTM(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int graphEndT, ComponentsType cctype);
	static void DFRTMOpt1(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int graphEndT, ComponentsType cctype);
	static void DFRTMPLus(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int graphEndT, ComponentsType cctype);
#pragma endregion

};

