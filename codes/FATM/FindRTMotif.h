#pragma once
#include "FindTMotif.h"

class FindRTMotif : public FindTMotif {

public:
	using ComputeRESForRTM = void (TGraph::*)(int intvB, int intvE,
		vec(int)*& selectedEdge, int& selectedNum, int& rightEndpoint, int k, bool*& fixLabel, bool isEdgeTypeFixed);
	using UpdateComponent = void (TGraph::*)(vec(CComponentsII*)& tempComponents, i2iHMap&  root2Comp,
		LinkedList<int>*& needCheckedCC, unordered_map<int, LinkedNode<int>*>& hasChecked,/* int& realMotifNum,*/
		veciter(int)& infoIter, int id, int root, int filterTime, int startTime, int endTime, int k, ComponentsD5Type type);
	using MaintainCC = void (TGraph::*)(LinkedNode<int>*& tempCCIter, vector<pair<int, int>>&tempRecordFromQueue,
		i2iHMap& vertex2Pos, DynamicConnectivity*& connectivity, i2iHMap& root2Comp, i2iHMap* subCCMap, i2iHMap& root2Id, CComponentsFRTM* generatedCC);

#pragma region static algorithm
	static void FRTM(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel);
	static void FRTMOpt1(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel);
	static void FRTMPlus(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel);
#pragma endregion 

#pragma region revised static algorithm
	static void FRTMDYN(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT);
	static void FRTMOpt1DYN(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT);
	static void FRTMPlusDYN(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int choiceStartT, int choiceEndT);
#pragma endregion 
#pragma region incremental algorithm
	static void DFRTM(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int graphEndT);
	static void DFRTMOpt1(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int graphEndT);
	static void DFRTMPLus(TGraph*& graph,
		vec(TMotifII*)*& result, bool*& fixLabel, int graphEndT);
#pragma endregion

};