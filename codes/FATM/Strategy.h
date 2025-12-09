#pragma once


//the type of algorithm
enum AlgorithmType : uint8_t{
	INIT = 0, //initialized
	TESTING = 6, //testing parameters
	
	FRTM = 12, //FRTM
	DFRTM = 19, //DFRTM
	FRTMFORDYN = 20, //FRTM revision for DFRTM
	
	FRTMOPT1 = 22, //FRTM with common edge identifying
	FRTMOPT1DYN = 23, //FRTMOPT1 revision for DFRTMOPT1
	DFRTMOPT1 = 24, //DFTM with common edge identifying
	
	FRTMPLUS = 31, //FRTM with two optimization strategies
	FRTMPLUSDYN = 33, //FRTMPLUS revision for DFRTMPLUS
	DFRTMPLUS = 34, //DFRTM with two optimization strategies
};

class Setting {
public:
	static AlgorithmType choice;
	static int nodes, edges, allNTimestamp;
	static double delta;//global relaxation bound
	static int c;//local relaxation bound
	//static bool compressMode;//save motif with less memory
};

enum EMaxIntvlChange : int {
	INTVINIT = -3,
	CHANGED = -1,
	UNCHANGED = -2
};