#include"stdafx.h"
#include"Util.h"
#include"FindRTMotif.h"
#include"TGraphUDEL.h"
#include"TMotif.h"

#pragma region parameters initialization
int FindTMotif::k = DEFAULT_K;
long long FindTMotif::motifNumber = 0;
int FindTMotif::output = 0;
bool FindTMotif::isEdgeTypeFixed = false;
char FindTMotif::outputSrc[FILE_NAME_LENGTH] = OUTPUT_FILE;
long long FindTMotif::motifMaxNum = 0;
long long FindTMotif::motifSum = 0;
//int FindTMotif::filterIntvSize = 0;
//int FindTMotif::filterEdgeNumber = 0;
//int FindTMotif::fixednode1 = -1;
//int FindTMotif::fixednode2 = -1;
int TGraph::posUsedForEdgeFilter = 0;
int TGraph::posUsedForEdgeFilterShortIntv = 0;
//int TGraph::saveEdgesControl = 0;
AlgorithmType Setting::choice = AlgorithmType::INIT;//choose which algorithm to run
long long Test::fproc, Test::compr, Test::gm, Test::gne;
long long Test::ckf, Test::gmni, Test::gmli;
clock_t Test::msTimer; 
//long long Test::comprp1t = 0, Test::comprp2t = 0, Test::comprp1n = 0, Test::comprp2n = 0;
long long Test::gne11 = 0, Test::gne121 = 0, Test::gne122 = 0;
long long Test::gne211 = 0, Test::gne2121 = 0, Test::gne2122 = 0, Test::gne221 = 0, Test::gne222 = 0, Test::gnenonoise = 0;
long long Test::gnefield = 0, Test::gnemaxfield = 0;
double Setting::delta = -1;
int Setting::c = 0x7fffffff;
int Setting::nodes = -1, Setting::edges = -1, Setting::allNTimestamp = -1;
int Test::testingMode = 0;//testing mode
long long Test::counter = 0;
long long Test::counter2 = 0;
long long Test::counter3 = 0;
long long Test::counter4 = 0;
long long Test::counter5 = 0;
long long Test::counter6 = 0;
long long Test::counter7 = 0;
long long Test::counter8 = 0;
double Test::counter9 = 0;
double Test::counter10 = 0;
long long Test::sumIntvLen, Test::maxIntvLen;
long long Test::allContainment = 0, Test::containment = 0, Test::selectedSum= 0;

SIZE_T Test::peakMemory = 0;
#pragma endregion 

#pragma region function declaration
#pragma region load temporal graph
void createTGraph(TGraph*& temporal_graph, const char* inputSrc, int indexId, int endT);
#pragma endregion

#pragma region which problem to run
void runStaticAlgorithm(TGraph*& temporal_graph, bool*& fixLabel);
void testD5(TGraph*& temporal_graph, bool*& fixLabel);
void runIncrementalAlgorithm(TGraph*& temporal_graph, const char * src
	, int graphEndT,bool*& fixLabel, int newEndT); 
//void testStaticAlgorithm(TGraph*& temporal_graph/*, unordered_map<int, bool>& fixLabel*/);
#pragma endregion
#pragma region test
void testEndpointD5(TGraph*& temporalgraph, vec(TMotifII*)& lis);
#pragma endregion
#pragma region output
//void outputResult(TGraph*& temporal_graph, vec(TMotifI*)*& result, int resultLen);
void outputResult(TGraph*& temporal_graph, vec(TMotifII*)*& result, int resultLen);
void compareResult(TGraph*& temporal_graph, vec(TMotifI*)*& result, vec(TMotifII*)*& result2, int resultLen, long long motifNum1, long long motifNum2);
void compareResult(TGraph*& temporal_graph, vec(TMotifII*)*& result, vec(TMotifII*)*& result2, int resultLen, long long motifNum1, long long motifNum2);
void countingResult(TGraph*& temporal_graph, vec(TMotifII*)*& result, int resultLen, long long motifNum2);
#pragma endregion

void readK(vec(int)& arr, const char* file);

#pragma region release memory
void releaseResult(vec(TMotifI*)*& result, int len);
void releaseResult(vec(TMotifII*)*& result, int len);
#pragma endregion
#pragma endregion 
/////////////////////////////////////////////////////////////////////////////
/*Input Parameters:
- v : The number of nodes
- e : The number of edges
- t : The number of snapshots
- i : File name of input temporal graph(default: example-temporal-graph.txt)
required: information of all edges are sorted by timestamp, and the interval of temporal graph is [0,t-1]
- f : File name of output(default: example-output.txt)
- k : File name of frequency threshold(default:example-k.txt,if not exists k=10)
		the program can test different frequent conditions for one time(untested)
	for static algorithm:
		file content format:
			10
			20
			30
- r : Algorithm Id (default:1)
		6:TESTING (testing parameters)

		12:FRTM
		19:DFRTM
		20:FRTMFORDYN (FRTM saving information for DFRTM)
	
		22:FRTMOPT1 (FRTM using first strategies)
		23:FRTMOPT1DYN (FRTM using first strategies for dynamic scenes)
		24:DFRTMOPT1 (DFRTM using first strategies)
	
		31:FRTMPLUS (FRTM using two strategies)
		33:FRTMPLUSDYN (FRTM using two strategies for dynamic scenes)
		34:DFRTMPLUS (DFRTM using two strategies)
- o : Output level of motif (default:0)
		three levels:
			0:only output motif number, the running time and memory use
			1:except for those outputs mentioned above, output the edges number of every motif
			2:except for those outputs mentioned above, output the detailed information of motif edges with their labels (output labels)
			//3:except for those outputs mentioned above, output the detailed information of motif edges with their labels (real labels)
- g : Index id (default:2)
		2:index table
- l : Limit the start time and end time of input data
		(default:doesn't limit) format:-l:500
		e.g.: if the interval of input temporal graph is [0,100]
		-l:50  means you only use the interval [0,50]
			of temporal graph when testing the algorithm  (used in dynamic algorithm)
- n : Limit the end time of input data when snapshot increasing
		(use - l at the same time) (default:doesn't limit)
		format:-l:500 -n:1600
		 means that the interval of temporal graph is [0,500] before snapshots increase,
		 and that the interval of temporal graph is [0,1600] after snapshots increase (used in dynamic algorithm)
- a : Global approximation parameters
	format: -a:value(double), control the total proposition of label mismatches
    e.g.: -a:0.04 means the total proposition of label mismatches in each edge of motifs is less than or equal to 4%

- b : Local relaxtion parameters
	format: -b:value(integer)control the number of continuous label mismatches
	e.g.: -b:3 means the number of continuous label mismatches in each edge of motifs is less than or equal to 3

//- z : Testing mode
	//0: not testing mode
	//1: testing mode (test definition)
    //2: counting mode (counting noise-tolerant motifs relation with original motifs)
	//3: counting mode2 (counting noise-tolerant motifs containment relations)
//- p : output part of output for testing, limit output motifs with interval length >= value 
//	default:0
//- m : output part of output for testing, limit output motifs with edge number >= value
//	default:0
//- q : output part of output for testing, limit output motifs with a specific edge
//	default:not limitation
//	format: -q:node1,node2
//- x : File name of fixed labels of motif edges (default: do not fix labels)
//		file content format:
//			1
//			2 
//			(means if the network has a label type set = {1,2,3}, motifs have all edges with label type in {1,2},
//			one line for one label type)
*/////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

#pragma region parameter settings
	char inputSrc[FILE_NAME_LENGTH] = INPUT_FILE;//input temporal graph file
	char kFile[FILE_NAME_LENGTH] = K_FILE;//file of frequent conditon k
	FindTMotif::k = DEFAULT_K;//frequent condition
	char fixedLabelFile[FILE_NAME_LENGTH]
		= FIXEDLABEL_FILE;//file of fixed labels
	int indexId = INDEXID;//index id
	//int runTimes = RUNTIMES;//the testing times
	int endT = -1;//end time of the temporal graph
	int limitEndT = -1;//limit end time when running special case 2
	int newEndT = -1;// new end time of the temporal graph(after snapshot increasing)
	vec(int) fixedK;// multiple frequent conditions
	FindTMotif::isEdgeTypeFixed = false;//whether the edge label is fixed
	//unordered_map<int, bool> fixedLabel;//record fixed labels for one edge
#pragma endregion
#pragma region load parameters
	for (int i = 1; i < argc; i++) {
		switch (argv[i][1]) {
		case 'a'://global approximation parameters
			Setting::delta = STR2DOU(argv[i] + 3);
			if (Setting::delta < 0) {
				WRONG_A_P
				return -1;
			}
			break;
		case 'b'://local approximation parameters
			Setting::c = STR2INT(argv[i] + 3);
			if (Setting::c < 0) {
				WRONG_A_P
				return -1;
			}
			break;
		case 'x'://Fixed the label of edge
			strcpy(fixedLabelFile, argv[i] + 3);
			FindTMotif::isEdgeTypeFixed = true;
			break;
		case 'v': //the number of nodes
			Setting::nodes = STR2INT(argv[i] + 3);
			break;
		case 'e': //the number of edges
			Setting::edges = STR2INT(argv[i] + 3);
			break;
		case 't': //the number of snapshots
			Setting::allNTimestamp = STR2INT(argv[i] + 3);
			break;
		case 'f'://output file name
			strcpy(FindTMotif::outputSrc, argv[i] + 3);
			break;
		case 'i'://File name of input graph
			strcpy(inputSrc, argv[i] + 3);
			break;
		case 'k'://Frequency threshold
			strcpy(kFile, argv[i] + 3);
			break;
		case 'l'://fixed endT before updating
			endT = STR2INT(argv[i] + 3);
			break;
		case 'n'://fixed newEndT after updating
			newEndT = STR2INT(argv[i] + 3);
			break;
		case 'o'://output the motif or not
			FindTMotif::output = STR2INT(argv[i] + 3);
			break;
		//case 's'://output the number of sub-connected components
		//	TGraph::saveEdgesControl = STR2INT(argv[i] + 3);
		//	break;
		//case 'p'://limit the minimum interval length of output motifs
		//	FindTMotif::filterIntvSize = STR2INT(argv[i] + 3);
		//	break;
		//case 'm'://limit the minimum edge number of output motifs
		//	FindTMotif::filterEdgeNumber = STR2INT(argv[i] + 3);
		//	break;
		//case 'q'://limlt output motifs with a specific edge
		//	sep = (int)(find(argv[i] + 3, argv[i] + 3 + LINE_LENGTH, SEP_CHAR) - (argv[i] + 3));
		//	FindTMotif::fixednode1 = STR2INT(argv[i] + 3);
		//	FindTMotif::fixednode2 = STR2INT(argv[i] + 4 + sep);
		//	break;
		case 'r'://Choose Algorithm Id
			Setting::choice = (AlgorithmType)STR2INT(argv[i] + 3);
			break; 
		case 'z'://testing mode
			Test::testingMode = STR2INT(argv[i] + 3);
			break;
		default://do nothing
			break;
		}
	}

	OUTPUT_FILE_OPEN
	//argument setting
	cout << "graph : " << inputSrc << endl;
	readK(fixedK, kFile);
	cout << "choice : " << Setting::choice << ",indexId : "<< indexId<< endl;
	cout  << "global relaxation bound : " << Setting::delta << endl;
	cout << "local relaxation bound : " << Setting::c << endl;
	if (endT != -1) {
		cout << "fixed end time : " <<  endT << endl;
	}
	if (Setting::nodes == -1 || Setting::edges == -1 || Setting::allNTimestamp == -1) {
		cout << "temporal graph information(the number of nodes, edges and snapshots) are needed" << endl;
		return -1;
	} 

#pragma endregion

#pragma region load graph
	TGraph* temporal_graph;
	createTGraph(temporal_graph, inputSrc, indexId, endT);
#pragma endregion

	int labelN = temporal_graph->getNumOfLabel();
	bool* fixedLabel = DBG_NEW bool[labelN];//record fixed labels for one edge
#pragma region load fixed label
	if (FindTMotif::isEdgeTypeFixed) {
		cout << "fixed label of edge : ";
		for (int i = 0; i < labelN; i++) {
			fixedLabel[i] = false;
		}
		temporal_graph->readFixedLabel(fixedLabel, fixedLabelFile, FindTMotif::isEdgeTypeFixed);
	}
	else {
		for (int i = 0; i < labelN; i++) {
			fixedLabel[i] = true;
		}
	}
#pragma endregion

#pragma region run the algorithm
	switch (Setting::choice) {
	case AlgorithmType::FRTM:
	case AlgorithmType::FRTMFORDYN:
	case AlgorithmType::FRTMOPT1:
	case AlgorithmType::FRTMOPT1DYN:
	case AlgorithmType::FRTMPLUS:
	case AlgorithmType::FRTMPLUSDYN:
		for (auto iter = fixedK.begin(); iter != fixedK.end(); ++iter) {
			FindTMotif::k = *iter;
			cout << "k : " << *iter << endl;
			runStaticAlgorithm(temporal_graph,
				fixedLabel);
		}
		FindTMotif::output = 0;
		break;
	case AlgorithmType::TESTING:
		for (auto iter = fixedK.begin(); iter != fixedK.end(); ++iter) {
			FindTMotif::k = *iter;
			cout << "k : " << *iter << endl;
			testD5(temporal_graph, fixedLabel);
		}
		break;
	case AlgorithmType::DFRTM: 
	case AlgorithmType::DFRTMOPT1: 
	case AlgorithmType::DFRTMPLUS: 
		for (auto iter = fixedK.begin(); iter != fixedK.end();
			++iter) {
			FindTMotif::k = *iter;
			cout << "k : " << *iter << endl;
			runIncrementalAlgorithm(
				temporal_graph, inputSrc,
				endT, fixedLabel, newEndT);
		}
		FindTMotif::output = 0;
		break;
	default:
		cout << "please choose a algorithm to run" << endl;
		return 0;
	}
#pragma endregion

	delete temporal_graph;
	delete[] fixedLabel;
	cout << "\n";

	_CrtDumpMemoryLeaks();//check memory leak
	OUTPUT_FILE_CLOSE
}

#pragma region function implementation
#pragma region load temporal graph
void createTGraph(TGraph*& temporal_graph, const char* inputSrc, int indexId,
	int endT) {
	switch (indexId) {
	case 2://index table
		temporal_graph = DBG_NEW TGraphUDEL();
		cout << "choose Index table" << endl;
		break;
	default:
		cout << "wrong index id" << endl;
		exit(-1);
	}
	if (endT != -1) {
		temporal_graph->constructGraph(inputSrc, FindTMotif::k, endT);
	}
	else {
		temporal_graph->constructGraph(inputSrc,FindTMotif::k);
	}

	cout << *temporal_graph;
	cout << "after create index: "; Test::showMemoryUse();
	//exit(0);
}
#pragma endregion

#pragma region which problem to run

#pragma region static algorithm
void runStaticAlgorithm(TGraph*& temporal_graph, bool*& fixLabel) {
#pragma region initialzation
	int nTimestamps = temporal_graph->getCurrNTimestamp();
	#pragma region final result
	int resultLen = (nTimestamps - FindTMotif::k + 2)
		*(nTimestamps - FindTMotif::k + 1) >> 1;
	vec(TMotifII*)* res = nullptr;//used for LIMITEDFIXNUM and PROPORTION
	//vec(TMotifLM*)* res2LM = nullptr;//used for LIMITEDFIXNUM and PROPORTION
	FindTMotif::motifNumber = 0;
	res = DBG_NEW vec(TMotifII*)[resultLen];
#pragma endregion 
#pragma endregion
	Test::counter = Test::counter2 = Test::counter3 = Test::counter4 = Test::counter5 = 0;
	Test::fproc = Test::gne = Test::gm = Test::compr = 0;
	Test::peakMemory = 0;
	cout << "before algorithm: "; Test::showMemoryUse();
	int startT = temporal_graph->getStartT(), endT = temporal_graph->getEndT();
	BEGIN_TIMER(startTime)
	switch (Setting::choice) {
		case AlgorithmType::FRTM:
			FindRTMotif::FRTM(temporal_graph, res, fixLabel); 
			break;
		case AlgorithmType::FRTMFORDYN: 
			FindRTMotif::FRTMDYN(temporal_graph, res, fixLabel, startT, endT); 
			break;
		case AlgorithmType::FRTMOPT1: 
			FindRTMotif::FRTMOpt1(temporal_graph, res, fixLabel); break;
		case AlgorithmType::FRTMOPT1DYN: 
			FindRTMotif::FRTMOpt1DYN(temporal_graph, res, fixLabel, startT, endT);
			break;
		case AlgorithmType::FRTMPLUS:
			FindRTMotif::FRTMPlus(temporal_graph, res, fixLabel);
			break;
		case AlgorithmType::FRTMPLUSDYN:
			FindRTMotif::FRTMPlusDYN(temporal_graph, res, fixLabel, startT, endT);
			break;
	}
	auto allTime = END_TIMER(startTime);

	//cout << Test::counter << " " << Test::counter2 <<  " " << Test::counter3 << " " << Test::counter4 << " " << Test::counter5 << endl;
	#ifdef _USENSTIMER
		cout << "fproc: "<< Test::fproc<<"ns, compr: " << Test::compr << "ns, gm: " << Test::gm << "ns, gne: " << Test::gne << "ns" << endl;
	#else
		cout << "fproc: " << Test::fproc << "ms, compr: " << Test::compr << "ms, gm: " << Test::gm << "ms, gne: " << Test::gne << "ms" << endl;
	#endif // _USENSTIMER

	cout << "gne11: " << Test::gne11 << " gne121: " << Test::gne121 << " gne122: " << Test::gne122
		<< " gne211: " << Test::gne211 << " gne2121: " << Test::gne2121 << " gne2122: " << Test::gne2122
		<< " gne221: " << Test::gne221 << " gne222: " << Test::gne222 << " gnenonoise: " << Test::gnenonoise
		<<  " gnefield: " << Test::gnefield*1.0/(Test::gne11 + Test::gne121 + Test::gne122 + Test::gne211+ Test::gne2121 +
			Test::gne2122+ Test::gne221 + Test::gne222 + Test::gnenonoise)
		<<  " gnemaxfield: " << Test::gnemaxfield << endl;
	//cout << Test::counter*1.0 /1000 << endl;
	//cout << Test::counter6 << " "<< Test::counter3 << " " << Test::counter7 << " " << Test::counter5 << endl;
	#ifdef _USENSTIMER
		cout << "time: " << allTime << "ns" << endl;
	#else
		cout << "time: " << allTime << "ms" << endl;
	#endif // _USENSTIMER

	outputResult(temporal_graph, res, resultLen);
	releaseResult(res, resultLen);//release the memory
	
	cout << "release result: ";
	Test::showMemoryUse(); cout << "\n";
	exit(0);
}
#pragma endregion

#pragma region dynamic algorithm
void runIncrementalAlgorithm(TGraph*& temporal_graph, const char * src
	, int graphEndT, bool*& fixLabel, int newEndT) {
	if (graphEndT == -1) {
		cout << "please use -l parameter to limit the snapshot of the original graph" << endl;
		exit(0);
	}
	int allNTimestamp = temporal_graph->getAllNTimestamp();
	if (newEndT == -1) newEndT = allNTimestamp - 1;
	int listSize = newEndT - temporal_graph->getStartT() - FindTMotif::k + 2;
	int resultLen = (listSize + 1)*(listSize) >> 1;//result size
	vec(TMotifII*)* result;
	result = DBG_NEW vec(TMotifII*)[resultLen];
	
	int startT = temporal_graph->getStartT(), endT = graphEndT;
#pragma region get original result of general case
	BEGIN_TIMER(startTime)
	switch (Setting::choice) {
		case AlgorithmType::DFRTM:
			FindRTMotif::FRTMDYN(temporal_graph, result, fixLabel, 0, graphEndT);
			break;
		case AlgorithmType::DFRTMOPT1:
			FindRTMotif::FRTMOpt1DYN(temporal_graph, result, fixLabel, 0, graphEndT);
			break;
		case AlgorithmType::DFRTMPLUS:
			FindRTMotif::FRTMPlusDYN(temporal_graph, result, fixLabel, 0, graphEndT);
			break;
	}
	#ifdef _USENSTIMER
			cout << "gctime: " << END_TIMER(startTime) << "ns" << endl;
	#else
			cout << "gctime: " << END_TIMER(startTime) << "ms" << endl;
	#endif // _USENSTIMER
	temporal_graph->showMidResult(result);
#pragma endregion
	cout << "middle memory: ";
	Test::showMemoryUse(); 
#pragma region snapshots increase
	temporal_graph->setEndT(newEndT);
	temporal_graph->setCurrNTimestamp(newEndT + 1);
	cout << *temporal_graph;
	temporal_graph->updateDS(src,graphEndT, newEndT);
	temporal_graph->resetStruct();
#pragma endregion
#pragma region preprocess
	FindTMotif::motifNumber = 0;
	Test::peakMemory = 0;
	cout << "before algorithm: "; Test::showMemoryUse();
#pragma endregion
	BEGIN_TIMER(incTime)
	switch (Setting::choice) {
	case AlgorithmType::DFRTM:
		FindRTMotif::DFRTM(temporal_graph, result, fixLabel, graphEndT/*, methodid*/);
		break;
	case AlgorithmType::DFRTMOPT1:
		FindRTMotif::DFRTMOpt1(temporal_graph, result, fixLabel, graphEndT/*, methodid*/);
		break;
	case AlgorithmType::DFRTMPLUS:
		FindRTMotif::DFRTMPLus(temporal_graph, result, fixLabel, graphEndT);
		break;
	}
	#ifdef _USENSTIMER
		cout << "inctime: " << END_TIMER(incTime) << "ns" << endl;
	#else
		cout << "inctime: " << END_TIMER(incTime) << "ms" << endl;
	#endif // _USENSTIMER

	outputResult(temporal_graph, result, resultLen);
	releaseResult(result, resultLen);

	
	cout << "release result: ";
	Test::showMemoryUse(); cout << "\n";
	exit(0);
}
#pragma endregion

#pragma endregion 

#pragma region test
void testD5(TGraph*& temporal_graph, bool*& fixLabel) {
#pragma region initialzation
	int nTimestamps = temporal_graph->getCurrNTimestamp();
#pragma region final result
	int resultLen = (nTimestamps - FindTMotif::k + 2)
		*(nTimestamps - FindTMotif::k + 1) >> 1;
	long long motifNum1, motifNum2;
	if (Test::testingMode == 2) {//compare with original motifs
		vec(TMotifI*)* res = DBG_NEW vec(TMotifI*)[resultLen];//used for exact label matches
		vec(TMotifII*)* res2 = DBG_NEW vec(TMotifII*)[resultLen];//used for label mismatches
		int startT = temporal_graph->getStartT(), endT = temporal_graph->getEndT();
		FindTMotif::motifNumber = 0;
#pragma endregion 
#pragma endregion
		Test::fproc = Test::gne = Test::gm = Test::compr = 0;
		Test::peakMemory = 0;
		BEGIN_TIMER(startTime)
		FindTMotif::FRTMExact(temporal_graph, res, fixLabel, startT, endT);
		motifNum1 = FindTMotif::motifNumber;
		#ifdef _USENSTIMER
			cout << "time: " << END_TIMER(startTime) << "ns, original motif number: " << motifNum1 << endl;
		#else
			cout << "time: " << END_TIMER(startTime) << "ms, original motif number: " << motifNum1 << endl;
		#endif // _USENSTIMER

		temporal_graph->resetStruct();
		FindTMotif::motifNumber = 0;
		BEGIN_TIMER(compareTime)
		FindRTMotif::FRTMOpt1(temporal_graph, res2, fixLabel);
		motifNum2 = FindTMotif::motifNumber;
		#ifdef _USENSTIMER
			cout << "time: " << END_TIMER(compareTime) << "ns, new motif number: " << motifNum2 << endl;
		#else
			cout << "time: " << END_TIMER(compareTime) << "ms, new motif number: " << motifNum2 << endl;
		#endif // _USENSTIMER
		compareResult(temporal_graph, res, res2, resultLen, motifNum1, motifNum2);
	}
	else if (Test::testingMode == 3) {//counting overlapped
		vec(TMotifII*)* res2 = DBG_NEW vec(TMotifII*)[resultLen];//used for LIMITEDFIXNUM and PROPORTION
		int startT = temporal_graph->getStartT(), endT = temporal_graph->getEndT();
		FindTMotif::motifNumber = 0;
#pragma endregion 
#pragma endregion
		Test::counter = Test::counter2 = Test::counter3 = Test::counter4 = Test::counter5 = 0;
		Test::fproc = Test::gne = Test::gm = Test::compr = 0;
		Test::peakMemory = 0;
		FindTMotif::motifNumber = 0;
		Setting::c = 0x7fffffff;//not limited
		FindRTMotif::FRTMOpt1(temporal_graph, res2, fixLabel); 
		motifNum2 = FindTMotif::motifNumber;
		cout << "motif number: " << FindTMotif::motifNumber << endl;
		//countingResult(temporal_graph, res, resultLen, motifNum2);
	}
	else if (Test::testingMode == 4) {//compare with c = inf
		vec(TMotifII*)* res = DBG_NEW vec(TMotifII*)[resultLen];//used for precise method and GLOBALFIXNUM
		vec(TMotifII*)* res2 = DBG_NEW vec(TMotifII*)[resultLen];//used for LIMITEDFIXNUM and PROPORTION
		int startT = temporal_graph->getStartT(), endT = temporal_graph->getEndT();
		FindTMotif::motifNumber = 0;
#pragma endregion 
#pragma endregion
		Test::fproc = Test::gne = Test::gm = Test::compr = 0;
		Test::peakMemory = 0;
		BEGIN_TIMER(startTime)
		int tempC = Setting::c;
		Setting::c = 0x7fffffff;//not limited
		FindRTMotif::FRTMPlus(temporal_graph, res, fixLabel);
		motifNum1 = FindTMotif::motifNumber;
		#ifdef _USENSTIMER
			cout << "time: " << END_TIMER(startTime) << "ns, original motif number: " << motifNum1 << endl;
		#else
			cout << "time: " << END_TIMER(startTime) << "ms, original motif number: " << motifNum1 << endl;
		#endif // _USENSTIMER
		temporal_graph->resetStruct();
		FindTMotif::motifNumber = 0;
		BEGIN_TIMER(compareTime)
		Setting::c = tempC;
		FindRTMotif::FRTMPlus(temporal_graph, res2, fixLabel);
		motifNum2 = FindTMotif::motifNumber;
		#ifdef _USENSTIMER
			cout << "time: " << END_TIMER(compareTime) << "ns, new motif number: " << motifNum2 << endl;
		#else
			cout << "time: " << END_TIMER(compareTime) << "ms, new motif number: " << motifNum2 << endl;
		#endif // _USENSTIMER
		compareResult(temporal_graph, res, res2, resultLen, motifNum1, motifNum2);
	}
	cout << "\n";
	exit(0);
}

#pragma endregion


#pragma region output
//void outputResult(TGraph*& temporal_graph, vec(TMotifI*)*& result, int resultLen) {
//	veciter(TMotifI*) iter;
//	TMotifI* motif;
//	size_t size;
//	FindTMotif::motifMaxNum = FindTMotif::motifSum = 0;
//	Test::maxIntvLen = 0, Test::sumIntvLen = 0;
//	for (int i = 0; i < resultLen; i++) {
//		size = result[i].size();
//		if (size == 0) continue;
//		iter = result[i].begin();
//		motif = (TMotifI*) *iter;
//		//if (motif->getEndT() - motif->getStartT() + 1 < FindTMotif::filterIntvSize) continue;
//		/*cout << "startT: " << motif->getStartT()
//			<< "\tendT: " << motif->getEndT() << endl;
//		cout << MOTIF_NUM << size << endl;
//		*/
//		//if (motif->getStartT() == 0 && motif->getEndT() == 119 )FindTMotif::output = 2;
//		if (FindTMotif::output >= 1) {
//			cout << MOTIF_NUM << size << endl;
//			cout << "startT: " << motif->getStartT()
//				<< "\tendT: " << motif->getEndT() << endl;
//		}
//		FindTMotif::print(temporal_graph,
//			result, i, false);
//		//if (motif->getStartT() == 0 && motif->getEndT() == 119) { FindTMotif::output = 0; exit(0); }
//	}
//
//	cout << "motif max edges num: " << FindTMotif::motifMaxNum <<
//		"\tmotif avg edges num: " << FindTMotif::motifSum * 1.0 / FindTMotif::motifNumber << endl;
//	cout << "motif max intverval length: " << Test::maxIntvLen <<
//		"\tmotif avg intverval length: " << Test::sumIntvLen * 1.0 / FindTMotif::motifNumber << endl;
//	cout << "sum: " << FindTMotif::motifNumber << endl;
//
//	//cout << "counter: " << Test::counter << endl;
//	/*if (FindTMotif::output == 3) {
//		cout << "-------------------------------------------------------" << endl;
//		for (int i = 0; i < resultLen; i++) {
//			size = result[i].size();
//			if (size == 0) continue;
//			iter = result[i].begin();
//			motif = *iter;
//			if (motif->getEndT() - motif->getStartT() + 1 < FindTMotif::filterIntvSize) continue;
//			cout << MOTIF_NUM << size << endl;
//			cout << "startT: " << motif->getStartT()
//				<< "\tendT: " << motif->getEndT() << endl;
//			FindTMotif::printRealLabel(temporal_graph, result[i]);
//		}
//	}*/
//
//	cout << "after algorithm: "; Test::showPeakMemoryUse();
//	Test::showRealPeakMemoryUse();
//
//	/*if (Test::testingMode) {
//		if (Test::globalApproType && Test::globalApproType != 1)
//			testOverlap(temporal_graph, res);
//		if(Test::globalApproType && (Test::globalApproType == 1 || Test::globalApproType == 3))
//			testLabel(temporal_graph, res);
//		showTesting();
//	}*/
//}
void outputResult(TGraph*& temporal_graph, vec(TMotifII*)*& result, int resultLen) {
	veciter(TMotifII*) iter;
	TMotifII* motif;
	size_t size;
	FindTMotif::motifMaxNum = FindTMotif::motifSum = 0;
	Test::maxIntvLen = 0, Test::sumIntvLen = 0;
	for (int i = 0; i < resultLen; i++) {
		size = result[i].size();
		if (size == 0) continue;
		iter = result[i].begin();
		motif = (TMotifII*) *iter;
		//if (motif->getEndT() - motif->getStartT() + 1 < FindTMotif::filterIntvSize) continue;
		/*cout << "startT: " << motif->getStartT()
			<< "\tendT: " << motif->getEndT() << endl;
		cout << MOTIF_NUM << size << endl;
		*/
		//if (motif->getStartT() == 736 && motif->getEndT() == 801)FindTMotif::output = 2;

		if (FindTMotif::output >= 1) {
			cout << MOTIF_NUM << size << endl;
			cout << "startT: " << motif->getStartT()
				<< "\tendT: " << motif->getEndT() << endl;
		}
		FindTMotif::print(temporal_graph,
			result, i, false);
		//if (motif->getStartT() == 736 && motif->getEndT() == 801) { FindTMotif::output = 0; }
		
		//if(motif->getStartT() == 1 )exit(0);
		
		//testing
		/*if (Test::testingMode == 1 && Setting::globalApproType == 5) {
			testEndpointD5(temporal_graph, result[i]);
			if (temporal_graph->getCurrNTimestamp()>2000) cout << "Complete testing: " << i <<" / "<< resultLen << endl;
		}*/
	}
	/*if (Test::testingMode == 1) {
		cout << "Complete testing" << endl;
	}*/
	cout << "motif max edges num: " << FindTMotif::motifMaxNum <<
		"\tmotif all edges num: " << FindTMotif::motifSum << endl;
	cout << "motif max intverval length: " << Test::maxIntvLen <<
		"\tmotif avg intverval length: " << Test::sumIntvLen * 1.0 / FindTMotif::motifNumber << endl;
	cout << "sum: " << FindTMotif::motifNumber << endl;

	/*if (FindTMotif::output == 3) {
		cout << "-------------------------------------------------------" << endl;
		for (int i = 0; i < resultLen; i++) {
			size = result[i].size();
			if (size == 0) continue;
			iter = result[i].begin();
			motif = *iter;
			if (motif->getEndT() - motif->getStartT() + 1 < FindTMotif::filterIntvSize) continue;
			cout << MOTIF_NUM << size << endl;
			cout << "startT: " << motif->getStartT()
				<< "\tendT: " << motif->getEndT() << endl;
			FindTMotif::printRealLabel(temporal_graph, result[i]);
		}
	}*/

	cout << "after algorithm: "; Test::showPeakMemoryUse();
	Test::showRealPeakMemoryUse();
}

void compareResult(TGraph*& temporal_graph, vec(TMotifI*)*& result1, vec(TMotifII*)*& result2, int resultLen,long long motifNum1, long long motifNum2) {
	veciter(TMotifI*) iter1, iter1End;
	veciter(TMotifII*) iter2, iter2End, newiter1,newiter1End;
	TMotifI* motif1;
	TMotifII* motif2, *newmotif1;
	size_t size1, size2;
	int graphStartT = temporal_graph->getStartT(), graphEndT = temporal_graph->getEndT();
	FindTMotif::motifMaxNum = FindTMotif::motifSum = 0;
	Test::maxIntvLen = 0, Test::sumIntvLen = 0;
	EdgeIdArray* edgesList2, *newedgesList1;
	vec(int) edgesList1;
	vec(TMotifII*)* res = DBG_NEW vec(TMotifII*)[resultLen];
	int allcontainedmotif = 0, allcontainmotif = 0, allcontain1motif = 0;
	for (int i = 0; i < resultLen; i++) {
		size1 = result1[i].size();
		if (size1 == 0) continue;
		iter1 = result1[i].begin();
		iter1End = result1[i].end();
		for (; iter1 != iter1End; ++iter1) {
			motif1 = (TMotifI*)*iter1;

			edgesList1.clear();
			temporal_graph->getMotifEdges(result1, motif1, edgesList1, FindTMotif::k);
			sort(edgesList1.begin(), edgesList1.end());
			res[i].emplace_back(DBG_NEW TMotifII(edgesList1, motif1->getStartT(), motif1->getEndT()));

			delete motif1;
		}
		result1[i].clear();
	}
	delete[] result1;

	for (int i = 0; i < resultLen; i++) {
		size2 = result2[i].size();
		if (size2 == 0) continue;
		iter2 = result2[i].begin();
		iter2End = result2[i].end();
		for (int num2 = 1; iter2 != iter2End; ++iter2,++num2) {
			motif2 = (TMotifII*)*iter2;

			//if (motif2->getSize() > 50) continue;//filter for case study

			motif2->sortEdges();
			edgesList2 = motif2->getMotifEdge();

			auto edgesBegin2 = edgesList2->begin(), edgesEnd2 = edgesList2->end();
			int startT = motif2->getStartT(), endT = motif2->getEndT();
			int stopST = endT - FindTMotif::k + 1;
			int contain = 0;

			cout << "#["<<startT<< ","<< endT<< "]:"<<num2 << " -> "<<endl;
			temporal_graph->printMotif(motif2, num2);

			for (int st = startT; st <= stopST; ++st) {
				int resultPos = resultPos(st, st + FindTMotif::k - 1, graphStartT, graphEndT, FindTMotif::k);
				for (int et = st + FindTMotif::k - 1; et <= endT; ++et, ++resultPos) {
					size1 = res[resultPos].size();
					if (size1 == 0) continue;
					newiter1 = res[resultPos].begin();
					newiter1End = res[resultPos].end();
					for (int num1 = 1; newiter1 != newiter1End; ++newiter1, ++num1) {
						newmotif1 = (TMotifII*)*newiter1;
						//if (newmotif1->getSize() > 50) continue;//filter for case study
						newedgesList1 = newmotif1->getMotifEdge();
						
						auto edgeiter1 = newedgesList1->begin(), edgesEnd1 = newedgesList1->end();
						auto edgeiter2 = edgesBegin2;
						long long edgesNum1 = newedgesList1->size(), edgesNum2 = edgesList2->size();
						while(edgeiter1 != edgesEnd1 && edgeiter2 != edgesEnd2) {
							if (*edgeiter1 == *edgeiter2) {
								edgeiter1++;
								edgesNum1--;
							}
							edgeiter2++;
							edgesNum2--;
							if (edgesNum1 > edgesNum2)
								break;
						}
						if (edgeiter1 == edgesEnd1) {//containment
							contain++;
							//if (edgeiter2 != edgesEnd2) {//not same
							cout << "$[" << st << "," << et << "]:" << num1 << " -> " << endl;
							temporal_graph->printMotif(newmotif1, num1);
							//}
							break;
						}
					}
				}
			}
			if (contain != 0) {
				allcontainmotif++;
				if (contain == 1) allcontain1motif++;
				allcontainedmotif += contain;
			}
			cout << "@contain " << contain << " motifs" << endl;
		}
	}
	cout << "new motifs have " << allcontainmotif << "/" << motifNum2 << "=" << allcontainmotif * 1.0 / motifNum2 << " which contains original motifs, each motif contains average " << allcontainedmotif * 1.0 / motifNum2<<" motifs, new motifs have " << allcontain1motif << " / " << motifNum2 << " = " << allcontain1motif * 1.0 / motifNum2 << " which contains only 1 original motif"<<endl;
	cout << "original motifs have " << allcontainedmotif << "/" << motifNum1 << "=" << allcontainedmotif * 1.0 / motifNum1 << " which are contained" << endl;
}

void compareResult(TGraph*& temporal_graph, vec(TMotifII*)*& result1, vec(TMotifII*)*& result2, int resultLen,long long motifNum1, long long motifNum2) {
	veciter(TMotifII*) iter1, iter1End;
	veciter(TMotifII*) iter2, iter2End, newiter1,newiter1End;
	TMotifII* motif1;
	TMotifII* motif2, *newmotif1;
	size_t size1, size2;
	int graphStartT = temporal_graph->getStartT(), graphEndT = temporal_graph->getEndT();
	FindTMotif::motifMaxNum = FindTMotif::motifSum = 0;
	Test::maxIntvLen = 0, Test::sumIntvLen = 0;
	EdgeIdArray* edgesList2, *newedgesList1;
	vec(int) edgesList1;
	vec(TMotifII*)* res = DBG_NEW vec(TMotifII*)[resultLen];
	int allcontainedmotif = 0, allcontainmotif = 0, allcontain1motif = 0;
	for (int i = 0; i < resultLen; i++) {
		size1 = result1[i].size();
		if (size1 == 0) continue;
		iter1 = result1[i].begin();
		iter1End = result1[i].end();
		for (; iter1 != iter1End; ++iter1) {
			motif1 = (TMotifII*)*iter1;

			edgesList1.clear();
			edgesList1.insert(edgesList1.begin(),motif1->getMotifEdge()->begin(), motif1->getMotifEdge()->end());
			sort(edgesList1.begin(), edgesList1.end());
			res[i].emplace_back(DBG_NEW TMotifII(edgesList1, motif1->getStartT(), motif1->getEndT()));

			delete motif1;
		}
		result1[i].clear();
	}
	delete[] result1;

	for (int i = 0; i < resultLen; i++) {
		size2 = result2[i].size();
		if (size2 == 0) continue;
		iter2 = result2[i].begin();
		iter2End = result2[i].end();
		for (int num2 = 1; iter2 != iter2End; ++iter2,++num2) {
			motif2 = (TMotifII*)*iter2;
			motif2->sortEdges();
			edgesList2 = motif2->getMotifEdge();

			auto edgesBegin2 = edgesList2->begin(), edgesEnd2 = edgesList2->end();
			int startT = motif2->getStartT(), endT = motif2->getEndT();
			int stopST = endT - FindTMotif::k + 1;
			int contain = 0;

			cout << "#["<<startT<< ","<< endT<< "]:"<<num2 << " -> "<<endl;
			
			for (int st = startT; st <= stopST; ++st) {
				int resultPos = resultPos(st, st + FindTMotif::k - 1, graphStartT, graphEndT, FindTMotif::k);
				for (int et = st + FindTMotif::k - 1; et <= endT; ++et, ++resultPos) {
					size1 = res[resultPos].size();
					if (size1 == 0) continue;
					newiter1 = res[resultPos].begin();
					newiter1End = res[resultPos].end();
					for (int num1 = 1; newiter1 != newiter1End; ++newiter1, ++num1) {
						newmotif1 = (TMotifII*)*newiter1;
						newedgesList1 = newmotif1->getMotifEdge();
						
						auto edgeiter1 = newedgesList1->begin(), edgesEnd1 = newedgesList1->end();
						auto edgeiter2 = edgesBegin2;
						long long edgesNum1 = newedgesList1->size(), edgesNum2 = edgesList2->size();
						while(edgeiter1 != edgesEnd1 && edgeiter2 != edgesEnd2) {
							if (*edgeiter1 == *edgeiter2) {
								edgeiter1++;
								edgesNum1--;
							}
							edgeiter2++;
							edgesNum2--;
							if (edgesNum1 > edgesNum2)
								break;
						}
						if (edgeiter1 == edgesEnd1) {//containment
							contain++;
							//if (edgeiter2 != edgesEnd2) {//not same
							cout << "$[" << st << "," << et << "]:" << num1 << endl;
							//}
							break;
						}
					}
				}
			}
			if (contain != 0) {
				allcontainmotif++;
				if (contain == 1) allcontain1motif++;
				allcontainedmotif += contain;
			}
			cout << "@contain " << contain << " motifs" << endl;
		}
	}
	cout << "new motifs have " << allcontainmotif << "/" << motifNum2 << "=" << allcontainmotif * 1.0 / motifNum2 << " which contains original motifs, each motif contains average " << allcontainedmotif * 1.0 / motifNum2<<" motifs, new motifs have " << allcontain1motif << " / " << motifNum2 << " = " << allcontain1motif * 1.0 / motifNum2 << " which contains only 1 original motif"<<endl;
	cout << "original motifs have " << allcontainedmotif << "/" << motifNum1 << "=" << allcontainedmotif * 1.0 / motifNum1 << " which are contained" << endl;
}

void countingResult(TGraph*& temporal_graph, vec(TMotifII*)*& result2, int resultLen, long long motifNum2) {
	veciter(TMotifII*) iter1, iter1End, iter2, iter2End;
	TMotifII* motif1, *motif2;
	size_t size1, size2;
	FindTMotif::motifMaxNum = FindTMotif::motifSum = 0;
	Test::maxIntvLen = 0, Test::sumIntvLen = 0;
	EdgeIdArray* edgesList1, *edgesList2;
	int graphStartT = temporal_graph->getStartT(), graphEndT = temporal_graph->getEndT(); 
	int allcontainedmotif = 0, allcontainmotif = 0;
	for (int i = 0; i < resultLen; i++) {
		size2 = result2[i].size();
		if (size2 == 0) continue;
		iter2 = result2[i].begin();
		iter2End = result2[i].end();
		for (; iter2 != iter2End; ++iter2) {
			motif2 = (TMotifII*)*iter2;
			motif2->sortEdges();
		}
	}	
	for (int i = 0; i < resultLen; i++) {
		size2 = result2[i].size();
		if (size2 == 0) continue;
		iter2 = result2[i].begin();
		iter2End = result2[i].end();
		for (int num2 = 1; iter2 != iter2End; ++iter2, ++num2) {
			motif2 = (TMotifII*)*iter2;
			edgesList2 = motif2->getMotifEdge();
			auto edgesBegin2 = edgesList2->begin(), edgesEnd2 = edgesList2->end();
			int startT = motif2->getStartT(), endT = motif2->getEndT();
			int snapshot2 = endT - startT + 1;
			double gap = snapshot2 * Setting::delta/(1 - Setting::delta);
			if (LESS(gap,FindTMotif::k)) break;
			edgesList2 = motif2->getMotifEdge();
			int stopST = endT - FindTMotif::k;
			int contain = 0;

			cout << "#[" << startT << "," << endT << "]:" << num2 << " -> " << endl;

			for (int st = startT+1; st <= stopST; ++st) {
				int resultPos = resultPos(st, st + FindTMotif::k - 1, graphStartT, graphEndT, FindTMotif::k);
				for (int et = st + FindTMotif::k - 1; et < endT; ++et, ++resultPos) {
					if (LESSEQ(et - st + 1, gap)) continue;
					size1 = result2[resultPos].size();
					if (size1 == 0) continue;
					iter1 = result2[resultPos].begin();
					iter1End = result2[resultPos].end();
					for (int num1 = 1; iter1 != iter1End; ++iter1, ++num1) {
						motif1 = (TMotifII*)*iter1;
						edgesList1 = motif1->getMotifEdge();

						auto edgeiter1 = edgesList1->begin(), edgesEnd1 = edgesList1->end();
						auto edgeiter2 = edgesBegin2;
						long long edgesNum1 = edgesList1->size(), edgesNum2 = edgesList2->size();
						while (edgeiter1 != edgesEnd1 && edgeiter2 != edgesEnd2) {
							if (*edgeiter1 == *edgeiter2) {
								edgeiter1++;
								edgesNum1--;
							}
							edgeiter2++;
							edgesNum2--;
							if (edgesNum1 > edgesNum2)
								break;
						}
						if (edgeiter1 == edgesEnd1) {//containment
							contain++;
							cout << "$[" << st << "," << et << "]:" << num1 << endl;
							break;
						}
					}
				}
			}
			if (contain != 0) {
				allcontainmotif++;
				allcontainedmotif += contain;
			}
			cout << "@contain " << contain << " motifs" << endl;
		}
	}
	cout << "motifs have " << allcontainmotif << "/" << motifNum2 << "=" << allcontainmotif * 1.0 / motifNum2 << " which contains other motifs, each motif contains average " << allcontainedmotif*1.0 / motifNum2 << " motifs" << endl;
	cout << "edges have "<< Test::allContainment << "/" << Test::selectedSum << "=" << Test::allContainment * 1.0 / Test::selectedSum << " intervals which contain other intervals" << endl;
}


//print the peak working set size in bytes
void showPeakMemoryUse() {
	//get the handle of the current process
	HANDLE currentProcess = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
	cout << "Peak Memory Use:" << fixed << setprecision(2) <<
		pmc.PeakWorkingSetSize*1.0 << "Byte" << endl;
}

//print the working set size in bytes
void showMemoryUse() {
	//get the handle of the current process
	HANDLE currentProcess = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
	cout << "Memory Use:" <<
		pmc.WorkingSetSize*1.0 << "Byte" << endl;
}
#pragma endregion 

#pragma region read

//load the setting of frequent condition k
void readK(vec(int)& arr, const char* file) {
	FILE* f;
	f = fopen(file, "r+");
	if (!f) {
		arr.emplace_back(DEFAULT_K);
	}
	else {
		char line[LINE_LENGTH];
		CLEARALL(line, 0, LINE_LENGTH, char);
		while (fgets(line, LINE_LENGTH, f)) {
			//FindTMotif::motifNumber = 0;
			if (strlen(line) == 0) continue;
			arr.emplace_back(STR2INT(line));
		}
		fclose(f);
	}
}


#pragma endregion

#pragma region release memory
void releaseResult(vec(TMotifI*)*& result, int len) {
	size_t size;
	for (int i = 0; i < len; i++) {
		size = result[i].size();
		for (size_t j = 0; j < size; j++) {
			if (result[i][j] != nullptr)
				delete result[i][j];
		}
		result[i].clear();
	}
	delete[] result;
}
void releaseResult(vec(TMotifII*)*& result, int len) {
	size_t size;
	for (int i = 0; i < len; i++) {
		size = result[i].size();
		for (size_t j = 0; j < size; j++) {
			if (result[i][j] != nullptr)
				delete result[i][j];
		}
		result[i].clear();
	}
	delete[] result;
}
#pragma endregion
#pragma endregion 


/*test whether motifs have all edges with same label in endpoints of interval
used for PROPORTION*/
void testEndpointD5(TGraph*& temporalgraph, vec(TMotifII*)& lis) {
	veciter(TMotifII*) listIter = lis.begin(),
		listEnd = lis.end();
	TMotifII* motif;
	int i = 0;
	for (; listIter != listEnd; ++listIter, ++i) {
		motif = *listIter;
		temporalgraph->checkMotifEndpointsD5(motif);
	}
}
