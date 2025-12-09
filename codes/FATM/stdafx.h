#pragma once

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#else
#define DBG_NEW new
#endif

#include <stdio.h>
#include<iostream>
#include<cctype>
#include<cstring>
#include<algorithm>
#include<vector>
#include<queue>
#include<limits>
#include<exception>
#include<unordered_map>
#include<unordered_set>
#include<string>
#include<fstream>
#include<sstream>
#include<stack>
#include<map>
#include<set>
#include<ctime>
#include<cstdlib>
#include<Windows.h>
#include<psapi.h>
#include<iomanip>
#include<functional>
#include <memory>
#include <cassert>
#include <chrono>
using namespace std;

#define EPS 1e-7
#define EQ(a,b) ((fabs((a) - (b)))<(EPS))
#define MORE(a,b) (((a) - (b)) > (EPS))
#define MOREEQ(a,b) (((a) - (b)) > (-EPS))
#define LESS(a,b) (((b) - (a)) > (EPS))
#define LESSEQ(a,b) (((a) - (b)) < (EPS))

#pragma region default input parameter
#define DEFAULT_K 10 //default frequent condition k
#define FILE_NAME_LENGTH 100 //maximum length of file name
#define OUTPUT_FILE "example-output.txt" //default output file name
#define INPUT_FILE "example-temporal-graph.txt" //default input temporal graph file name
#define K_FILE "example-k.txt" //default file name for frequent condition k
#define FIXEDLABEL_FILE "example-fixedlabel.txt" //default file name for fixed label
#define ALGORITHM_ID 1 //default algorithm id
#define RUNTIMES 1 //default times of running the algorithm
#define INDEXID 2 //default index id (del table)
#pragma endregion

#pragma region util function
#define CLEARALL(a,value,num,type) memset(a,value,sizeof(type)*num)
#define OUTPUT_FILE_OPEN freopen(FindTMotif::outputSrc, "a", stdout); ios::sync_with_stdio(false);
#define OUTPUT_FILE_CLOSE fclose(stdout);
#define VEC_RELEASE(type,vec) vector<type>().swap(vec)
#pragma endregion

#pragma region output parameter
#define EDGE_NUM "edges: "
#define INTV_NUM "maxVe: "
#define MEAN_INTV_NUM "averVe: "
#define MOTIF_ID "No: "
#define MIDRESULT_MEMORY "Test intermediate result memory:\n"
#define RUN_TIMES "run times: "
#define MOTIF_NUM "motifs: "
#define SELECT_EDGE "maxSm: "
#define MEAN_SELECT_EDGE "averSm: "
#pragma endregion

#pragma region testing
#define TIMES_PER_SEC (1.0e9)
#pragma endregion


#pragma region type
#define vec(a) vector<a>
#define veciter(a) vector<a>::iterator
#define vecriter(a) vector<a>::reverse_iterator

using ibPairVec = vector<pair<int, bool>>;
using ibPairVec_Iter = vector<pair<int, bool>>::iterator;

using i2tupHMap = unordered_map<int, tuple<int,int,bool>>;
using i2tupHMap_Iter = unordered_map<int, tuple<int, int, bool>>::iterator;
using i2iHMap = unordered_map<int, int>;
using iSet = unordered_set<int>;
using i2bHMap = unordered_map<int, bool>;
using i2iHMap_Iter = unordered_map<int, int>::iterator;
using i2bHMap_Iter = unordered_map<int, bool>::iterator;
#pragma endregion

//information saved for R edge sets and maxIntv
struct SaveInfo {
	int edgeId; //edge in R edge sets
	int startT; //interval for maxIntv[e], only need to save the starting timestamp
	int labelid; //label of e in interval
	SaveInfo(int id, int start, int w) {
		edgeId = id;
		startT = start;
		labelid = w;
	}
};
using SAVEINFO_Vec = vec(SaveInfo);
using SAVEINFO_VIter = veciter(SaveInfo);
