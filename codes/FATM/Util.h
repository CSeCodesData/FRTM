#pragma once
#include "stdafx.h"
#include "Strategy.h"
#define LINE_LENGTH 1000
#define SEP_CHAR ','
#define STR2INT(a) Util::stringToInt(a)
#define STR2DOU(a) Util::stringToDouble(a)
#define STR2BOOL(a) Util::stringToBool(a)
#define EXIT exit(-1);
//#define _USENSTIMER

#ifdef _USENSTIMER
	#define BEGIN_TIMER(a) auto a = std::chrono::steady_clock::now();
	#define END_TIMER(a) std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count()
	#define OUTPUT_TIMER(a) cout<<END_TIMER(a)<<endl;
#else
	#define BEGIN_TIMER(a) auto a = clock();
	#define END_TIMER(a) clock() - a
	#define OUTPUT_TIMER(a) cout<<END_TIMER(a)<<endl;
#endif // _USENSTIMER


#define NEWLINE cout << endl;
#define OUTPUT(...) Util::output(__VA_ARGS__); NEWLINE

#define MYINFINITE 0x3f3f3f3f

#define INTV(a,b) cout<<"("<<a<<", "<<b<<") ";

using namespace std;
class Util {
public:
	//transform char* to double
	inline static double stringToDouble(const char* str) {
		return atof(str);
	}

	//transform char* to int
	inline static int stringToInt(const char* str) {
		return atoi(str);
	}

	//transform char* to bool
	inline static bool stringToBool(const char* str) {
		return 0!=atoi(str);
	}
	
	//transform long to string
	inline static string longToString(const long long num) {
		ostringstream os;
		os << num;
		string result;
		istringstream is(os.str());
		is >> result;
		return result;
	}

	
	//print exception
	inline static void printError(const char* hint) {
		cout << hint << endl;
		EXIT
	}
	
	//return the maximum number among a,b,c
	inline static int getMax(int a,int b,int c) {
		int max = a > b ? a : b;
		return max > c ? max : c;
	}
	
	//return the minimum number among a,b,c
	inline static int getMin(int a, int b, int c) {
		int min = a < b ? a : b;
		return min < c ? min : c;
	}

	//whether item is in arr
	static bool findItem(vec(int)& arr, int item) {
		if (arr.size() == 0) return false;
		veciter(int) arrEnd = arr.end();
		for (auto arrIter = arr.begin(); arrIter != arrEnd; ++arrIter) {
			if (*arrIter == item) {
				return true;
			}
		}
		return false;
	}

	template<typename T>
	static void releaseVec(vector<T>*& arr) {
		auto arrEnd = arr->end();
		for (auto iter = arr->begin(); iter != arrEnd; ++iter) {
			delete *iter;
		}
		delete arr;
	}

	//output
	template<typename T, typename... Args>
	static void output(const T& t) {
		cout << t << " ";
	}

	//output
	template<typename T>
	static void printVec(vector<T>& t) {
		auto iterEnd = t.end();
		for (auto iter = t.begin(); iter != iterEnd; ++iter) {
			cout << *iter << endl;
		}
	}

	template<typename T, typename... Args>
	static void output(const T& t, const Args&... rest) {
		cout << setprecision(4) << fixed<< t << " ";
		output(rest...);
	}
};

#define MYERROR(a) Util::printError(a);
#define FILE_N_E MYERROR("file not exists")
#define WRONG_A_P MYERROR("wrong approximation parameter")
#define LOAD_ERROR MYERROR("load error")

#pragma region testing
class Test {
public:
	//static unordered_map<int, int> hist;
	//static unordered_map<int, int> vertexHist;

	static long long fproc;
	static long long compr, gm, gne;
	static long long ckf, gmni, gmli;

	static long long gne11, gne121, gne122;
	static long long gne211, gne2121, gne2122, gne221, gne222;
	static long long gnenonoise;
	static long long gnefield, gnemaxfield;

	static int testingMode;//testing mode

	static long long allContainment, containment, selectedSum;

	static clock_t msTimer;
	
	static long long sumIntvLen, maxIntvLen;

	static long long counter, counter2, counter3, counter4, counter5, counter6, counter7, counter8;
	static double counter9, counter10;

	static SIZE_T peakMemory;

	static void updateMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		if (peakMemory < pmc.WorkingSetSize) peakMemory = pmc.WorkingSetSize;

		if (peakMemory / 1073741824 > 64) {
			cout << "out of memory" << endl;
			exit(0);
		}
	}

	static void showIncMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		if (peakMemory < pmc.WorkingSetSize) { 
			cout << "Inc Peak Memory Use:" <<
				pmc.WorkingSetSize - peakMemory << "Byte" << endl;
			peakMemory = pmc.WorkingSetSize;
		}
	}

	static void showRealPeakMemoryUse() {
		cout << "Real Peak Memory Use:" <<
			peakMemory << "Byte" << endl;
	}

	//print the working set size in bytes
	static void showMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		cout << "Memory Use:" <<
			pmc.WorkingSetSize << "Byte" << endl;
	}

	//print the peak working set size in bytes
	static void showPeakMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		cout << "Peak Memory Use:" <<
			pmc.PeakWorkingSetSize << "Byte" << endl;
	}
};

#define TIMES_PER_SEC (1.0e9)
#pragma endregion

