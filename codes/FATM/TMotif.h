#pragma once
#include "stdafx.h"
#include "DisjointSet.h"
#include "LinkedList.h"
#include "Edge.h"

struct QueueItem {
	int startTime, endTime, id;
	QueueItem(int s, int e, int p) :startTime(s), endTime(e), id(p) {}
	bool operator < (const QueueItem& b) const {
		return this->endTime < b.endTime || (this->endTime == b.endTime && this->id < b.id);
	}
};

#pragma region connected components
struct SaveCCInfo {
	int savePos;//cc's saving position in the result
	int saveEndTime;//the end time of the interval of cc
	int saveCCInfoPos;//cc's saveCCInfo position
	SaveCCInfo(int saveP, int saveET) {
		savePos = saveP;
		saveEndTime = saveET;
		saveCCInfoPos = -1;
	}
	SaveCCInfo(int saveP, int saveET, int saveCCInfoP) {
		savePos = saveP;
		saveEndTime = saveET;
		saveCCInfoPos = saveCCInfoP;
	}

	SaveCCInfo() {
		savePos = saveEndTime = saveCCInfoPos = -1;
	}
};
/*used in precise method and GLOBALFIXNUM*/
class CComponents {
public:
	vec(TEdge) edges;//cc's edges with labels
	vec(SaveCCInfo) saveInfo;//cc's saving information in the result
	int startT;//the starting timestamp of cc.intvl
	int root;//the root of cc in disjoint set
	//bool overlapFlag;//whether overlap with other components
	//int LastNoisePos;// the condition when motif is left expandable
	CComponents(int start, int root) :startT(start),
		root(root)/*, overlapFlag(true)*/{}
	CComponents() = default;
	~CComponents() = default;
};

class CComponentsShortIntv {
public:
	vec(int) edges;//cc's edges with labels
	int startT;//the starting timestamp of cc.intvl
	int root;//the root of cc in disjoint set
	int scopeL, scopeR;// scope[e] for each e in cc
	bool haveNonExpand;//whether has non-expandable edges

	CComponentsShortIntv(int start, int root) :startT(start),
		root(root), haveNonExpand(false){}
	CComponentsShortIntv() = default;
	~CComponentsShortIntv() = default;


	virtual void combineFrom(CComponentsShortIntv*& other) {
		vec(int)& tempEdges = other->edges;//edges of currentCComponents
		int size = (int)edges.size();
		edges.resize(size + tempEdges.size());
		copy(tempEdges.begin(), tempEdges.end(), edges.begin() + size);

		if (startT < other->startT) {
			startT = other->startT;
		}
		haveNonExpand |= other->haveNonExpand;
		scopeL = max(scopeL, other->scopeL);
		scopeR = min(scopeR, other->scopeR);
	}
};


class CComponentsII {
public:
	vec(int) edges;//cc's edges
	int root;//the root of cc in disjoint set
	int newInsert;//the start position of new inserted edges
	
	CComponentsII(int root) : root(root) {
		newInsert = 0;
	}
	
	CComponentsII() { 
		newInsert = 0;
	}
	virtual ~CComponentsII() {
	}
};


struct NoisePos {
	int startTime;
	int pos;
	NoisePos(int startT, int p) {
		startTime = startT;
		pos = p;
	}
	bool operator<(const NoisePos & nosieP) const {
		return (startTime < nosieP.startTime ||
			(startTime == nosieP.startTime && pos < nosieP.pos));
	}
};
/*cc set for FRTM*/
class CComponentsFRTM : public CComponentsII {
public:
	CComponentsFRTM(int root) :CComponentsII(root){
		subCCs = nullptr;
		subCCsaved = nullptr;
	}
	CComponentsFRTM() :CComponentsII() {
		subCCs = nullptr;
		subCCsaved = nullptr;
	}
	virtual ~CComponentsFRTM() {
		if (subCCsaved != nullptr) delete[] subCCsaved;
		if (subCCs != nullptr) delete[] subCCs;
	}

	virtual void combineFrom(CComponentsFRTM*& other) {
		vec(int)& tempEdges = other->edges;//edges of currentCComponents
		int size = (int)edges.size();
		edges.resize(size + tempEdges.size());
		copy(tempEdges.begin(), tempEdges.end(), edges.begin() + size);

		while (!other->noisePosQueue.empty()) {//combine noise position
			auto noiseEdge = other->noisePosQueue.top();
			noisePosQueue.emplace(noiseEdge.first, NoisePos(noiseEdge.second.startTime, noiseEdge.second.pos + size));
			other->noisePosQueue.pop();
		}
		
		while (!other->preEMaxIntvlChangePos.empty()) {
			auto EMaxIntvlChangePos = other->preEMaxIntvlChangePos.top();
			preEMaxIntvlChangePos.emplace(EMaxIntvlChangePos.startTime, EMaxIntvlChangePos.endTime, EMaxIntvlChangePos.id);
			other->preEMaxIntvlChangePos.pop();
		}

		while (!other->maxEMaxIntvlStartTQueue.empty()) {
			auto item = other->maxEMaxIntvlStartTQueue.top();
			maxEMaxIntvlStartTQueue.emplace(item);
			other->maxEMaxIntvlStartTQueue.pop();
		}

		minEMaxIntvlEndT = min(other->minEMaxIntvlEndT, minEMaxIntvlEndT);
		/*while (!other->minEMaxIntvlEndTQueue.empty()) {
			auto item = other->minEMaxIntvlEndTQueue.top();
			minEMaxIntvlEndTQueue.emplace(item);
			other->minEMaxIntvlEndTQueue.pop();
		}*/

		newInsert = (int)edges.size();
		tabuTChangePos = -1;
	}
	vec(int)* subCCs;//sub components
	bool* subCCsaved;//whether subCCs need to be saved  id -> saved/not saved
	int subCCNum;//the number of subCCs

	priority_queue<pair<int, NoisePos>> noisePosQueue; //the position of noise on each edge, pair=<timepos,<timeEnd,edgepos>> (maximum heap, key = timepos)
	int minEMaxIntvlEndT;
	priority_queue<QueueItem> preEMaxIntvlChangePos;//preMaxIntv[e] change when preMaxIntv[e][current-1].endT >= motifEndT

	priority_queue<pair<int, int>, vector<pair<int, int>>, less<pair<int, int>> > maxEMaxIntvlStartTQueue; //maxEMaxIntvlStartT on each edge, pair=<maxEMaxIntvlStartT, edgeid> (maximum heap, key = maxEMaxIntvlStartT)
	int tabuTChangePos;//the position that there exists a tabu time of an edge stop at t=nextNoiseChangePos+1
};


class CComponentsFRTMOPT1 : public CComponentsII {
public:
	CComponentsFRTMOPT1(int root) :CComponentsII(root) {
		haveNewEdges = 0;
		subCCs = nullptr;
		subCCsaved = nullptr;
	}
	CComponentsFRTMOPT1() :CComponentsII() { 
		haveNewEdges = 0;
		subCCs = nullptr;
		subCCsaved = nullptr;
	}
	virtual ~CComponentsFRTMOPT1() {
		if (subCCsaved != nullptr) delete[] subCCsaved;
		if (subCCs != nullptr) delete[] subCCs;
	}

	virtual void combineFrom(CComponentsFRTMOPT1*& other) {
		vec(int)& tempEdges = other->edges;//edges of currentCComponents
		int size = (int)edges.size();
		edges.resize(size + tempEdges.size());
		copy(tempEdges.begin(), tempEdges.end(), edges.begin() + size);

		while (!other->noisePosQueue.empty()) {//combine noise position
			auto noiseEdge = other->noisePosQueue.top();
			noisePosQueue.emplace(noiseEdge.first, NoisePos(noiseEdge.second.startTime, noiseEdge.second.pos + size));
			other->noisePosQueue.pop();
		}

		while (!other->preEMaxIntvlChangePos.empty()) {
			auto EMaxIntvlChangePos = other->preEMaxIntvlChangePos.top();
			preEMaxIntvlChangePos.emplace(EMaxIntvlChangePos.startTime, EMaxIntvlChangePos.endTime, EMaxIntvlChangePos.id);
			other->preEMaxIntvlChangePos.pop();
		}

		while (!other->maxEMaxIntvlStartTQueue.empty()) {
			auto item = other->maxEMaxIntvlStartTQueue.top();
			maxEMaxIntvlStartTQueue.emplace(item);
			other->maxEMaxIntvlStartTQueue.pop();
		}

		minEMaxIntvlEndT = min(other->minEMaxIntvlEndT, minEMaxIntvlEndT);
		
		newInsert = (int)edges.size();
		tabuTChangePos = -1;

		haveNewEdges |= other->haveNewEdges;

	}
	vec(int)* subCCs;//sub components
	bool* subCCsaved;//whether subCCs need to be saved  id -> saved/not saved
	int subCCNum;//the number of subCCs

	priority_queue<pair<int, NoisePos>> noisePosQueue; //the position of noise on each edge, pair=<timepos,<timeEnd,edgepos>> (maximum heap, key = timepos)
	priority_queue<QueueItem> preEMaxIntvlChangePos;//preMaxIntv[e] change when preMaxIntv[e][current-1].endT >= motifEndT

	priority_queue<pair<int, int>, vector<pair<int, int>>, less<pair<int, int>> > maxEMaxIntvlStartTQueue; //maxEMaxIntvlStartT on each edge, pair=<maxEMaxIntvlStartT, edgeid> (maximum heap, key = maxEMaxIntvlStartT)
	int minEMaxIntvlEndT;
	int tabuTChangePos;//the position that there exists a tabu time of an edge stop at t=nextNoiseChangePos+1
	bool haveNewEdges;//false: does not need check; true: need check
};

enum ComponentsD5Type : char {
	CCFRTM,
	CCPLUS
};

class CComponentsIID5Factory {
	
public:
	static CComponentsII* instance(int root, ComponentsD5Type type) {
		switch (type) {
			case CCFRTM:
				return DBG_NEW CComponentsFRTM(root);
			case CCPLUS:
				return DBG_NEW CComponentsFRTMOPT1(root);
			default:
				return nullptr;
		}
	}
};

#pragma endregion 

#pragma region temporal motifs

/*used in exact label matches*/
class TMotifI {
public:
	TMotifI() :startT(0), endT(0), motifEdge(nullptr), otherEdge(nullptr) {}

	TMotifI(vec(TEdge)& selectedEdge, int start, int end)
		: otherEdge(nullptr) {
		startT = start;
		endT = end;
		motifEdge = nullptr;
		//maskEdge = nullptr;
		setMotifEdge(selectedEdge);
	}

	TMotifI(int start, int end) : startT(start), endT(end), motifEdge(nullptr),
		otherEdge(nullptr) {}

	// copy-constructed
	TMotifI(const TMotifI& motif) {
		this->motifEdge = DBG_NEW vec(TEdge)(*motif.motifEdge);
		this->otherEdge = DBG_NEW vec(SaveCCInfo)(*motif.otherEdge);
		//this->maskEdge = DBG_NEW vec(int)(*motif.maskEdge);
		startT = motif.startT;
		endT = motif.endT;
		edgeNumber = motif.edgeNumber;
	}
	/* add edge to the motif*/
	inline void addEdge(int id, int edgeW) {
		motifEdge->emplace_back(id, edgeW);
		edgeNumber++;
	}

	/* add edge to the motif*/
	inline void addEdge(TEdge& edge) {
		motifEdge->emplace_back(edge);
	}

	inline void updateEdgeNumber() {
		edgeNumber = (int)motifEdge->size();
	}

	void sortEdges() {
		sort(motifEdge->begin(), motifEdge->end());
	}

	//size
	inline size_t getSize() {
		return edgeNumber;
	}

	vec(TEdge)* getMotifEdge() { return motifEdge; }

	//vec(int)* getMaskEdge() { return maskEdge; }

	inline int getStartT() const { return startT; }

	//inline int getMotifEdgeNumber() const { return edgeNumber; }

	inline int getEndT() const { return endT; }
	void setMotifEdge(vec(TEdge)& motifEdge) {
		if (this->motifEdge)
			delete this->motifEdge;
		this->motifEdge = DBG_NEW vec(TEdge)(motifEdge);
		this->edgeNumber = (int)this->motifEdge->size();
	}
	/* linked to other motifs*/
	/*inline int linkToMotifs(SaveCCInfo& saveInfo, TMotif*& motif) {
		if (this->otherEdge == nullptr)
			this->otherEdge = DBG_NEW vec(SaveCCInfo)();
		otherEdge->emplace_back(saveInfo);
		this->edgeNumber += motif->getSize();
		return otherEdge->size() - 1;
	}*/
	inline int linkToMotifs(SaveCCInfo& saveInfo, TMotifI*& motif) {
		if (this->otherEdge == nullptr)
			this->otherEdge = DBG_NEW vec(SaveCCInfo)();
		otherEdge->emplace_back(saveInfo);
		this->edgeNumber += (int)motif->getSize();
		return (int)otherEdge->size() - 1;
	}

	/* linked to other motifs for motifs with interval = [.,endT], in order to update SaveCCInfo of those motifs when they (interval = [.,endT]) update in the incremental algorithm*/
	inline void tempLinkToMotifs(SaveCCInfo& saveInfo) {
		if (this->otherEdge == nullptr)
			this->otherEdge = DBG_NEW vec(SaveCCInfo)();
		otherEdge->emplace_back(saveInfo);
	}

	vec(SaveCCInfo)* getOtherEdge() { return otherEdge; }

	~TMotifI() {
		if (otherEdge) {
			otherEdge->clear();
			delete otherEdge;
		}
	}

private:
	vec(TEdge)* motifEdge; //motif's edges with labels (new inserted)
	vec(SaveCCInfo)* otherEdge; //reuse other motif's edges (already inserted)
	int startT, endT; // starting timestamp(include) and ending timestamp(include)
	int edgeNumber;
	//for memory compression
};

/*used in RTMs*/
class TMotifII {
public:
	TMotifII() :startT(0), endT(0), motifEdge(nullptr)/*, maskEdge(nullptr) */{}

	TMotifII(vec(int)& selectedEdge, int start, int end) {
		startT = start;
		endT = end;
		motifEdge = nullptr;
		//maskEdge = nullptr;
		setMotifEdge(selectedEdge);
	}

	TMotifII(int start, int end) :startT(start), endT(end)
	{
		motifEdge = DBG_NEW vec(int)();
		//maskEdge = nullptr;
	}

	// copy-constructed
	TMotifII(const TMotifII& motif) {
		this->motifEdge = DBG_NEW vec(int)(*motif.motifEdge);
		//this->maskEdge = DBG_NEW vec(int)(*motif.maskEdge);
		startT = motif.startT;
		endT = motif.endT;
	}


	/* add edge to the motif*/
	/*inline void addEdge(int id, Label edgeW) {
		motifEdge->emplace_back(id, edgeW);
		edgeNumber++;
	}*/

	/* add edge to the motif*/
	inline void addEdge(int edgeId) {
		motifEdge->emplace_back(edgeId);
	}

	inline void copyEdges(vec(int)& motif) {
		motifEdge->resize(motif.size());
		copy(motif.begin(), motif.end(), motifEdge->begin());
	}

	/*inline void copyEdges(vec(int)& motif, pair<int,int>*tabel, int pos, int stime, int etime) {
		motifEdge->resize(motif.size());
		copy(motif.begin(), motif.end(), motifEdge->begin());

		for (auto eid : (*motifEdge)) {
			tabel[pos + eid].first = stime;
			tabel[pos + eid].second = etime;
		}
	}*/

	void sortEdges() {
		sort(motifEdge->begin(), motifEdge->end());
	}

	//size
	inline size_t getSize() {
		return motifEdge->size();
	}

	vec(int)* getMotifEdge() { return motifEdge; }

	//vec(int)* getMaskEdge() { return maskEdge; }

	inline int getStartT() const { return startT; }

	//inline int getMotifEdgeNumber() const { return edgeNumber; }

	inline int getEndT() const { return endT; }

	void setMotifEdge(vec(int)& motifEdge) {
		if (this->motifEdge)
			delete this->motifEdge;
		this->motifEdge = DBG_NEW vec(int)(motifEdge);
		//this->edgeNumber = this->motifEdge->size();
	}

	/*void setMaskEdge(vec(int)& maskEdge) {
		if (this->maskEdge)
			delete this->maskEdge;
		this->maskEdge = DBG_NEW vec(int)(maskEdge);
		this->edgeNumber -= this->maskEdge->size();
	}*/

	virtual ~TMotifII() {
		if (motifEdge) {
			motifEdge->clear();
			delete motifEdge;
		}
		/*if (maskEdge) {
			maskEdge->clear();
			delete maskEdge;
		}*/
	}

	inline void setInterval(int startT, int endT) {
		this->startT = startT;
		this->endT = endT;
	}

	inline void setEndT(int endT) {
		this->endT = endT;
	}

protected:
	vec(int)* motifEdge; //motif's edges with labels (new inserted)
	int startT, endT; // starting timestamp(include) and ending timestamp(include)
	//vec(int)* maskEdge;//edges which can not be saved 
};

#pragma endregion