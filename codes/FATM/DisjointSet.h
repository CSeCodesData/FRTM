#pragma once
#include "stdafx.h"
#include "Util.h"
#include <absl/container/flat_hash_map.h>

/*disjoint set used for maintaining connected components
	(use path compression and union by rank)*/
class DisjointSet {
protected:
		int* parent;
		int* sizeOfT;//the size of subtree
		int size;
public:
	DisjointSet(){}
	DisjointSet(int size)
		:size(size){
		parent = DBG_NEW int[size];
		sizeOfT = DBG_NEW int[size];
		reset();
	}
	virtual ~DisjointSet() {
		delete[]parent;
		delete[]sizeOfT;
	}

	/*find the root of the node whose position is num*/
	virtual int findRoot(int num);
	virtual int findRootDC(int num);

	/*union two trees where two nodes locate
	and return the old root of tree which is combined to other tree*/
	virtual int addE(int a, int b);

	virtual inline int getMaxRound() const { return -1; }

	inline int getSize()const { return size; }
	
	inline int getSizeOfTree(int num) const { return sizeOfT[num]; }

	inline int getSizeOfTreeForRoot(int num) const { 
		int root = num, child;
		while (root != parent[root]) {
			child = root;
			root = parent[root];
			if (root > size) {
				return sizeOfT[child];
			}
		}
		return sizeOfT[root];
	}

	virtual void reset() {
		CLEARALL(parent, -1, size, int);
		for (int i = 0; i < size; i++) {
			//parent[i] = -1;
			sizeOfT[i] = 1;
		}
	}

	void checkValid() {
		for (int i = 0; i < size; i++) {
			if (parent[i]< -1 || parent[i] >= size) {
				cout << "#" << endl;
				exit(0);
			}
		}
	}

	void removeE(int uid, int vid) {//only reset information of (uid,vid)
		parent[uid] = -1;
		parent[vid] = -1;
		sizeOfT[uid] = 1;
		sizeOfT[vid] = 1;
	}

	void reset(int resetSize) {
		CLEARALL(parent, -1, resetSize, int);
		for (int i = 0; i < resetSize; i++) {
			//parent[i] = -1;
			sizeOfT[i] = 1;
		}
	}
	inline void resetPos(int pos) {
		parent[pos] = -1;
		sizeOfT[pos] = 1;
	}
	inline void resetPos(int pos1, int pos2) {
		parent[pos1] = -1;
		sizeOfT[pos1] = 1;
		parent[pos2] = -1;
		sizeOfT[pos2] = 1;
	}
	inline int getSizeValue(int pos) {
		return  sizeOfT[pos];
	}
	inline void setValue(int pos1, int pos2, int root, int size1, int size2) {
		parent[pos1] = root;
		sizeOfT[pos1] = size1;
		parent[pos2] = root;
		sizeOfT[pos2] = size2;
	}

	DisjointSet(const DisjointSet& ufset) {
		size = ufset.getSize();
		parent = DBG_NEW int[size];
		sizeOfT = DBG_NEW int[size];
		copy(ufset.sizeOfT, ufset.sizeOfT + size, sizeOfT);
		copy(ufset.parent, ufset.parent + size, parent);
	}

	virtual void print() {
		for (int i = 0; i < size; i++) {
			cout<<i<<"(p:"<<parent[i]<<
				",d:"<<sizeOfT[i]<<")"<<endl;
		}
	}
	
};


class DisjointSetWithTime : public DisjointSet {
public:
	int* updateT;//update time
	DisjointSetWithTime(int size)
		:DisjointSet(size) {
		updateT = DBG_NEW int[size];
		CLEARALL(updateT, -1, size, int);
	}
	virtual ~DisjointSetWithTime() {
		delete[]updateT;
	}
	virtual void print(int pos) {
		cout << "p:" << parent[pos] <<
			",s:" << sizeOfT[pos] <<",t:"<<updateT[pos]<< endl;
	}
	/*union two trees where two nodes locate
	and return the old root of tree which is combined to other tree*/
	int addE(int a, int b, int time);

	inline void setParent(int n, int p, int time) {
		parent[n] = p;
		updateT[n] = time;
	}
	inline void setParent(int n, int p) {
		parent[n] = p;
	}
	inline int getParent(int pos, int time) {
		if (updateT[pos] < time) return -1;
		return parent[pos];
	}
};
