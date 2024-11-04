#pragma once
#include "stdafx.h"
#include "Util.h"
#include "DynamicConnectivity.h"

/*disjoint set used for maintaining connected components
	(use path compression and union by rank)*/
class DisjointSet : public DynamicConnectivity {
private:
		int* parent;
		int* sizeOfT;//the size of subtree
		int size;
public:
	DisjointSet(int size)
		:size(size){
		parent = DBG_NEW int[size];
		sizeOfT = DBG_NEW int[size];
		reset();
	}
	~DisjointSet() {
		delete[]parent;
		delete[]sizeOfT;
	}

	/*find the root of the node whose position is num*/
	int findRoot(int num);
	
	/*union two trees where two nodes locate
	and return the old root of tree which is combined to other tree*/
	int addE(int e, int a, int b);

	inline int getSize()const { return size; }

	int* getParent()const { return parent; }

	int* getDepth()const { return sizeOfT; }

	void reset() {
		CLEARALL(parent, -1, size, int);
		for (int i = 0; i < size; i++) {
			//parent[i] = -1;
			sizeOfT[i] = 1;
		}
	}

	void resetAll(vec(int)*& edges, pair<int,int>*& edgeList, i2iHMap& vertex2Pos) {
		for (auto eid : (*edges)) {//O(Em)
			auto tempE = &edgeList[eid];
			int vertex = vertex2Pos[tempE->first];
			parent[vertex] = -1;
			sizeOfT[vertex] = 1;
			vertex = vertex2Pos[tempE->second];
			parent[vertex] = -1;
			sizeOfT[vertex] = 1;
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

	void removeE(int eid, int uid, int vid) {//only reset information of (uid,vid)
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
	inline void setValue(int pos, int root, int size) {
		parent[pos] = root;
		sizeOfT[pos] = size;
	}


	DisjointSet(const DisjointSet& ufset) {
		size = ufset.getSize();
		parent = DBG_NEW int[size];
		sizeOfT = DBG_NEW int[size];
		int* d = ufset.getDepth();
		int* p = ufset.getParent();
		for (int i = 0; i < size; ++i) {
			sizeOfT[i] = d[i];
			parent[i] = p[i];
		}
	}

	void print() {
		for (int i = 0; i < size; i++) {
			cout<<i<<"(p:"<<parent[i]<<
				",d:"<<sizeOfT[i]<<")"<<endl;
		}
	}

};