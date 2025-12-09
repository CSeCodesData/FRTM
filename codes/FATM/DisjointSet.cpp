#include"DisjointSet.h"
#include"stdafx.h"

/*find the root of the node whose position is num*/
int DisjointSet::findRoot(int num) {
	if (num >= size) return -1;
	if (parent[num] == -1) return -1;
	int root = num;
	while (root != parent[root]) {
		root = parent[root];
		assert(root < size);
	}
	int temp=num,p;
	while (temp != root) {
		p = parent[temp];
		parent[temp] = root;
		temp = p;
	}
	
	return parent[num];
}
int DisjointSet::findRootDC(int num) {
	if (num >= size) return -1;
	if (parent[num] == -1) return -1;
	int root = num;
	while (root != parent[root]) {
		root = parent[root];
		if (root >= size) {
			break;
		}
	}
	int temp = num, p;
	while (temp != root) {
		p = parent[temp];
		parent[temp] = root;
		temp = p;
	}

	return parent[num];
}
/*union two trees where two nodes locate
	and return the old root of tree which is combined to other tree*/
int DisjointSet::addE(int a, int b) {
	
	if (parent[a] == -1 && parent[b] == -1) {
		parent[a] = parent[b] = b;
		sizeOfT[b] += 1;
	}
	else if (parent[a] == -1) {
		int p_b = findRoot(b);
		parent[a] = parent[p_b];
		sizeOfT[p_b] += 1;
	}
	else if (parent[b] == -1) {
		int p_a = findRoot(a); 
		parent[b] = parent[p_a];
		sizeOfT[p_a] += 1;
	}
	else {
		int p_a = findRoot(a);
		int p_b = findRoot(b);
		if (p_a == p_b)return -1;
		int rtn;
		if (sizeOfT[p_a] > sizeOfT[p_b]) {
			rtn = parent[p_b];
			parent[p_b] = parent[p_a]; 
			sizeOfT[p_a] += sizeOfT[p_b];
		}
		else {
			rtn = parent[p_a];
			parent[p_a] = parent[p_b]; 
			sizeOfT[p_b] += sizeOfT[p_a];
		}
		
		return rtn;
	}
	
	return -1;
}

int DisjointSetWithTime::addE(int a, int b, int time) {
	//not updated -> parent== -1 
	if ((updateT[a] < time || parent[a] == -1) && (updateT[b] < time || parent[b] == -1)) {
		parent[a] = parent[b] = b;
		sizeOfT[a] = 1;
		sizeOfT[b] = 2;
		updateT[a] = updateT[b] = time;
	}
	else if (updateT[a] < time || parent[a] == -1) {
		int p_b = findRoot(b);
		parent[a] = parent[p_b];
		sizeOfT[p_b] ++;
		sizeOfT[a] = 1;
		updateT[a] = time;
	}
	else if (updateT[b] < time || parent[b] == -1) {
		int p_a = findRoot(a);
		parent[b] = parent[p_a];
		sizeOfT[p_a] ++;
		sizeOfT[b] = 1;
		updateT[b] = time;
	}
	else {
		int p_a = findRoot(a);
		int p_b = findRoot(b);
		if (p_a == p_b)return -1;
		int rtn;
		if (sizeOfT[p_a] > sizeOfT[p_b]) {
			rtn = parent[p_b];
			parent[p_b] = parent[p_a];
			sizeOfT[p_a] += sizeOfT[p_b];
		}
		else {
			rtn = parent[p_a];
			parent[p_a] = parent[p_b];
			sizeOfT[p_b] += sizeOfT[p_a];
		}

		return rtn;
	}

	return -1;
}
