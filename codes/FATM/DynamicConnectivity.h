#pragma once
class DynamicConnectivity {
public:
	DynamicConnectivity() = default;
	virtual ~DynamicConnectivity() {}
	virtual int findRoot(int uid) = 0; // find root of the tree
	virtual int addE(int eid, int uid, int vid) { return 0xffffffff; }
	virtual void addE(int uid, int vid) {}
	virtual void reset() = 0; // remove all information
	virtual void reset(int length) {} // remove part of information
};
