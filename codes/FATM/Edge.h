#pragma once
#include"stdafx.h"

using NodePair = pair<uint32_t, uint32_t>;

class Edge {
public:
	Edge() :s(-1), t(-1), id(-1) {}

	Edge(uint32_t _s, uint32_t _t) :id(-1) {
		//if (_s < _t) {//ensure s < t
		s = _s;
		t = _t;
		//}
		//else {
		//	t = _s;
		//	s = _t;
		//}
	}

	Edge(const Edge& e)
		:s(e.s), t(e.t), id(e.id) {}

	bool operator<(const Edge& edge1) const {
		return s < edge1.s || (s == edge1.s && t < edge1.t);
	}

	bool operator!=(const Edge& edge1) const {
		return (s == edge1.s && t == edge1.t) ||
			(s == edge1.t&&t == edge1.s);
	}

	friend ostream& operator<<(ostream& output, const Edge& e) {
		//output << "edge{startid: " << e.s << " endid: " << e.t << " }";
		output << e.s << "," << e.t;
		return output;
	}

	bool operator==(const Edge& edge1) const {
		return (s == edge1.s && t == edge1.t) ||
			(s == edge1.t && t == edge1.s);
	}

	template <class H>
	friend H AbslHashValue(H h, const Edge& k) {
		return H::combine(move(h), k.s, k.t);
	}

	uint32_t s, t;//two endpoints, ensure s < t
	uint32_t id;//id of an edge
};


struct TEdge {//edge's id and label
	int id;
	int labelid;
	TEdge(int i, int lab) {
		id = i;
		labelid = lab;
	}
	TEdge() = default;
	bool operator==(const TEdge & tedge) const {
		return (id == tedge.id && labelid == tedge.labelid);
	}
	bool operator!=(const TEdge & tedge) const {
		return (id != tedge.id || labelid != tedge.labelid);
	}
	bool operator<(const TEdge& tedge) const {
		return id < tedge.id;
	}
	bool operator>(const TEdge& tedge) const {
		return id > tedge.id;
	}
};