#pragma once
#include "stdafx.h"
#include "Util.h"
#include "Strategy.h"

template<typename T>
struct CircularQueue {
public:
	CircularQueue()
		:front(0),rear(0){
		q = DBG_NEW T[queueSize];
	}
	~CircularQueue() {
		delete[] q;
	}
	void reset() {
		CLEARALL(q,0, queueSize,T);
		front = 0;
		rear = 0;
	}
	void push(T&& item) {
		if ((rear + 1) % queueSize == front) {
			/*for (int temp = front; temp != rear;
				temp = (temp + 1) % queueSize) {
				cout << q[temp] << endl;
			}
			cout << Test::counter << " " << Test::counter2 << endl;*/
			Util::printError("queue full");
		}
		q[rear] = item;
		//rear = (rear + 1) % queueSize;
		rear++;
	}
	void pop() {
		assert(rear != front);
		//front = (front + 1) % queueSize;
		front++;
	}
	void pop_back() {
		assert(rear != front);
		//rear = (rear - 1 + queueSize) % queueSize;
		rear--;
	}
	inline T& top() {//need check empty
		return q[front];
	}
	inline T& last() {//need check empty
		return q[(rear+ queueSize -1)% queueSize];
	}
	inline bool empty() {
		return rear == front;
	}
	inline void swapToTop(int swapPos) {
		T temp = q[swapPos];
		q[swapPos] = q[front];
		q[front] = temp;
	}

	void copyTo(CircularQueue<T>& cp) {
		int temp = front;
		while (temp != rear) {
			cp.q[temp] = q[temp];
			//temp = (temp + 1) % queueSize;
			temp++;
		}
		cp.front = front;
		cp.rear = rear;
	}
	T* q;
	int front, rear;
	static int queueSize;
};