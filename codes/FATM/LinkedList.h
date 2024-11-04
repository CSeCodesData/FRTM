#pragma once
#include "stdafx.h"
#include "Strategy.h"

/*one node of a link list*/
template<typename T>
struct LinkedNode {
public:
	LinkedNode(T vec)
		:item(vec),
		next(nullptr) {}
	
	~LinkedNode() {
		//item->clear();
		next = nullptr;
	}

	LinkedNode<T>* next;//next node
	T item;
};

/*one header of a link list*/
template<typename Y>
struct LinkedList {
public:
	LinkedList():/*length(0),*/first(nullptr),tail(nullptr) {}

	void addNodeAtFirst(LinkedNode<Y>*& node) {
		//if (length == 0)tail = node;
		if (first == nullptr) { 
			first = tail = node;
			return;
		}
		node->next = first;
		first = node;
		//length++;
	}

	void addNodeAtLast(LinkedNode<Y>*& node) {
		//if (length == 0)tail = node;
		if (tail == nullptr) first = node;
		else tail->next = node;
		tail = node;
		//length++;
	}

	void addItemAtFirst(Y item) {
		LinkedNode<Y>* node = DBG_NEW LinkedNode<Y>(item);
		//if (length == 0)tail = node;
		if (!first) {
			first = tail = node;
			return;
		}
		node->next = first;
		first = node;
		//length++;
	}

	void addItemAtLast(Y item) {
		LinkedNode<Y>* node = DBG_NEW LinkedNode<Y>(item);
		//if (length == 0)tail = node;
		if (tail == nullptr) first = node;
		else tail->next = node;
		tail = node;
		//length++;
	}

	//add the link list with header head(remove the head after adding)
	void addChain(LinkedList<Y>*& head) {
		if (head->tail == nullptr)return;
		//if (head->length == 0)return;
		if (tail == nullptr)tail = head->tail;
		//if (length == 0)tail = head->tail;
		head->tail->next = first;
		first = head->next;
		//length += head->length;
		head->next = nullptr;//for delete
	}
	//delete item
	void deleteItem(Y item) {
		LinkedNode<Y>* one = first;
		LinkedNode<Y>* pre = nullptr;
		while (one) {
			if (one->item == item) {
				if (!pre) {
					first = first->next;
					if (!first) tail = nullptr;
				}
				else {
					pre->next = one->next;
					if (!one->next) tail = pre;
				}
				delete one;
				return;
			}
			pre = one;
			one = one->next;
		}
	}
	//delete node after the current node
	void deleteNextNode(LinkedNode<Y>*& node) {
		LinkedNode<Y>* one = node->next;
		if (one == nullptr)return;
		node->next = one->next;
		if (one->next == nullptr)tail = node;
		delete one;
		//length -= 1;
	}

	//delete first node
	void deleteFirstNode() {
		LinkedNode<Y>* one = first;
		if (one == nullptr) return;
		first = one->next;
		if (first == nullptr) tail = nullptr;
		delete one;
		//length -= 1;
	}

	void release() {
		LinkedNode<Y>* one;
		while (first != nullptr) {
			one = first;
			first = first->next;
			delete one;
		}
		first = tail = nullptr;
		//length = 0;
	}

	~LinkedList() {
		release();
	}

	
	void output() {
		LinkedNode<Y>* now = first;
		while (now != nullptr) {
			cout << now->item << " "; 
			now = now->next;
		}
		cout << endl;
	}

	//int length;//length of linked list
	//Y key;// save item in first node
	LinkedNode<Y>* first;//first node
	LinkedNode<Y>* tail;//tail node
};


/*one node of a link list*/
template<typename T>
struct DoublyLinkedNode {
public:
	DoublyLinkedNode(T vec)
		:item(vec),
		next(nullptr),
		pre(nullptr) {}
	~DoublyLinkedNode() {
	}
	DoublyLinkedNode<T>* next;//next node
	DoublyLinkedNode<T>* pre;//pre node
	T item;
};

/*one header of a link list*/
template<typename Y>
struct DoublyLinkedList {
public:
	DoublyLinkedList() :/*length(0),*/first(nullptr), tail(nullptr) {}

	void addNodeAtFirst(DoublyLinkedNode<Y>*& node) {
		//if (length == 0)tail = node;
		if (first == nullptr) tail = node;
		else first->pre = node;
		node->next = first;
		first = node;
		//length++;
	}

	void addNodeAft(DoublyLinkedNode<Y>*& pos, Y&& item) {
		DoublyLinkedNode<Y>* node = DBG_NEW DoublyLinkedNode<Y>(item);
		node->next = pos->next;
		if (pos->next != nullptr) pos->next->pre = node;
		pos->next = node;
		node->pre = pos;
		if (tail == pos) tail = node;
	}

	void moveTo(DoublyLinkedNode<Y>*& pos, DoublyLinkedNode<Y>*& moved) {
		
		auto pre = moved->pre;
		if (pre != nullptr) pre->next = moved->next;
		else first = moved->next;
		if (moved->next != nullptr) moved->next->pre = pre;

		auto next = pos->next;
		if (next != nullptr) next->pre = moved;
		else tail = moved;
		moved->next = next;
		pos->next = moved;
		moved->pre = pos;
	}

	void addNodeAtLast(DoublyLinkedNode<Y>*& node) {
		//if (length == 0)tail = node;
		if (tail == nullptr) first = node;
		else tail->next = node;
		node->pre = tail;
		tail = node;
		//length++;
	}

	virtual void addItemAtFirst(Y&& item) {
		DoublyLinkedNode<Y>* node = DBG_NEW DoublyLinkedNode<Y>(item);
		//if (length == 0)tail = node;
		if (first == nullptr) tail = node;
		else first->pre = node;
		node->next = first;
		first = node;
		//length++;
	}

	virtual void addItemAtLast(Y&& item) {
		DoublyLinkedNode<Y>* node = DBG_NEW DoublyLinkedNode<Y>(item);
		//if (length == 0)tail = node;
		if (tail == nullptr) first = node;
		else tail->next = node;
		node->pre = tail;
		tail = node;
		//length++;
	}

	//add the link list with header head(remove the head after adding)
	void addChain(DoublyLinkedList<Y>*& head) {
		if (head->tail == nullptr)return;
		//if (head->length == 0)return;
		if (tail == nullptr)tail = head->tail;
		else //first!=nullptr
			first->pre = head->tail;
		//if (length == 0)tail = head->tail;
		head->tail->next = first;
		first = head->next;
		//length += head->length;
		head->next = nullptr;//for delete
	}

	//delete node after the current node
	void deleteNextNode(DoublyLinkedNode<Y>*& node) {
		DoublyLinkedNode<Y>* one = node->next;
		if (one == nullptr)return;
		node->next = one->next;
		if (one->next == nullptr)tail = node;
		else one->next->pre = node;
		delete one;
		//length -= 1;
	}

	//delete first node
	void deleteFirstNode() {
		DoublyLinkedNode<Y>* one = first;
		if (one == nullptr) return;
		first = one->next;
		if (first == nullptr) tail = nullptr;
		else first->pre = nullptr;
		delete one;
		//length -= 1;
	}

	void release() {
		DoublyLinkedNode<Y>* one;
		while (first != nullptr) {
			one = first;
			first = first->next;
			delete one;
		}
		first = tail = nullptr;
	}

	~DoublyLinkedList() {
		release();
	}

	void output() {
		DoublyLinkedNode<Y>* now = first;
		while (now != nullptr) {
			cout << now->item << " ";
			now = now->next;
		}
		cout << endl;
	}



	//int length;//length of linked list
	//Y key;// save item in first node
	DoublyLinkedNode<Y>* first;//first node
	DoublyLinkedNode<Y>* tail;//tail node
};




using ForbidTimeNode = LinkedNode<int>;
struct ForbidTimeList : public LinkedList<int> {
	void removeNodeLessThan(int filterNum) {
		LinkedNode<int>* now = first, *pre = nullptr;
		while (now != nullptr) {//O(delta)
			if (now->item < filterNum) {
				if (pre == nullptr) {
					deleteFirstNode();
					now = first;
				}
				else {
					deleteNextNode(pre);
					now = pre->next;
				}
			}
			else break;
		}
	}

	void removeNodeMoreThan(int filterNum) {
		LinkedNode<int>* now = first, *pre = nullptr;
		while (now != nullptr) {//O(delta)
			if (now->item <= filterNum) {
				pre = now;
				now = now->next;
			}
			else {
				if (pre == nullptr) {
					deleteFirstNode();
					now = first;
				}
				else {
					deleteNextNode(pre);
					now = pre->next;
				}
			}
		}
	}

	void copyToArr(int*& array) {
		LinkedNode<int>* now = first, *pre = nullptr;
		int pos = 0;
		while (now != nullptr) {//O(delta)
			array[pos++] = now->item;
			now = now->next;
		}
	}
};

/*used for the array vioT
*/
using Triad = tuple<int, int, int>; 
using ForbidTriadNode = DoublyLinkedNode<Triad>; 
struct ForbidTriadList : public DoublyLinkedList<Triad> {
	void removeNodeLessThan(int filterNum) {
		DoublyLinkedNode<Triad>* now = first, *temp = nullptr;
		int intvE, intvS;
		while (now != nullptr) {//O(delta)
			tie(intvS, intvE, std::ignore) = now->item;
			if (intvE < filterNum) {
				if (now->pre == nullptr) {
					deleteFirstNode();
					now = first;
				}
				else {
					temp = now->next;
					deleteNextNode(now->pre);
					now = temp;
				}
			}
			else if (intvS < filterNum) {
				get<0>(now->item) = filterNum;
				break;
			}
			else break;
		}
	}

	void removeNodeMoreThan(int filterNum) {
		DoublyLinkedNode<Triad>* now = tail, *temp = nullptr;
		int intvE, intvS; 
		while (now != nullptr) {//O(delta)
			tie(intvS, intvE, std::ignore) = now->item; 
			if (intvS > filterNum) {
				if (now->pre == nullptr) {
					deleteFirstNode();
					break;
				}
				else {
					temp = now->pre;
					deleteNextNode(temp);
					now = temp;
				}
			}
			else if (intvE > filterNum) {
				get<1>(now->item) = filterNum;
				break;
			}
			else break;
		}
	}
	void addItemAtLastWithCheck(Triad item) {
		DoublyLinkedNode<Triad>* node = DBG_NEW DoublyLinkedNode<Triad>(item);
		int firstTailItem, secondTailItem, thirdTailItem;
		int firstItem, secondItem, thirdItem;
		//if (length == 0)tail = node;
		if (tail == nullptr) { 
			first = node;
			tail = node;
		}
		else { 
			tie(firstTailItem, secondTailItem, thirdTailItem) = tail->item;
			tie(firstItem, secondItem, thirdItem) = item;
			if (thirdTailItem == -1 && thirdItem == -1) {
				if (secondTailItem >= firstItem && firstTailItem <= secondItem) {
					get<0>(tail->item) = min(firstTailItem, firstItem);
					get<1>(tail->item) = max(secondTailItem, secondItem);
				}
				else if (firstTailItem == secondItem + 1) {
					get<0>(tail->item) = firstItem;
				}
				else {
					tail->next = node;
					tail = node;
				}
			}
			else if (thirdTailItem != -1 && thirdItem != -1) {
				if (firstTailItem == firstItem) {
					if (thirdTailItem >= secondItem && secondTailItem <= thirdItem) {
						get<1>(tail->item) = min(secondTailItem, secondItem);
						get<2>(tail->item) = max(thirdTailItem, thirdItem);
					}
					else if (secondTailItem == thirdItem + 1) {
						get<1>(tail->item) = secondItem;
					}
					else {
						tail->next = node;
						tail = node;
					}
				}
				else {
					tail->next = node;
					tail = node;
				}
			}
			else {
				tail->next = node;
				tail = node;
			}
		}
		//length++;
	}
	void output() {
		ForbidTriadNode* now = first;
		while (now != nullptr) {
			cout << "("<< get<0>(now->item) << ", " << get<1>(now->item)<< ", " << get<2>(now->item) << ") ";
			now = now->next;
		}
		cout << endl;
	}

	void outputExceptLast() {
		ForbidTriadNode* now = first;
		while (now != nullptr&& now->next != nullptr) {
			cout << "(" << get<0>(now->item) << ", " << get<1>(now->item) << ", " << get<2>(now->item) << ") ";
			now = now->next;
		}
		cout << endl;
	}

	void output2() {
		ForbidTriadNode* now = first;
		int intvS, intvE;
		while (now != nullptr) {
			tie(intvS, intvE, std::ignore) = now->item;
			for (int start = intvS; start <= intvE; start++) {
				cout << start << " ";
			}
			now = now->next;
		}
		cout << endl;
	}


	//for testing
	bool check(vec(int) lis) {
		ForbidTriadNode* now = first;
		int intvS, intvE;
		if (lis.size() == 0 && now != nullptr) return false;
		for (auto iter = lis.begin(); iter != lis.end(); ++iter) {
			tie(intvS, intvE, std::ignore) = now->item;
			if (now != nullptr) {
				if (*iter < intvS || *iter > intvE) {
					return false;
				}
				else if (*iter == intvE) {
					now = now->next;
				}
			}
			else return false;
		}
		return true;
	}

};

using Pair = pair<int, int>;
using ForbidPairNode = DoublyLinkedNode<Pair>;
struct ForbidPairList : public DoublyLinkedList<Pair> {
	void copyTo(ForbidPairList& list) {
		DoublyLinkedNode<Pair>* temp = list.first;
		while (temp != nullptr) {
			this->addItemAtLast(make_pair(temp->item.first, temp->item.second));
			temp = temp->next;
		}
	}

	void removeNodeLessThan(int filterNum) {
		DoublyLinkedNode<Pair>* now = first, *temp = nullptr;
		//int intvE, intvS;
		while (now != nullptr) {//O(delta)
			//tie(intvS, intvE, std::ignore) = now->item;
			if (now->item.second < filterNum) {
				if (now->pre == nullptr) {
					assert(first != nullptr);
					deleteFirstNode();
					now = first;
				}
				else {
					temp = now->next;
					deleteNextNode(now->pre);
					now = temp;
				}
			}
			else if (now->item.first < filterNum) {
				get<0>(now->item) = filterNum;
				break;
			}
			else break;
		}
	}

	void removeNodeMoreThan(int filterNum) {
		DoublyLinkedNode<Pair>* now = tail, *temp = nullptr;
		//int intvE, intvS;
		while (now != nullptr) {//O(delta)
			//tie(intvS, intvE, std::ignore) = now->item;
			if (now->item.first > filterNum) {
				if (now->pre == nullptr) {
					deleteFirstNode();
					break;
				}
				else {
					temp = now->pre;
					deleteNextNode(temp);
					now = temp;
				}
			}
			else if (now->item.second > filterNum) {
				get<1>(now->item) = filterNum;
				break;
			}
			else break;
		}
	}
	void output() {
		ForbidPairNode* now = first;
		while (now != nullptr) {
			cout << "(" << now->item.first << ", " << now->item.second << ") ";
			now = now->next;
		}
		cout << endl;
	}

	void outputExceptLast() {
		ForbidPairNode* now = first;
		while (now != nullptr&& now->next != nullptr) {
			cout << "(" << now->item.first << ", " << now->item.second << ") ";
			now = now->next;
		}
		cout << endl;
	}

	void output2() {
		ForbidPairNode* now = first;
		int intvS, intvE;
		while (now != nullptr) {
			tie(intvS, intvE) = now->item;
			for (int start = intvS; start <= intvE; start++) {
				cout << start << " ";
			}
			now = now->next;
		}
		cout << endl;
	}


	//for testing
	bool check(vec(int) lis) {
		ForbidPairNode* now = first;
		int intvS, intvE;
		if (lis.size() == 0 && now != nullptr) return false;
		for (auto iter = lis.begin(); iter != lis.end(); ++iter) {
			tie(intvS, intvE) = now->item;
			if (now != nullptr) {
				if (*iter < intvS || *iter > intvE) {
					return false;
				}
				else if (*iter == intvE) {
					now = now->next;
				}
			}
			else return false;
		}
		return true;
	}

};

using EdgeIdArray = vector<int>;
using EdgeIdArrayNode = LinkedNode<EdgeIdArray>;
using EdgeIdArrayList = LinkedList<EdgeIdArray>;

