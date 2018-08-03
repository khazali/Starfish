#ifndef ELEMENTLIST_H
#define ELEMENTLIST_H
#include "Globals.h"

class ElementList {
private:
	unsigned int ListLength;
	unsigned int *List;
public:
	ElementList();
	~ElementList();
	void AddReadyToFillElement(unsigned int, ElementType);
	void SortElements(SortType);		//Batcher Sort
	bool SwapFlag(unsigned int, unsigned int, SortType);	
	void DestroyList(void);
	unsigned int GetListContent(unsigned int);
	unsigned int GetListLength(void);
	void AddElement(unsigned int);
	void RemoveElement(unsigned int);
};
#endif