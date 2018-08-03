#include "ElementList.h"
#include "Pore.h"
#include "Throat.h"
#include <stdlib.h>

ElementList::ElementList(void) {
	ListLength=0;
	List=NULL;
}

void ElementList::AddReadyToFillElement(unsigned int ElementIndex, ElementType EType) {
	ListLength++;
	if ((List=(unsigned int*) realloc(List, ListLength*sizeof(int)))==NULL) TerM("Can not allocate memory for element list");
	switch (EType) {
		case PORE:
			List[ListLength-1]=ElementIndex;
			break;
		case THROAT:
			List[ListLength-1]=ElementIndex+PoreNO;
			break;
	}	
}

ElementList::~ElementList(void) {
	if (List!=NULL) {
		free(List);
		List=NULL;
	}
}

void ElementList::SortElements(SortType SType) {
	register unsigned int k, p, r, d, q, i, C2UMSize, u, TempVar;

	u=ListLength;
	k=0;
	while (u) {
		u>>=1;
		k++;
	}
	u=1;
	u<<=(k-1);

	if (u==ListLength) C2UMSize=ListLength;
	else C2UMSize=u<<1;

	i=C2UMSize>>1;
	for (p=i; p>=1; p>>=1) {
		r=0;
		d=p;
		for (q=i; q>=p; q>>=1) {
			for (k=0; k<(C2UMSize-d); k++) {
				if ((k & p)!=r) continue;				
				if ((k<ListLength) && ((k+d)<ListLength) && (SwapFlag(k, k+d, SType))) {					
					TempVar=List[k];
					List[k]=List[k+d];
					List[k+d]=TempVar;
				}
			}			
			d=q-p;
			r=p;
		}
	}
}

inline bool ElementList::SwapFlag(unsigned int First, unsigned int Second, SortType SType) {
	unsigned int EType1, EType2;
	FloatType Pc1, Pc2;

	EType1=List[First];
	EType2=List[Second];

	if (EType1<PoreNO) {											//Pore
		Pc1=pores[EType1].GetPc();
	}
	else if ((EType1>=PoreNO) && (EType1<(PoreNO+ThroatNO))) {		//Throat
		Pc1=throats[EType1-PoreNO].GetPc();
	}

	if (EType2<PoreNO) {											//Pore
		Pc2=pores[EType2].GetPc();
	}
	else if ((EType2>=PoreNO) && (EType2<(PoreNO+ThroatNO))) {		//Throat
		Pc2=throats[EType2-PoreNO].GetPc();
	}

	if (SType==ASCENDING) return (Pc1>Pc2);
	else return (Pc1<Pc2);
}

void ElementList::DestroyList(void) {
	if (List!=NULL) {
		free(List);
		List=NULL;
	}
	ListLength=0;
}

unsigned int ElementList::GetListContent(unsigned int ListIndex) {
	return List[ListIndex];
}

unsigned int ElementList::GetListLength(void){
	return ListLength;
}

void ElementList::AddElement(unsigned int ElementIndex) {
	ListLength++;
	if ((List=(unsigned int*) realloc(List, ListLength*sizeof(int)))==NULL) TerM("Can not allocate memory for element list");
	List[ListLength-1]=ElementIndex;
}

void ElementList::RemoveElement(unsigned int ElementIndex) {
	register unsigned int i, j;
	register bool cond = false;

	for (i = 0; i < ListLength; i++) if (List[i] == ElementIndex) {
		ListLength--;
		cond = true;
		break;
	}
	for (j = i; j < ListLength; j++) List[i] = List[i + 1];
	if (cond) if ((List = (unsigned int*) realloc(List, ListLength*sizeof(int))) == NULL) TerM("Can not allocate memory for element list");
}