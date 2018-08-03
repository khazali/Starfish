#include "Globals.h"
#include "Pore.h"
#include "Throat.h"
#include "ElementList.h"

void RecursiveSweepForDeletion(void) {
	register unsigned int i, ThroatIndex, ThroatsListLengh;

	for (i = 0; i < PoreNO; i++) {
		pores[i].SetIsVisited(false);
		pores[i].SetIsConnectedToOutlet(WATER, false);

	}
	for (i = 0; i < ThroatNO; i++) {
		throats[i].SetIsVisited(false);
		throats[i].SetIsConnectedToOutlet(WATER, false);
	}

	ThroatsListLengh = OutletThroats.GetListLength();
	for (i = 0; i < ThroatsListLengh; i++) {
		ThroatIndex = OutletThroats.GetListContent(i);
		if (!throats[ThroatIndex].GetIsVisited()) {
			throats[ThroatIndex].SetIsVisited(true);
			throats[ThroatIndex].SetIsConnectedToOutlet(WATER, true);
			throats[ThroatIndex].SweepAdjacentPoresForDeletion(OUTLET);
		}
	}

	
	for (i = 0; i < PoreNO; i++) {
		pores[i].SetIsVisited(false);
		pores[i].SetIsConnectedToInlet(false);

	}
	for (i = 0; i < ThroatNO; i++) {
		throats[i].SetIsVisited(false);
		throats[i].SetIsConnectedToInlet(false);
	}

	ThroatsListLengh = InletThroats.GetListLength();
	for (i = 0; i < ThroatsListLengh; i++) {
		ThroatIndex = InletThroats.GetListContent(i);
		if (!throats[ThroatIndex].GetIsVisited()) {
			throats[ThroatIndex].SetIsVisited(true);
			throats[ThroatIndex].SetIsConnectedToInlet(true);
			throats[ThroatIndex].SweepAdjacentPoresForDeletion(INLET);
		}
	}
}

void FilterIsolatedElements(void) {	
	register unsigned int i, j, n, Nulls;
	unsigned int PorePointerCounter, ThroatPointerCounter;

	RecursiveSweepForDeletion();

	PorePointer = new unsigned int[PoreNO];
	ThroatPointer = new unsigned int[ThroatNO];

	PorePointerCounter = 0;
	for (i = 0; i < PoreNO; i++) {
		PorePointer[i] = PorePointerCounter;
		if (pores[i].IsIsolated()) {
			PorePointerCounter++;
			//TotalVolume -= pores[i].GetReducedVolume();
			/*n = pores[i].GetCoordinationNumber();
			Nulls = 0;
			for (j = 0; j < n; j++) if (pores[i].IsConnectedThroatNull(j)) Nulls++;
			TotalMSize -= (n - Nulls + 1);*/
		}
	}
	
	ThroatPointerCounter = 0;
	for (i = 0; i < ThroatNO; i++) {
		ThroatPointer[i] = ThroatPointerCounter;
		if (throats[i].IsIsolated()) {
			ThroatPointerCounter++;
			//TotalVolume -= throats[i].GetReducedVolume();
			if ((throats[i].GetIOStat()) == 0) OutletThroats.RemoveElement(throats[i].GetIndex());
		}
	}

	PoreNO -= PorePointerCounter;
	ThroatNO -= ThroatPointerCounter;
	for (i = 0; i < PoreNO; i++) {
		pores[i] = pores[i + PorePointer[i]];
	}
	for (i = 0; i < ThroatNO; i++) {
		throats[i] = throats[i + ThroatPointer[i]];
	}

	

	InletThroats.DestroyList();
	delete[] PorePointer;
	delete[] ThroatPointer;	
}