#include "Globals.h"
#include "Pore.h"
#include "Throat.h"
#include "ElementList.h"

void RecursiveSweepForConnection(FluidType Fluid) {
	register unsigned int i, OutletThroatIndex, OutletThroatsListLengh;

#pragma omp parallel for
	for (i=0; i<PoreNO; i++) {
		pores[i].SetIsVisited(false);
		pores[i].SetIsConnectedToOutlet(Fluid, false);

	}
#pragma omp parallel for
	for (i=0; i<ThroatNO; i++) {
		throats[i].SetIsVisited(false);
		throats[i].SetIsConnectedToOutlet(Fluid, false);
	}

	OutletThroatsListLengh=OutletThroats.GetListLength();
#pragma omp parallel for
	for (i=0; i<OutletThroatsListLengh; i++) {
		OutletThroatIndex=OutletThroats.GetListContent(i);
		if ((!throats[OutletThroatIndex].GetIsVisited()) && (((Fluid==OIL) && (throats[OutletThroatIndex].GetIsOilFilled())) || ((Fluid==WATER) && (throats[OutletThroatIndex].GetWaterCoatingExist())))) {
			throats[OutletThroatIndex].SetIsVisited(true);
			throats[OutletThroatIndex].SetIsConnectedToOutlet(Fluid, true);
			throats[OutletThroatIndex].SweepAdjacentPores(Fluid);
		}
	}
}