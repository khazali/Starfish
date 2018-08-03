#include <iostream>
#include <fstream>
#include "ElementList.h"
#include "Pore.h"
#include "Throat.h"
#include <iomanip>

void Imbibition(std::ofstream& OutFile) {
	ElementList EList;
	register unsigned int i, j;
	FloatType Pc;	

	//Start Water Flooding
#ifdef PETSC_ENABLED	
	if (!rank) {
#endif
		Rmin = IFT / CAPILLARYLIMIT;
		//std::cout << std::setprecision(20) << "Absolute Permeability = " << AbsPerm << std::endl << std::endl << std::endl;
		OutFile << std::setprecision(20) << "Absolute Permeability = " << AbsPerm << std::endl << std::endl << std::endl;
#ifdef PETSC_ENABLED	
	}
#endif
	Pc=CAPILLARYLIMIT;
	do {
#ifdef PETSC_ENABLED	
		if (!rank) {
#endif
			do {
				RecursiveSweepForConnection(OIL);
				RecursiveSweepForConnection(WATER);

				j = 0;
#pragma omp parallel for
				for (i = 0; i < PoreNO; i++) {
					j += pores[i].CalculateImbibitionPc(Pc);
				}
#pragma omp parallel for
				for (i = 0; i < ThroatNO; i++) {
					j += throats[i].CalculateImbibitionPc(Pc);
				}
			} while (j);


			if (Pc < 0) {
#pragma omp parallel for
				for (i = 0; i < PoreNO; i++) {
					pores[i].OilLayerExist(Pc);
				}
#pragma omp parallel for
				for (i = 0; i < ThroatNO; i++) {
					throats[i].OilLayerExist(Pc);
				}
			}

			TotalWaterSaturation = 0;
#pragma omp parallel for reduction(+:TotalWaterSaturation)
			for (i = 0; i < PoreNO; i++) {
				TotalWaterSaturation += pores[i].GetTotalWaterSaturation();
			}
#pragma omp parallel for reduction(+:TotalWaterSaturation)
			for (i = 0; i < ThroatNO; i++) {
				TotalWaterSaturation += throats[i].GetTotalWaterSaturation();
			}
#ifdef PETSC_ENABLED	
		}
#endif
		CalcRelPerm();
#ifdef PETSC_ENABLED	
		if (!rank) {
#endif
			OutFile << std::setprecision(20) << TotalWaterSaturation << '\t' << Pc << '\t' << WaterRelPerm << '\t' << OilRelPerm << std::endl;
			std::cout << std::setprecision(20) << TotalWaterSaturation << '\t' << Pc << '\t' << WaterRelPerm << '\t' << OilRelPerm << std::endl;
			OutFile.flush();
#ifdef PETSC_ENABLED	
		}
#endif

		Pc-=CAPILLARYINCREMENT;		
	//} while (Pc>=0);	
	} while (Pc>=(-CAPILLARYLIMIT));

	//TotalOilSaturation=1-TotalWaterSaturation;
}