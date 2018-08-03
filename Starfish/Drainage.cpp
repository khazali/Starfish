#include <iostream>
#include <fstream>
#include "ElementList.h"
#include "Pore.h"
#include "Throat.h"
#include <iomanip>


void Drainage(std::ofstream& OutFile) {
	ElementList EList;
	register unsigned int i, j;
	FloatType Pc;
	
	//Start Oil Flooding
	Pc=0;
	
	CalcAbsPerm();
	
#ifdef PETSC_ENABLED	
	if (!rank) {
#endif
		std::cout << std::endl << std::setprecision(20) << "Absolute Permeability = " << AbsPerm << std::endl << std::endl << std::endl;
		OutFile << std::setprecision(20) << "Absolute Permeability = " << AbsPerm << std::endl << std::endl << std::endl;
#ifdef PETSC_ENABLED	
	}
#endif

		do {
#ifdef PETSC_ENABLED	
			if (!rank) {
#endif
				do {
					RecursiveSweepForConnection(WATER);
					RecursiveSweepForConnection(OIL);

					j = 0;
#pragma omp parallel for
					for (i = 0; i < PoreNO; i++) {
						j += pores[i].TestFillReadiness(Pc);
					}
#pragma omp parallel for
					for (i = 0; i < ThroatNO; i++) {
						j += throats[i].TestFillReadiness(Pc);
					}

					TotalWaterSaturation = 0;
#pragma omp parallel for reduction(+:TotalWaterSaturation)
					for (i = 0; i < PoreNO; i++) {
						pores[i].CalculateOilSaturation(Pc);
						TotalWaterSaturation += pores[i].GetTotalWaterSaturation();
					}
#pragma omp parallel for reduction(+:TotalWaterSaturation)
					for (i = 0; i < ThroatNO; i++) {
						throats[i].CalculateOilSaturation(Pc);
						TotalWaterSaturation += throats[i].GetTotalWaterSaturation();
					}

				} while (j);

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

		Pc+=CAPILLARYINCREMENT;

	} while (Pc<=CAPILLARYLIMIT);

#ifdef PETSC_ENABLED	
	if (!rank) {
#endif	
		TotalOilSaturation = 1 - TotalWaterSaturation;
#ifdef PETSC_ENABLED	
	}
#endif
}