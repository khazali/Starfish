#include <fstream>

#include "Globals.h"
#include "ElementList.h"
#include "Pore.h"
#include "Throat.h"

void StoreToBinFile(void) {
	std::ofstream TempBin;
	register unsigned int i;
	bool tempBool;
	FloatType tempFloatType;

	TempBin.open("data\\TemBin",  std::ios::out|std::ios::binary);

	TempBin.write(reinterpret_cast<char *>(&TotalVolume), sizeof(TotalVolume));
	TempBin.write(reinterpret_cast<char *>(&Rmin), sizeof(Rmin));
	TempBin.write(reinterpret_cast<char *>(&TotalOilSaturation), sizeof(TotalOilSaturation));

	for (i=0; i<PoreNO; i++) {
		tempFloatType=pores[i].temp_GetOilSaturation();
		TempBin.write(reinterpret_cast<char *>(&tempFloatType), sizeof(tempFloatType));

		tempBool=pores[i].GetIsOilFilled();
		TempBin.write(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));

		tempBool=pores[i].GetIsWaterFilled();
		TempBin.write(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));
	}

	for (i=0; i<ThroatNO; i++) {
		tempFloatType=throats[i].temp_GetOilSaturation();
		TempBin.write(reinterpret_cast<char *>(&tempFloatType), sizeof(tempFloatType));

		tempBool=throats[i].GetIsOilFilled();
		TempBin.write(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));

		tempBool=throats[i].GetIsWaterFilled();
		TempBin.write(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));
	}

	TempBin.close();
}

void ReadFromBinFile(void) {
	std::ifstream TempBin;
	register unsigned int i;
	bool tempBool;
	FloatType tempFloatType;

	TempBin.open("data\\TemBin",  std::ios::in|std::ios::binary);

	TempBin.read(reinterpret_cast<char *>(&TotalVolume), sizeof(TotalVolume));
	TempBin.read(reinterpret_cast<char *>(&Rmin), sizeof(Rmin));
	TempBin.read(reinterpret_cast<char *>(&TotalOilSaturation), sizeof(TotalOilSaturation));

	for (i=0; i<PoreNO; i++) {
		TempBin.read(reinterpret_cast<char *>(&tempFloatType), sizeof(tempFloatType));
		pores[i].temp_SetOilSaturation(tempFloatType);

		TempBin.read(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));
		pores[i].temp_SetIsOilFilled(tempBool);

		TempBin.read(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));
		pores[i].SetIsWaterFilled(tempBool);
	}

	for (i=0; i<ThroatNO; i++) {
		TempBin.read(reinterpret_cast<char *>(&tempFloatType), sizeof(tempFloatType));
		throats[i].temp_SetOilSaturation(tempFloatType);

		TempBin.read(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));
		throats[i].temp_SetIsOilFilled(tempBool);

		TempBin.read(reinterpret_cast<char *>(&tempBool), sizeof(tempBool));
		throats[i].SetIsWaterFilled(tempBool);
	}

	TempBin.close();
}