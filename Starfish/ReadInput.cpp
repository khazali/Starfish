#include "MIfstream.h"
#include "Globals.h"
#include "Pore.h"
#include "Throat.h"
#include <cstdio>

#ifdef PETSC_ENABLED
#include <petscksp.h>
#include <mpi.h>
#endif

void ReadStatoilFormat(char *FilePath, char *Prefix, std::ofstream& DOut, std::ofstream& IOut) {
	MIfstream ThroatData1, ThroatData2, PoreData1, PoreData2, FluidData;
	char ThroatData1File[MAX_PATH_LENGTH], ThroatData2File[MAX_PATH_LENGTH], PoreData1File[MAX_PATH_LENGTH], PoreData2File[MAX_PATH_LENGTH], FluidDataFile[MAX_PATH_LENGTH], DOutFile[MAX_PATH_LENGTH], IOutFile[MAX_PATH_LENGTH];
	char str[MAX_STRING_LENGTH];
	unsigned int sLen, i;
#ifdef PETSC_ENABLED
	int MaxCN, j;
	int Cdata[2];
#endif



#ifdef PETSC_ENABLED	
	if (!rank) {
#endif

		sLen = (unsigned int)strlen(FilePath);
		if (FilePath[sLen - 1] != '\\') {		//MS Windows Path for now
			FilePath[sLen] = '\\';
			sLen++;
		}

		i = 0;
		while (Prefix[i]) {
			FilePath[sLen] = Prefix[i];
			sLen++;
			i++;
		}

		FilePath[sLen] = '\0';
		//sLen++;

		//std::cout<<FilePath;

		strcpy_s(ThroatData1File, MAX_PATH_LENGTH - 12, FilePath);
		strcpy_s(ThroatData2File, MAX_PATH_LENGTH - 12, FilePath);
		strcpy_s(PoreData1File, MAX_PATH_LENGTH - 12, FilePath);
		strcpy_s(PoreData2File, MAX_PATH_LENGTH - 12, FilePath);
		strcpy_s(FluidDataFile, MAX_PATH_LENGTH - 12, FilePath);
		strcpy_s(DOutFile, MAX_PATH_LENGTH - 12, FilePath);
		strcpy_s(IOutFile, MAX_PATH_LENGTH - 12, FilePath);

		ThroatData1File[sLen] = '_';
		ThroatData1File[sLen + 1] = 'l';
		ThroatData1File[sLen + 2] = 'i';
		ThroatData1File[sLen + 3] = 'n';
		ThroatData1File[sLen + 4] = 'k';
		ThroatData1File[sLen + 5] = '1';
		ThroatData1File[sLen + 6] = '.';
		ThroatData1File[sLen + 7] = 'd';
		ThroatData1File[sLen + 8] = 'a';
		ThroatData1File[sLen + 9] = 't';
		ThroatData1File[sLen + 10] = '\0';

		ThroatData2File[sLen] = '_';
		ThroatData2File[sLen + 1] = 'l';
		ThroatData2File[sLen + 2] = 'i';
		ThroatData2File[sLen + 3] = 'n';
		ThroatData2File[sLen + 4] = 'k';
		ThroatData2File[sLen + 5] = '2';
		ThroatData2File[sLen + 6] = '.';
		ThroatData2File[sLen + 7] = 'd';
		ThroatData2File[sLen + 8] = 'a';
		ThroatData2File[sLen + 9] = 't';
		ThroatData2File[sLen + 10] = '\0';

		PoreData1File[sLen] = '_';
		PoreData1File[sLen + 1] = 'n';
		PoreData1File[sLen + 2] = 'o';
		PoreData1File[sLen + 3] = 'd';
		PoreData1File[sLen + 4] = 'e';
		PoreData1File[sLen + 5] = '1';
		PoreData1File[sLen + 6] = '.';
		PoreData1File[sLen + 7] = 'd';
		PoreData1File[sLen + 8] = 'a';
		PoreData1File[sLen + 9] = 't';
		PoreData1File[sLen + 10] = '\0';

		PoreData2File[sLen] = '_';
		PoreData2File[sLen + 1] = 'n';
		PoreData2File[sLen + 2] = 'o';
		PoreData2File[sLen + 3] = 'd';
		PoreData2File[sLen + 4] = 'e';
		PoreData2File[sLen + 5] = '2';
		PoreData2File[sLen + 6] = '.';
		PoreData2File[sLen + 7] = 'd';
		PoreData2File[sLen + 8] = 'a';
		PoreData2File[sLen + 9] = 't';
		PoreData2File[sLen + 10] = '\0';

		FluidDataFile[sLen] = '_';
		FluidDataFile[sLen + 1] = 'f';
		FluidDataFile[sLen + 2] = 'l';
		FluidDataFile[sLen + 3] = 'u';
		FluidDataFile[sLen + 4] = 'i';
		FluidDataFile[sLen + 5] = 'd';
		FluidDataFile[sLen + 6] = '.';
		FluidDataFile[sLen + 7] = 'd';
		FluidDataFile[sLen + 8] = 'a';
		FluidDataFile[sLen + 9] = 't';
		FluidDataFile[sLen + 10] = '\0';

		DOutFile[sLen] = '_';
		DOutFile[sLen + 1] = 'D';
		DOutFile[sLen + 2] = 'r';
		DOutFile[sLen + 3] = 'a';
		DOutFile[sLen + 4] = 'i';
		DOutFile[sLen + 5] = 'n';
		DOutFile[sLen + 6] = 'a';
		DOutFile[sLen + 7] = 'g';
		DOutFile[sLen + 8] = 'e';
		DOutFile[sLen + 9] = '.';
		DOutFile[sLen + 10] = 'd';
		DOutFile[sLen + 11] = 'a';
		DOutFile[sLen + 12] = 't';
		DOutFile[sLen + 13] = '\0';

		IOutFile[sLen] = '_';
		IOutFile[sLen + 1] = 'I';
		IOutFile[sLen + 2] = 'm';
		IOutFile[sLen + 3] = 'b';
		IOutFile[sLen + 4] = 'i';
		IOutFile[sLen + 5] = 'b';
		IOutFile[sLen + 6] = 'i';
		IOutFile[sLen + 7] = 't';
		IOutFile[sLen + 8] = 'i';
		IOutFile[sLen + 9] = 'o';
		IOutFile[sLen + 10] = 'n';
		IOutFile[sLen + 11] = '.';
		IOutFile[sLen + 12] = 'd';
		IOutFile[sLen + 13] = 'a';
		IOutFile[sLen + 14] = 't';
		IOutFile[sLen + 15] = '\0';


		ThroatData1.open(ThroatData1File, std::fstream::in);
		ThroatData2.open(ThroatData2File, std::fstream::in);
		PoreData1.open(PoreData1File, std::fstream::in);
		PoreData2.open(PoreData2File, std::fstream::in);
		FluidData.open(FluidDataFile, std::fstream::in);

		DOut.open(DOutFile, std::fstream::out);
		IOut.open(IOutFile, std::fstream::out);


		TotalVolume = 0;
		FluidData.ReadWord(str);
		RhoW = atof(str);
		FluidData.ReadWord(str);
		RhoO = atof(str);
		DeltaRho = RhoW - RhoO;
		FluidData.ReadWord(str);
		IFT = atof(str);
		FluidData.ReadWord(str);
		RecedingContactAngle = atof(str);
		FluidData.ReadWord(str);
		AdvancingContactAngle = atof(str);
		FluidData.ReadWord(str);
		WaterViscosity = atof(str);
		FluidData.ReadWord(str);
		OilViscosity = atof(str);


		if (!ThroatData1.ReadWord(str)) TerM("Incorrect link1 file format!");
		ThroatNO = atoi(str);
		if (!PoreData1.ReadWord(str)) TerM("Incorrect node1 file format!");
		PoreNO = atoi(str);
		if (!PoreData1.ReadWord(str)) TerM("Incorrect node1 file format!");
		Dx = atof(str);
		if (!PoreData1.ReadWord(str)) TerM("Incorrect node1 file format!");
		Dy = atof(str);
		if (!PoreData1.ReadWord(str)) TerM("Incorrect node1 file format!");
		Dz = atof(str);

		throats = new Throat[ThroatNO];
		pores = new Pore[PoreNO];



		std::srand((unsigned int)time(NULL));
#ifdef PETSC_ENABLED
		MaxCN = 0;
#endif
		for (i = 0; i < PoreNO; i++) {
			pores[i].SetIndex(i);
			pores[i].ReadNode1(PoreData1);
			pores[i].ReadNode2(PoreData2);
#ifdef PETSC_ENABLED
			if (MaxCN < pores[i].GetCoordinationNumber()) MaxCN = pores[i].GetCoordinationNumber();
#endif
		}
		for (i = 0; i < ThroatNO; i++) {
			throats[i].SetIndex(i);
			throats[i].ReadLink1(ThroatData1);
			throats[i].ReadLink2(ThroatData2);
		}
		
#ifdef PETSC_ENABLED
		Cdata[0] = PoreNO;
		Cdata[1] = MaxCN;
	}
	MPI_Bcast(Cdata, 2, MPI_INT, 0, PETSC_COMM_WORLD);
	if (rank) PoreNO = Cdata[0];
	ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, Cdata[0], Cdata[0], 2 * Cdata[1], NULL, 2 * Cdata[1], NULL, &PetscA); CHKERRXX(ierr);
	ierr = MatZeroEntries(PetscA); CHKERRXX(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Cdata[0], &Petscb); CHKERRXX(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Cdata[0], &Petscx); CHKERRXX(ierr);
	ierr = MatGetOwnershipRange(PetscA, &Istart, &Iend); CHKERRXX(ierr);
	j = Iend - Istart;
	MPI_Allgather(&j, 1, MPI_INT, VectorMap, 1, MPI_INT, PETSC_COMM_WORLD);
	DisplacementMap[0] = 0;
	for (i = 0; i < wsize - 1; i++) DisplacementMap[i + 1] = DisplacementMap[i] + VectorMap[i];
#else
	TotalMSize = 0;
	Row = new int[PoreNO + 1];
#endif
#ifdef PETSC_ENABLED	
	if (!rank) {
#endif
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].CalcProps();
			pores[i].CalculateDrainagePc();
#ifndef PETSC_ENABLED	
			Row[i] = TotalMSize;
			pores[i].SetMatrixPlace(TotalMSize);
			TotalMSize += pores[i].NumberOfConnections() + 1;
#endif	
		}
#ifndef PETSC_ENABLED	
		Row[PoreNO] = TotalMSize;
#endif
#pragma omp parallel for
		for (i = 0; i < ThroatNO; i++) {
			throats[i].CalcProps();
			throats[i].CalculateDrainagePc();
		}

#ifndef PETSC_ENABLED
		Ans = new FloatType[PoreNO];
		Col = new int[TotalMSize];
		CoeffMatrix = new FloatType[TotalMSize];
#endif	

		ThroatData1.close();
		ThroatData2.close();
		PoreData1.close();
		PoreData2.close();
		FluidData.close();
#ifdef PETSC_ENABLED	
	}
#endif
}